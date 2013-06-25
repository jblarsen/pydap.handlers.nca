"""A Pydap handler for aggregating NetCDF files."""

import os
from glob import glob
import re
import time
from stat import ST_MTIME
from email.utils import formatdate
from ast import literal_eval as eval

import numpy as np
from configobj import ConfigObj

from pupynere import netcdf_file

from pydap.model import *
from pydap.handlers.lib import BaseHandler
from pydap.exceptions import OpenFileError
from pydap.lib import fix_slice


class FileAttributeCache(object):

    """Cache for storing and retrieving values tied to files which can change.

    This is used so we don't need to open NetCDF files for each request to get
    data from the aggregation axis.

    """

    def __init__(self):
        self.cache = {}
        self.mtime = {}

    def __setitem__(self, filepath, value):
        """Assign a value to the file `filepath`."""
        self.cache[filepath] = value
        self.mtime[filepath] = os.stat(filepath)[ST_MTIME]

    def __getitem__(self, filepath):
        """Try to retrieve the value associated with `filepath`.

        Return the cached valud or raise and exception if the file has been
        modified.

        """
        if filepath not in self.cache:
            raise KeyError

        mtime = os.stat(filepath)[ST_MTIME]
        if mtime > self.mtime[filepath]:
            del self.cache[filepath]
            raise KeyError

        return self.cache[filepath]


class PseudoNetCDFVariable(object):

    """Create a pseudo NetCDF variable.

    This is used to create the aggregation axis as if it were a variable in the
    original NetCDF files. Basically we need just to return data and
    attributes.

    """

    def __init__(self, data, attributes):
        self.data = data
        self._attributes = attributes

    def __getitem__(self, index):
        return self.data[index]


class NetCDFAggregatorHandler(BaseHandler):

    """A NetCDF 3 aggregator."""

    extensions = re.compile(r"^.*\.nca$", re.IGNORECASE)

    def __init__(self, filepath):
        BaseHandler.__init__(self)

        self.filepath = filepath
        self.cache = FileAttributeCache()

    def parse(self, projection, selection):
        """Build aggregrated dataset.

        Aggregation must be done per request, since new files may be added.
        Once the dataset is built, we return a call to the superclass.

        """
        try:
            config = ConfigObj(self.filepath)
        except Exception, exc:
            message = 'Unable to open file %s: %s' % (self.filepath, exc)
            raise OpenFileError(message)

        # get list of files
        files = glob(config['dataset']['match'])
        files.sort()

        # add last-modified header, taking the last modified file from the list
        mtime = max(os.stat(file)[ST_MTIME] for file in files)
        self.additional_headers = [
            (k, v) for (k, v) in self.additional_headers
            if k != 'Last-modified']
        self.additional_headers.append(
            ('Last-modified', (
                formatdate(time.mktime(time.localtime(mtime))))))

        # get a template
        template = netcdf_file(files[0])
        vars = template.variables
        dims = template.dimensions

        # build aggregation axis
        axis = config['dataset']['axis']
        if axis in config:
            # create a new axis
            values = config[axis]['values']
            data = parse_values(values, len(files))
            attributes = {
                k: v for k, v in config[axis].items() if k != 'values'
            }
            count = [(file, 0) for file in files]
        elif axis in dims:
            # aggregate along existing axis
            count = []
            data = []
            for file in files:
                try:
                    values = self.cache[file]
                except KeyError:
                    f = netcdf_file(file)
                    values = f.variables[axis][:]
                    f.close()
                    self.cache[file] = values
                data.extend(values)
                count.append((file, len(values)))
            data = np.array(data)
            attributes = vars[axis]._attributes
        else:
            raise Exception('Invalid axis definition "%s".' % name)

        # add new axis, or replace if existing
        dims[axis] = len(data)
        vars[axis] = PseudoNetCDFVariable(data, attributes)

        # build dataset
        name = os.path.split(self.filepath)[1]
        self.dataset = DatasetType(
            name, attributes=dict(NC_GLOBAL=template._attributes))

        # add grids
        grids = [var for var in vars if var not in dims]
        for grid in grids:
            # should we add the axis to the dimensions?
            if axis not in vars[grid].dimensions:
                vars[grid].dimensions = (axis,) + vars[grid].dimensions
                shape = (dims[axis],) + vars[grid].shape
                index = None
            else:
                shape = list(vars[grid].shape)
                index = vars[grid].dimensions.index(axis)
                shape[index] = dims[axis]
                shape = tuple(shape)
            dtype = np.dtype(vars[grid].typecode())

            # create grid
            self.dataset[grid] = GridType(grid, vars[grid]._attributes)

            # add array
            self.dataset[grid][grid] = BaseType(
                grid, AggregateNetcdfData(grid, index, count, dtype, shape),
                vars[grid].dimensions, vars[grid]._attributes)

            # add maps
            for dim in vars[grid].dimensions:
                self.dataset[grid][dim] = BaseType(
                    dim, vars[dim][:], None, vars[dim]._attributes)

        # add dims
        for dim in dims:
            self.dataset[dim] = BaseType(
                dim, vars[dim][:], None, vars[dim]._attributes)

        template.close()

        return BaseHandler.parse(self, projection, selection)


class AggregateNetcdfData(object):

    """Aggregates variables from NetCDF files along a given axis."""

    def __init__(self, name, axis, count, dtype, shape):
        self.name = name
        self.axis = axis
        self.count = count
        self.dtype = dtype
        self.shape = shape

    def __getitem__(self, index):
        index = fix_slice(index, self.shape)

        # create a new axis
        if self.axis is None:
            # get the slice along the aggregation axis, and leave the rest for
            # the variable itself
            slice_, index = index[0], index[1:]
            data = []
            for file, n in self.count[slice_]:
                f = netcdf_file(file)
                data.append(f.variables[self.name][index])
                f.close()
            return np.array(data).astype(self.dtype)

        # concatenate along an existing axis
        else:
            # convert index to list so we can change it
            index = list(index)

            # get the slice along the aggregation axis and store it in a
            # boolean array that we'll map to the files
            slice_ = index[self.axis]
            indexes = np.zeros(self.shape[self.axis], bool)
            indexes[slice_] = 1

            offset = 0
            data = []
            for file, n in self.count:
                selected_here = indexes[offset:offset+n]
                if any(selected_here):
                    index[self.axis] = selected_here
                    f = netcdf_file(file)
                    data.append(f.variables[self.name][tuple(index)])
                    f.close()
                offset += n
            return np.concatenate(data, axis=self.axis).astype(self.dtype)

    # Comparisons are passed to the data.
    def __eq__(self, other):
        return self[:] == other

    def __ne__(self, other):
        return self[:] != other

    def __ge__(self, other):
        return self[:] >= other

    def __le__(self, other):
        return self[:] <= other

    def __gt__(self, other):
        return self[:] > other

    def __lt__(self, other):
        return self[:] < other

    # Implement the sequence and iter protocols.
    def __len__(self):
        return self.shape[0]

    def __iter__(self):
        return iter(self[:])


def parse_values(input, size):
    """Parse fuzzy values for the aggregation axis.

    The input comes from ConfigObj as a list of strings::

        >>> print parse_values(["10", "20", "..."], 5)
        [ 10.  20.  30.  40.  50.]
        >>> print parse_values(["1", "...", "10"], 5)
        [  1.     3.25   5.5    7.75  10.  ]
        >>> print parse_values(["1", "1", "2", "3", "5"], 5)
        [1 1 2 3 5]

    Returns a list with the values.

    """
    if len(input) == size and "..." not in input:
        return np.asarray(map(eval, input))
    start, next, stop = input[:]
    if next == '...':
        return np.linspace(eval(start), eval(stop), size)
    elif stop == '...':
        dx = eval(next)-eval(start)
        return np.arange(eval(start), eval(start)+size*dx, dx)
    else:
        raise Exception('Unable to parse: %s' % input)


if __name__ == "__main__":
    import sys
    from werkzeug.serving import run_simple

    application = NetCDFAggregatorHandler(sys.argv[1])
    run_simple('localhost', 8001, application, use_reloader=True)
