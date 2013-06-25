"""
A NetCDF aggregator.

Suppose we have several model outputs in files output1.nc, output2.nc,
etc., and we want to combine them into a single dataset by creating
a new axis called "ensemble". This would be done by creating the
following file with the extension ``.nca``::

    [dataset]
    ; obligatory parameters
    match = output*.nc
    axis = ensemble
    ; optional metadata
    name = MOM4_ensemble
    history = Ensemble members from a simulation using MOM4
    
    [ensemble]
    values = 1, 2, ...
    ; optional metadata
    units = member

Values can be specified using 3 different syntaxes::

    values = 1, 2, ...
    values = 1.0, ..., 10.0
    values = 1, 2, 3, 4, 5, 6

The first example assigns a value of 1 for the first point, 2 for
the second, and so on, with a delta of 1 up to the last point. The
second example assigns a value of 1.0 for the first point, 10.0 for
the last, and evenly distributed values in between, depending on
the number of files. The last one explicitly assigns values to each
point.

Another example: we want to combine monthly output in a single
dataset. This is done by aggregating along an existing axis (say,
"time")::

    [dataset]
    name = annual
    match = output2009??.nc
    axis = time

The aggregation is done lazily, ie, it avoids having to read data
from all files into memory in order to build the new dataset. Data
is only read as necessary.

"""

import re
import os.path
from glob import glob
import string

from numpy import asarray, swapaxes, linspace, concatenate, arange
from configobj import ConfigObj
from pupynere import netcdf_file

from pydap.model import *
from pydap.lib import combine_slices, fix_slice
from pydap.handlers.lib import BaseHandler
from pydap.handlers.helper import constrain
from pydap.responses.das import get_type
from pydap.util.safeeval import expr_eval


def lazy_eval(s):
    """
    Try to evalute expression or fallback to string.
    
        >>> lazy_eval("1")
        1
        >>> lazy_eval("a")
        'a'

    """
    try:
        s = expr_eval(s)
    except:
        pass
    return s


class Handler(BaseHandler):

    extensions = re.compile(r"^.*\.nca$", re.IGNORECASE)

    def __init__(self, filepath):
        config = ConfigObj(filepath)
        files = glob(config['dataset']['match'])
        template = netcdf_file(files[0])  # use first file as a template
        to_close = [ template ]           # handlers we need to close later

        # Create dataset and populate with metadata from the INI file.
        self.dataset = DatasetType(
                name=config['dataset'].get('name', os.path.split(filepath)[1]),
                attributes=dict((k, v) for (k, v) in config['dataset'].items()
                    if k not in ['name', 'match', 'axis']))
        self.dataset.attributes['NC_GLOBAL'] = template._attributes.copy()

        # Build aggregation axis.
        name = config['dataset']['axis']
        if name in config:
            # Create a new axis.
            data = parse_values(config[name]['values'], len(files))
            type_ = get_type(data)
            attributes = dict((k, v) for (k, v) in config[name].items()
                    if k != 'values')
        elif name in template.variables:
            # Combine existing axes into a single one.
            data = []
            for file in files:
                f = netcdf_file(file)
                to_close.append(f)
                data.extend(f.variables[name][:])
            data = asarray(data)
            type_ = template.variables[name].typecode()
            attributes = template.variables[name]._attributes.copy()
        else:
            raise Exception('Invalid axis definition "%s".' % name)
        agg_axis = BaseType(name=name, data=data, type=type_, shape=data.shape, attributes=attributes)

        # Add grids to the dataset.
        for name in [name for name in template.variables
                if name not in template.dimensions]:
            var = template.variables[name]
            self.dataset[name] = grid = GridType(name=name, attributes=var._attributes.copy())
            grid[name] = BaseType(name=name, type=type_, attributes=var._attributes.copy())

            # Prepare aggregated data.
            if agg_axis.name not in var.dimensions:
                axis = None
                self.dataset[agg_axis.name] = grid[agg_axis.name] = agg_axis
                grid[name].dimensions = (agg_axis.name,) + var.dimensions
            else:
                axis = list(var.dimensions).index(agg_axis.name)
                grid[name].dimensions = var.dimensions
            fs = [netcdf_file(file) for file in files]
            to_close.extend(fs)
            grid[name].data = LazyConcatenator(
                    [f.variables[name] for f in fs],
                    axis=axis)
            grid[name].shape = grid[name].data.shape

            # Add axes for the dimensions.
            for dim in var.dimensions:
                if dim == agg_axis.name:
                    # replace axis with aggregated axis
                    self.dataset[dim] = grid[dim] = agg_axis
                else:
                    # add a regular axis
                    self.dataset[dim] = grid[dim] = BaseType(name=dim,
                            data=template.variables[dim][:],
                            type=template.variables[dim].typecode(),
                            shape=template.variables[dim].shape,
                            attributes=template.variables[dim]._attributes.copy())

        # Close all open files.
        self.dataset.close = lambda: [f.close() for f in to_close]

    def parse_constraints(self, environ):
        # Let's use the automatic constraint evaluator.
        new_dataset = constrain(self.dataset, environ.get('QUERY_STRING', ''))
        if hasattr(self.dataset, 'close'): self.dataset.close()
        return new_dataset


class LazyConcatenator(object):
    """
    A lazy concatenator.

    This class is similar to the ``concatenate`` function from Numpy,
    except that the new object is a "view" of the concatenated arrays,
    and data is only read when this object is sliced.

    """
    def __init__(self, arrays, axis=None):
        if axis is None:
            # Create a new axis
            self.shape = (len(arrays),) + arrays[0].shape
            self._arrays = asarray(arrays, 'O') 
        else:
            # 'cat along a given axis
            shape = list(arrays[0].shape)
            shape[axis] = sum([array.shape[axis] for array in arrays])
            self.shape = tuple(shape)
            self._arrays = asarray(
                    [LazySlice((slice(None),)*axis+(i,), arrays[j]) 
                        for j in range(len(arrays))
                        for i in range(arrays[j].shape[axis])],
                    'O')
        self._axis = axis

    def __getitem__(self, index):
        # Extract the slice for the object array
        index = fix_slice(index, self.shape)
        if self._axis is None:
            select, index = index[0], index[1:]
            out = asarray(
                    [obj[index] for obj in self._arrays[select]])
        else:
            select = index[self._axis]
            out = concatenate(
                    [obj[index] for obj in self._arrays[select]],
                    axis=self._axis)
        return out

    @property
    def __array_interface__(self):
        data = self[:]
        return {
                'version': 3,
                'shape': data.shape,
                'typestr': data.dtype.str,
                'data': data,
                }


class LazySlice(object):
    """
    A lazy slice.

    This is a lazy slice. The initial slice is applied only when
    the object is sliced again; both slices are then combined
    and applied to the data.

    """
    def __init__(self, slice, obj):
        self.slice = fix_slice(slice, obj.shape)
        self.obj = obj

    def __getitem__(self, index):
        index = fix_slice(index, self.obj.shape)
        return self.obj[combine_slices(index, self.slice)]


def parse_values(input, size):
    """
    Parse fuzzy values for the aggregation axis. The input comes from
    ConfigObj as a list of strings::

        >>> print parse_values(["10", "20", "..."], 5)
        [ 10.  20.  30.  40.  50.]
        >>> print parse_values(["1", "...", "10"], 5)
        [  1.     3.25   5.5    7.75  10.  ]
        >>> print parse_values(["1", "1", "2", "3", "5"], 5)
        [1 1 2 3 5]

    """
    if len(input) == size and "..." not in input:
        return asarray(map(lazy_eval, input))
    start, next, stop = input[:]
    if next == '...':
        return linspace(float(start), float(stop), size)
    elif stop == '...':
        dx = float(next)-float(start)
        return arange(float(start), float(start)+size*dx, dx)
    else:
        raise Exception('Unable to parse: %s' % input)


if __name__ == '__main__':
    import sys
    from paste.httpserver import serve

    application = Handler(sys.argv[1])
    serve(application, port=8001)
