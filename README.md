pydap.handlers.nca
==================

A Pydap handler that allows aggregating multiple NetCDF into a single dataset. 
The configuration is extremely simple. As an example, to aggregate model output
in different files (say, `output1.nc`, `output2.nc`, etc.) along a new axis
"ensemble" just create an `INI` file with the extension `.nca`:

```ini
[dataset]
match = /path/to/output*.nc
axis = ensemble
; below optional metadata:
history = Test for NetCDF aggregator

[ensemble]
values = 1, 2, ...
long_name = Ensemble members
```

This will assign the values 1, 2, and so on to each ensemble member. Another 
example, suppose we have monthly data in files `data01.nc`, `data02.nc`, ...,
`data12.nc`, and we want to aggregate along the `TIME` axis:

```ini
[dataset]
match = /path/to/data*.nc
axis = TIME  # existing axis
```

The handler currently works with NetCDF 3 files, but in the future it will be
changed to work with any other Pydap-supported data format. 
