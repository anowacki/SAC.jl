# SAC.jl

[![Build Status](https://img.shields.io/travis/anowacki/SAC.jl.svg?style=flat-square&label=linux)](https://travis-ci.org/anowacki/SAC.jl)
[![Build status](https://img.shields.io/appveyor/ci/AndyNowacki/sac-jl.svg?style=flat-square&label=windows)](https://ci.appveyor.com/project/AndyNowacki/sac-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/anowacki/SAC.jl/badge.svg?branch=master)](https://coveralls.io/github/anowacki/SAC.jl?branch=master)

## What is SAC.jl?
A [Julia](http://julialang.org) package for dealing with seismic data in the
[SAC](http://ds.iris.edu/files/sac-manual/manual/file_format.html) format, and
processing that data in a similar way to the SAC program:
either the [SAC/BRIS](http://www1.gly.bris.ac.uk/~george/sac-bugs.html) or
[IRIS](http://ds.iris.edu/ds/nodes/dmc/software/downloads/sac/) versions.


## How to install
Although not registered as an official package, SAC.jl can be added to your
Julia install like so:

```julia
Pkg.clone("https://github.com/anowacki/SAC.jl")
```

You then need only do

```julia
using SAC
```

and if that works, you're ready to go.


## How to use
### SACtr type
SAC.jl represents SAC files with the `SACtr` type, which is exported.  Methods
are written expecting and dispatched on this type.  Methods are also defined
for arrays of `SACtr`, `Array{SACtr}`, which allows for easy operations on
multiple traces.

The `SACtr` type has fields whose names and types correspond to SAC headers.
Do to access or set a header variable, just manipulate the object directly:

```julia
julia> t = SAC.sample();

julia> typeof(t)
SAC.SACtr

julia> t.delta
0.01f0

julia> t.delta = 0.02
0.02f0
```

(Note that SAC floating point headers are `Float32`s.)

The field `t` contains the trace as an `Array{Float32,1}`; manipulate this when
performing operations on the trace:

```julia
julia> t.t += 1 # Add 1 to all points, like SAC's "add" command
1000-element Array{Float32,1}:
 0.90272
 ⋮      
 0.9232 
```

You can also get or modify several header values at once by indexing via a
symbol; basically, using the name of the header variable preceded by a colon:

```julia
julia> typeof(T)
Array{SAC.SACtr,1}

julia> T[:t0]
5-element Array{Float32,1}:
 1.0
 2.0
 3.0
 4.0
 5.0

julia> T[:t0] += 2 # Move all time markers back by 2 s
5-element Array{Float32,1}:
 3.0
 4.0
 5.0
 6.0
 7.0

julia> T[:kstnm] = ["A1", "B2", "C3", "D4", "E5"]
5-element Array{ASCIIString,1}:
 "A1"
 "B2"
 "C3"
 "D4"
 "E5"
```

### Reading files
The `read` function is not exported to avoid name clashes, so one must call
`SAC.read()`.   For example, to load a single file, do

```julia
t = SAC.read("XM.A01E.HHZ.SAC")
```

This loads the file `XM.A01E.HHZ.SAC` into the `SACtr` object `t`.

### Reading with wildcards
As with SAC, one can use wildcards to read a set of files, e.g.:

```julia
T, filenames = read_wild("*Z.SAC")
```

An array of `SACtr`, `T`, is returned as well as a list of matching file names
in `filenames`.

### Writing files
Use the exported `write` method, passing a `SACtr` object and the file name, or
an array of `SACtr` and filenames

```julia
write(t, "file.SAC")
write(T, filenames)
```

### Processing
A number of common processing steps are already implemented as methods, such as
`lowpass!`, `taper!`, `envelope!` and so on.  In many cases, methods which have
similar names to SAC commands can also be used with the SAC short forms.  For
instance, `bandpass!` and `bp!` are the same.  Note that as is convention in
Julia, these commands end with an exclamation mark (!) and modify the trace.
Create a copy of the trace with `copy` if this is not desired.

### File endianness
[SAC/BRIS](http://www1.gly.bris.ac.uk/~george/sac-bugs.html) expects files to
always be in big-endian format;
[SAC/IRIS](http://ds.iris.edu/ds/nodes/dmc/software/downloads/sac/) expects them
in the same endianness as the machine.  SAC.jl is agnostic and will both read
and write in either endianness, but generally prefers to stick to big-endian,
for compatibilty with SAC/BRIS.

### Plotting
A companion repo, [`SACPlot`](https://github.com/anowacki/SACPlot.jl)
can be used to perform some of the plotting that SAC can do.


## Getting help
Functions are documented, so at the REPL type `?` to get a `help?>` prompt,
and type the name of the function:

```julia
help?> bandpass!
search: bandpass!

  bandpass!(::SACtr, c1, c2; ftype="butterworth", npoles=2, passes=1)

  Perform a bandpass filter on the SAC trace, between frequency corners c1 and c2.

  Select type of filter with ftype: current options are: butterworth. Set number of
  poles with npoles.

  passes may be 1 (forward) or 2 (forward and reverse).
```

### Documentation
Documentation is a work in progress, but all useful commands are documented.
To see the list of commands, check the code, or in the REPL type `SAC.` then
press tab a couple of times to see all the module methods and variables.
Calling up the interactive help will give a useful description of each.


## Dependencies
- [Glob.jl](https://github.com/vtjnash/Glob.jl)
- [DSP.jl](https://github.com/JuliaDSP/DSP.jl)
- [GreatCircle.jl](https://github.com/acrosby/GreatCircle.jl)

Install these using by doing
```julia
Pkg.add("Glob")
Pkg.add("DSP")
Pkg.clone("https://github.com/anowacki/GreatCircle.jl", "an/precompile")
```

Note that we need a version of GreatCircle which support precompilation, which
is why we need the `Pkg.clone` command.

## Other software

* If you use Fortran, then you should investigate the following modules:
  - [`sacio90`](https://github.com/jwookey/sacio90)
  - `f90sac` in the [seismo-fortran](https://github.com/anowacki/seismo-fortran)
    repo, inspired by the James Wookey original above.
* If you use Python, use [ObsPy](https://github.com/obspy/obspy/wiki).
* If you use MATLAB, use [`msac`](https://github.com/jwookey/msac).
