# SAC.jl

[![Build Status](https://img.shields.io/travis/anowacki/SAC.jl.svg?style=flat-square&label=linux)](https://travis-ci.org/anowacki/SAC.jl)
[![Build status](https://img.shields.io/appveyor/ci/AndyNowacki/sac-jl.svg?style=flat-square&label=windows)](https://ci.appveyor.com/project/AndyNowacki/sac-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/anowacki/SAC.jl/badge.svg?branch=master)](https://coveralls.io/github/anowacki/SAC.jl?branch=master)

## Status

SAC.jl is now **deprecated** in favour of
[Seis.jl](https://github.com/anowacki/Seis.jl), and no new features will
be added to this package.

All SAC.jl functionality exists in Seis.jl, which also includes much
better [documentation](https://anowacki.github.io/Seis.jl/stable/)
and more IO options, and integrates with the wider Seis.jl
ecosystem (including [SeisSplit](https://github.com/anowacki/SeisSplit.jl),
[SeisTau](https://github.com/anowacki/SeisTau.jl) and
[Beamforming](https://github.com/anowacki/Beamforming.jl)).

## What is SAC.jl?
A [Julia](http://julialang.org) package for dealing with seismic data in the
[SAC](http://ds.iris.edu/files/sac-manual/manual/file_format.html) format, and
processing that data in a similar way to the SAC program:
either the [SAC/BRIS](http://www1.gly.bris.ac.uk/~george/sac-bugs.html) or
[IRIS](http://ds.iris.edu/ds/nodes/dmc/software/downloads/sac/) versions.


## How to install
Although not registered as an official package, SAC.jl can be added to your
Julia install like so:

On Julia v0.7 or v1.0 (press `]` to get to `pkg` mode first):

```julia
(v1.0) pkg> add https://github.com/anowacki/SAC.jl
```

(Or you can also do `import Pkg; Pkg.add("https://github.com/anowacki/SAC.jl")`.)

On Julia v0.6, versions before 0.3 can be installed like so:

```julia
Pkg.clone("https://github.com/anowacki/SAC.jl")
```

This will automatically install the depndencies you need.  You then need only do

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
These are accessed via `Symbol`s which are the name of header.  (To get a `Symbol`,
just write the name of the header with a colon in front.)

```julia
julia> t = SAC.sample()
SAC.SACtr:
    delta: 0.01
   depmin: -1.56928
   depmax: 1.52064
        b: 9.459999
		  ⋮
    kevnm: K8108838

julia> typeof(t)
SAC.SACtr

julia> t[:delta]
0.01f0

julia> t[:delta] = 0.02
0.02
```

(Note that SAC floating point headers are `Float32`s.)

The field `t` contains the trace as an `Array{Float32,1}`.  To change the trace,
just alter the `:t` index:

```julia
julia> t[:depmax]
1.52064f0

julia> t[:t] += 1;

julia> t.t
1000-element Array{Float32,1}:
 0.90272
 0.90272
 0.90144
 ⋮      
 0.92832
 0.9232 
 0.9232 

julia> t[:depmax]
2.52064f0
```

You can use the methods `+`, `-`, `*` and `/` to modify the traces without
needing to access `:t` directly, too:

```julia
julia> t == SAC.sample() + 1
true

julia> t == 1*t
true
```

You can also get or modify several header values at once:

```julia
julia> T = [SAC.sample() for _ in 1:5]
5-element Array{SAC.SACtr,1}:
 SAC.SACtr(delta=0.01, b=9.459999, npts=1000, kstnm=CDV, gcarc=3.357463, az=88.14708, baz=271.8529)
 SAC.SACtr(delta=0.01, b=9.459999, npts=1000, kstnm=CDV, gcarc=3.357463, az=88.14708, baz=271.8529)
 SAC.SACtr(delta=0.01, b=9.459999, npts=1000, kstnm=CDV, gcarc=3.357463, az=88.14708, baz=271.8529)
 SAC.SACtr(delta=0.01, b=9.459999, npts=1000, kstnm=CDV, gcarc=3.357463, az=88.14708, baz=271.8529)
 SAC.SACtr(delta=0.01, b=9.459999, npts=1000, kstnm=CDV, gcarc=3.357463, az=88.14708, baz=271.8529)

julia> typeof(T)
Array{SAC.SACtr,1}

julia> T[:t0] = 1:5 # Set time markers in t0
1:5

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

julia> T[:kstnm] = ["A1", "B2", "C3", "D4", "E5"] # Set station names
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
instance, `bandpass!` and `bp!` are the same.

Note that as is convention in Julia, these commands end with an exclamation
mark (!) and modify the trace in-place.  Copying versions of these commands are
available and do not have an exclamation mark (e.g., `lowpass`, `taper`, etc.).

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
search: bandpass! bandpass

  bandpass!(s::SACtr, c1, c2; ftype=:butterworth, npoles=2, passes=1) -> s

  Perform a bandpass filter on the SAC trace s, between frequency corners c1 and c2,
  returning the modified trace.

  Select type of filter with ftype: current options are: Symbol[:butterworth]. Set
  number of poles with npoles.

  passes may be 1 (forward) or 2 (forward and reverse).
```

### Documentation
Documentation is a work in progress, but all useful commands are documented.
To see the list of commands, check the code, or in the REPL type `SAC.` then
press tab a couple of times to see all the module methods and variables.
Calling up the interactive help will give a useful description of each.


## Other software

* If you use Fortran, then you should investigate the following modules:
  - [`sacio90`](https://github.com/jwookey/sacio90)
  - `f90sac` in the [seismo-fortran](https://github.com/anowacki/seismo-fortran)
    repo, inspired by the James Wookey original above.
* If you use Python, use [ObsPy](https://github.com/obspy/obspy/wiki).
* If you use MATLAB, use [`msac`](https://github.com/jwookey/msac).
