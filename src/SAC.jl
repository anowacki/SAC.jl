__precompile__()

"""
SAC.jl provides routines for dealing with SAC-formatted time series files,
including reading, writing, filtering, mean removal, rotating, and so on.
Sister library SACPlot.jl can be used for plotting.
"""
module SAC

using Compat
using Compat.LinearAlgebra

@static if VERSION >= v"0.7-"
    import StatsBase: linreg
end

import DSP, Dierckx
import Glob
import Base: ==, copy, getindex, fft, setindex!, time, write

export
    SACtr,
    add!,
    add,
    bandpass!,
    bp!,
    copy,
    cut!,
    differentiate!,
    divide!,
    divide,
    envelope!,
    fft,
    highpass!,
    hp!,
    integrate!,
    interpolate!,
    lp!,
    lowpass!,
    modify!,
    modify,
    multiply!,
    multiply,
    normalise!,
    normalise,
    read_wild,
    rmean!,
    rotate_through!,
    rotate_through,
    rotate_to_gcp!,
    rotate_to_gcp,
    rtrend!,
    sample,
    stalta,
    subtract!,
    subtract,
    taper!,
    time,
    tshift!,
    window_average,
    write

include("constants.jl")
include("types.jl")
include("great_circle.jl")
include("io.jl")
include("operations.jl")
include("filtering.jl")
include("util.jl")
include("statistics.jl")

# Submodules
include("Stack.jl")
using .Stack
export stack

# Build all copying routines
"""Dict with keys given by name of each function to have a copying version.
   Where an abbreviated version exists, that is given as the value; otherwise
   the value is `nothing`"""
const copying_funcs = Dict(
    :bandpass! => :bp!,
    :cut! => nothing,
    :differentiate! => nothing,
    :envelope! => nothing,
    :flip_component! => :flip!,
    :highpass! => :hp!,
    :integrate! => :int!,
    :interpolate! => nothing,
    :lowpass! => :lp!,
    :rmean! => nothing,
    :rtrend! => nothing,
    :taper! => nothing,
    :tshift! => nothing,
    )
for (name, abbrev) in copying_funcs
    new_name = Symbol(string(name)[1:end-1])
    new_abbrev = Symbol(string(abbrev)[1:end-1])
    @eval begin
        function ($new_name)(s::Union{SACtr,AbstractArray{SACtr}}, args...; kwargs...)
            s_new = deepcopy(s)
            $(name)(s_new, args...; kwargs...)
            s_new
        end
        @doc """
        $(@doc $name)
            $($new_name)(s::Union{SACtr,AbstractArray{SACtr}}, args...; kwargs...) -> s_new

        Copying version of `$($name)` which returns modified version(s) of the trace(s)
        in `s`, leaving the originals unaltered.
        """ $new_name
        export $new_name
    end
    if abbrev != nothing
        @eval begin
            const $new_abbrev = $new_name
            export $new_abbrev
        end
    end
end

# Build methods which can be used with chaining (`|>`)
"""Dict with keys given by name of each function to have a chaining version.
   Where an abbreviated version exists, that is given as the value; otherwise
   the value is nothing."""
const chaining_funcs = Dict(
    copying_funcs...,
    :normalise! => nothing
    )

for (name, abbrev) in copying_funcs
    new_name = Symbol(string(name)[1:end-1])
    @eval begin
        ($name)(args...; kwargs...) = s -> ($name)(s, args...; kwargs...)
        ($new_name)(args...; kwargs...) = s -> ($new_name)(s, args...; kwargs...)
    end
end

include("precompile.jl")

end # module SAC
