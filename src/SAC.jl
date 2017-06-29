"""
SAC.jl provides routines for dealing with SAC-formatted time series files,
including reading, writing, filtering, mean removal, rotating, and so on.
Sister library SACPlot.jl can be used for plotting.
"""
module SAC

__precompile__()

import DSP
import Glob
import Base: ==, copy, diff, getindex, fft, setindex!, time, write

export
    SACtr,
    bandpass!,
    bp!,
    copy,
    cut!,
    diff!,
    differentiate!,
    envelope!,
    fft,
    highpass!,
    hp!,
    integrate!,
    interpolate!,
    lp!,
    lowpass!,
    read_wild,
    rmean!,
    rotate_through!,
    rotate_through,
    rotate_to_gcp!,
    rotate_to_gcp,
    rtrend!,
    taper!,
    time,
    tshift!,
    write

include("constants.jl")
include("types.jl")
include("great_circle.jl")
include("io.jl")
include("operations.jl")
include("filtering.jl")
include("util.jl")

# Submodules
include("Stack.jl")
using .Stack
export stack

# Build all copying routines
#FIXME: Our definition of SAC.diff(::Array{SACtr}) is ambiguous with
#       Base.diff(::AbstractArray), meaning we get a
#           WARNING: imported binding for diff overwritten in module SAC
#       every time we `import SAC`.
"""Dict with keys given by name of each function to have a copying version.
   Where an abbreviated version exists, that is given as the value; otherwise
   the value is `nothing`"""
const copying_funcs = Dict(
    :add! => nothing,
    :bandpass! => :bp!,
    :cut! => nothing,
    :differentiate! => :diff!,
    :divide! => :div!,
    :envelope! => nothing,
    :flip_component! => :flip!,
    :highpass! => :hp!,
    :integrate! => :int!,
    :interpolate! => nothing,
    :lowpass! => :lp!,
    :multiply! => :mul!,
    :rmean! => nothing,
    :rtrend! => nothing,
    :taper! => nothing,
    :tshift! => nothing,
    )
for (name, abbrev) in copying_funcs
    new_name = Symbol(string(name)[1:end-1])
    new_abbrev = Symbol(string(abbrev)[1:end-1])
    @eval begin
        function ($new_name)(s::Union{SACtr,Array{SACtr}}, args...; kwargs...)
            s_new = deepcopy(s)
            $(name)(s_new, args...; kwargs...)
            s_new
        end
        @doc """
            $($new_name)(s::Union{SACtr,Array{SACtr}}, args...; kwargs...) -> s_new

        Copying version of `$($name)` which returns modified version(s) of the trace(s)
        in `s`, leaving the originals unaltered.  See docs of `$($name)` for details.
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

include("precompile.jl")

end # module SAC
