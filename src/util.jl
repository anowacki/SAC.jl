# Utility routines

"""
    isundefined(x) -> ::Bool

If the SAC value `x` is undefined, return `true`.

    isundefined(s::SACtr, x::Symbol) -> :: Bool

Test whether header `x` is undefined for trace `s`.
"""
isundefined(x::SACFloat) = x == sac_rnull
isundefined(x::SACInt) = x == sac_inull
isundefined(x::SACChar) = strip(x) == strip(sac_cnull)
isundefined(x::SACBool) = false
isundefined(s::SACtr, x::Symbol) = isundefined(getfield(s, x))

"""
    sample() -> ::SACtr

Return some sample SAC data, a regional earthquake arrival.
"""
function sample()
    file = joinpath(dirname(@__FILE__()), "../data/seis.sac")
    return SAC.read(file)
end

"""
    time(::SACtr) -> t

Return a FloatRange `t` which contains the times for each sample of the SAC trace.
"""
time(s::SACtr) = s.b + (0:s.npts-1)*s.delta
