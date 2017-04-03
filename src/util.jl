# Utility routines

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
