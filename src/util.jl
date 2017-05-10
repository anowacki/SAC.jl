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

"Directory containing sample data"
const SAMPLE_DATA_DIR = joinpath(dirname(@__FILE__()), "..", "data")

"""
    sample() -> ::SACtr
    sample(kind::Symbol) -> ::Array{SACtr}

Return some sample SAC data.

With no arguments, `sample` gives one trace from a local earthquake recorded in
California.

In the second form, a set of traces is returned according to the table below:

|`kind`|Description|
|:-----|:----------|
|`:local`|Livermore Valley, CA.  9 3-component stations|
|`:regional`|Nevada.  4 3-component stations|
|`:teleseism`|**Mid-period** recording of Eureka, CA event.  4 3-c stations|
|`:teleseisl`|**Long-period** recording of Eureka, CA event.  4 3-c stations|
|`:array`|Deep Fiji event.  60 vertical stations in the UK|
"""
function sample()
    file = joinpath(SAMPLE_DATA_DIR, "seis.sac")
    SAC.read(file)
end
const SAMPLE_DATA_KINDS = [:local, :regional, :teleseism, :teleseisl, :array]
function sample(kind::Symbol)
    file_pattern = "*"
    kind in SAMPLE_DATA_KINDS ||
        throw(ArgumentError("Sample data kind must be one of: $SAMPLE_DATA_KINDS"))
    dir = joinpath(SAMPLE_DATA_DIR, string(kind))
    if kind in [:teleseism, :teleseisl]
        file_pattern = "*" * dir[end] * ".*"
        dir = dir[1:end-1]
    end
    A, f = read_wild(file_pattern, dir, echo=false)
    A
end


"""
    time(::SACtr) -> t

Return a FloatRange `t` which contains the times for each sample of the SAC trace.
"""
time(s::SACtr) = s.b + (0:s.npts-1)*s.delta
