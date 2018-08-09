# Utility routines

"""
    angle_difference(α, β, degrees::Bool=true) -> Δ

Return the angular difference `Δ` between two angles `α` and `β`, in the direction
α → β (i.e., Δ = (β - α)).  This means that Δ is positive if the direction β is
clockwise of α, and negative otherwise.

Angles are assumed to be degreees unless `degrees` is `false`.
"""
function angle_difference(a::Real, b::Real, degrees::Bool=true)
    whole_circ = degrees ? 360. : 2π
    half_circ = whole_circ/2.
    mod(b - a + half_circ, whole_circ) - half_circ
end

"""
    flip_component!(s::SACtr) -> s

Flip a component so that it points in the opposite direction.
"""
function flip_component!(s::SACtr)
    s[:cmpaz] = mod(s[:cmpaz] + 180, 360)
    s[:cmpinc] = 180 - s[:cmpinc]
    s[:kcmpnm] = sacstring(Compat.round(s[:cmpaz], digits=2, base=10))
    multiply!(s, -1)
end

"""
    linear_regression(x, y)

Perform simple linear regression using Ordinary Least Squares. Returns `a` and `b` such
that `a + b*x` is the closest straight line to the given points `(x, y)`, i.e., such that
the squared error between `y` and `a + b*x` is minimized.
"""
function linear_regression(x::AbstractVector, y::AbstractVector)
    size(x) == size(y) || throw(DimensionMismatch("x and y must be the same size"))
    mx, my = mean(x), mean(y)
    b = covm(x, mx, y, my)/varm(x, mx)
    a = my - b*mx
    return a, b
end

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
        file_pattern = "*" * dir[end:end] * ".*"
        dir = dir[1:end-1]
    end
    A, filenames = read_wild(file_pattern, dir, echo=false)
    A[:kcmpnm] = [split(filename, ".")[end] for filename in filenames]
    A
end

"""
    swap_traces!(s1, s2) -> s2, s1

Swap the contents of two traces, `s1` and `s2`, by an in-place transposition of all
header and trace values.
"""
function swap_traces!(s1::SACtr, s2::SACtr)
    for f in fieldnames(s1)
        s1[f], s2[f] = s2[f], s1[f]
    end
    s1, s2
end

"""
    traces_are_orthogonal(s1::SACtr, s2::SACtr, tol=eps(SACFloat)) -> ::Bool

Return `true` if the two traces `s1` and `s2` have component azimuths
90° apart.
"""
traces_are_orthogonal(s1::SACtr, s2::SACtr, tol=eps(SACFloat)) =
    isapprox(SACFloat(abs(angle_difference(s1[:cmpaz], s2[:cmpaz]))), SACFloat(90), atol=tol)

"""
    time(::SACtr) -> t

Return a FloatRange `t` which contains the times for each sample of the SAC trace.
"""
time(s::SACtr) = s.b .+ (0:s.npts-1)*s.delta
