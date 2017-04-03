"""
# SAC.Stack

SAC.Stack is a submodule of SAC which implements simple stacking routines
to create stacks of SAC traces.

### Using

Pass an array of `SACtr` objects to the `stack` routine, with a range of times
over which to sum, and optionally times on which to align the traces, weights
to apply and the method to use.

### Available stacking methods

- `:linear` [default].  A linear (optionally weighted) stack.
"""
module Stack

import SAC

export
    stack

"List of available stacking routines"
const available_stack_methods = [:linear]

"""
    stack(S::Array{SACtr}, time_range, align=zeros(length(S)), weight; method=:linear) -> s::SACtr

Return the linear stack of all traces in `S`, aligned in time on the
value in `align`, between `times[1]` and `times[end]`.

`times` can be a Tuple, Range or Array; in all cases, only the first and last values
are used.  **N.B. If using a non-integer Range, this may not include the end value
you expect.  E.g., `0.1:0.9` by default contains only the value `0.1` because the
default step is `1`.**

`align` can be a header symbol, or an array of equal length to `S`.  If it is a Symbol,
then the values in the headers are used.

All traces must have the same sampling frequency, and an error is thrown
if not all traces can be included in the stack.

If provided, traces are weighted by the values in `weight`, which are normalised
to sum to unity.
"""
function stack(S::Array{SAC.SACtr}, time_range,
               align::Union{Symbol,AbstractArray}=zeros(SAC.SACFloat, length(S)),
               weight=ones(SAC.SACFloat, length(S));
               method=:linear)
    method in available_stack_methods ||
        throw(ArgumentError("Stacking method '$(method)' is not valid.  " *
                            "Choose from: $(available_stack_methods)"))
    # Convert to a vector of times if given as a symbol
    if typeof(align) == Symbol
        align in SAC.sac_float_hdr || error("Header '$(align)' is not a float header")
        any(S[align] .== SAC.sac_rnull) && error("Not all traces have header $align")
        align = S[align]
    end
    all(S[1].delta .== S[:delta]) || error("All traces must have same delta")
    delay = -align
    t1, t2 = time_range[1], time_range[end]
    weight /= sum(weight) # Normalise weights
    if lowercase(string(method)) == "linear"
        stack_trace = _stack_linear([w.*@view(s.t[:]) for (s, w) in zip(S, weight)],
            S[1].delta, S[:b], delay, t1, t2)
    end
    out = SAC.SACtr(S[1].delta, length(stack_trace), t1)
    out.t[:] .= stack_trace[:]
    SAC.update_headers!(out)
    out
end

"""
    _stack_linear(a, delta, b, delay, t1, t2) -> stacked_trace::Vector{eltype(a)}

Return an array containing the linear stack of the vector of vectors `a`, whose
start times are given in `b`, and where the stack is returned between times `t1` and `t2`
(which can be negative if delay is also negative).  `delta` is the sampling interval
for all the traces.

The trace starts at time `t1` and has sample spacing `delta`.
"""
function _stack_linear{A<:AbstractArray}(a::Vector{A}, delta::Real, b::Vector,
        delay::Vector, t1::Real, t2::Real)
    N, npts, b_shift, stack_trace = _stack_array(a, delta, b, delay, t1, t2)
    # stack_trace checks that we are always in array bounds
    @inbounds for j in 1:N
        ip1 = round(Int, (t1 - b_shift[j])/delta) + 1
        @simd for i in ip1:(ip1+npts-1)
            stack_trace[i-ip1+1] += a[j][i]
        end
    end
    stack_trace # Do not normalise as weightings are applied in `stack`
end

"""
    _stack_array(a::Vector{Vector{<:Real}}, delta, b::Vector, delay::Vector, t1, t2) -> N, npts, b_shifted, stack_trace

Return the number of traces `N`, the number of points in the output stack trace,
`npts`, an array of shifted trace start times `b_shifted` and a zeroed array
for the stack itself, `stack_trace`.  This routine also checks that the input
arrays are the required length and no points are missing from the desired stack.
"""
function _stack_array{A<:AbstractArray}(a::Vector{A}, delta, b, delay, t1, t2)
    N = length(a)
    length(b) == length(delay) == N ||
        throw(ArgumentError("Lengths of `a`, `b` and `delay` must all be the same"))
    t1 < t2 || throw(ArgumentError("`t2` is before or the same as `t1`"))
    b_shift = b .+ delay
    e_shift = b_shift .+ (length.(a) .- 1).*delta
    all(b_shift .<= t1) && all(t2 .<= e_shift) ||
        error("Not all traces have data within time range $(t1):$(t2) s
        range in b: $(extrema(b_shift))
        range in e: $(extrema(e_shift))")
    npts = round(Int, (t2 - t1)/delta) + 1
    stack_trace = zeros(eltype(eltype(a)), npts)
    N, npts, b_shift, stack_trace
end

end # module
