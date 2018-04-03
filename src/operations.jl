# Trace operations

# Linear operations on the trace
"""
    add!(s::SACtr, value) -> s

Add a constant `value` to a SAC trace `s` in place.
"""
function add!(s::SACtr, val)
    s.t .= s.t .+ val
    update_headers!(s)
end

"""
    add(s::SACtr, value) -> s_new
    s + value
    value + s

Add a constant `value` to a SAC trace `s` and return a new trace `s_new`.
"""
add(s::SACtr, val) = add!(deepcopy(s), val)

Base.:+(s::SACtr, val) = add(s, val)
Base.:+(val, s::SACtr) = add(s, val)

"""
    subtract!(s::SACtr, value) -> s

Substract a constant `value` to a SAC trace `s` in place.
"""
subtract!(s::SACtr, val) = add!(s, -val)

"""
    subtract(s::SACtr, value) -> s_new
    s - value
    value - s

Subtract a constant `value` from a SAC trace `s` (or vice versa) and
return a new trace `s_new`.
"""
subtract(s::SACtr, val) = subtract!(deepcopy(s), val)

Base.:-(s::SACtr, val) = add(s, -val)
Base.:-(val, s::SACtr) = add(-s, val)
Base.:-(s::SACtr) = multiply(s, -1)

"""
    multiply!(s::SACtr, value) -> s
    mul!(s::SACtr, value) -> s

Multiply the values in a SAC trace `s` by `value`.
"""
function multiply!(s::SACtr, val)
    s.t .= s.t .* val
    update_headers!(s)
end
const mul! = multiply!

"""
    multiply(s::SACtr, value) -> s_new
    s * value
    value * s

Multiply the trace `s` by `value` and return a new trace `s_new`.
"""
multiply(s::SACtr, val) = multiply!(deepcopy(s), val)

Base.:*(s::SACtr, val) = multiply(s, val)
Base.:*(val, s::SACtr) = multiply(s, val)

"""
    divide!(s::SACtr, value) -> s

Divide the values in a SAC trace `s` by `value` in place.
"""
function divide!(s::SACtr, value)
    value != 0 || error("SAC.divide!: Cannot divide by 0")
    multiply!(s, 1.0/value)
end

"""
    divide(s::SACtr, value) -> s_new
    s / value
    value + s

Divide a SAC trace `s` by `value` and a return a new trace `s_new`.
"""
function divide!(value, s::SACtr)
    s.t .= value ./ s.t
    update_headers!(s)
end

divide(s::SACtr, value) = divide!(deepcopy(s), value)
divide(value, s::SACtr) = divide!(value, deepcopy(s))

Base.:/(s::SACtr, val) = divide(s, val)
Base.:/(val, s::SACtr) = divide(val, s)

"""
    modify!(s::SACtr, f) -> s

Modify a SACtr trace `s` in place with the function `f`, ensuring
that the headers are updated in the process.  `f` is applied to each
element of `s`.

For example, to square each value of a trace:

```
julia> s = SACtr(1:10, 1);

julia> modify!(s, x->x^2);

julia> s[:depmax]
100.0f0
```
"""
function modify!(s::SACtr, f)
    s.t .= f.(s.t)
    update_headers!(s)
end

"""
    modify(s::SACtr, f) -> s_new

Modify a copy of the SACtr trace `s` by the function `f`.

For example, to square each value of a trace:

```julia
julia> s = SACtr(1:10, 1);

julia> modify(s, x->x^2)[:depmax]
100.0f0
```
"""
modify(s::SACtr, f) = modify!(deepcopy(s), f)

# Linear operations on multiple traces
check_traces_combine(s1, s2) = begin s1.b ≈ s2.b && s1.npts == s2.npts ||
    throw(ArgumentError("Traces must have same strart time and length")) end

"""
    add(s1::SACtr, s2::SACtr) -> total::SACtr
    s1 + s2

Return a new SACtr `total` which is the sum of the two traces `s1` and `s2`.
"""
function add(s1::SACtr, s2::SACtr)
    check_traces_combine(s1, s2)
    out = deepcopy(s1)
    out.t .= s1.t .+ s2.t
    update_headers!(out)
end
Base.:+(s1::SACtr, s2::SACtr) = add(s1, s2)

"""
    subtract(s1::SACtr, s2::SACtr) -> difference::SACtr
    s1 - s2

Return a new SACtr `difference` which is the second trace `s2` subtraced
from the first `s1`.
"""
function subtract(s1::SACtr, s2::SACtr)
    check_traces_combine(s1, s2)
    out = deepcopy(s1)
    out.t .= s1.t .- s2.t
    update_headers!(out)
end
Base.:-(s1::SACtr, s2::SACtr) = subtract(s1, s2)

"""
    multiply(s1::SACtr, s2::SACtr) -> product::SACtr
    s1 * s2

Return a new SACtr `product` which is the result of multiplying two traces `s1`
and `s2` together.
"""
function multiply(s1::SACtr, s2::SACtr)
    check_traces_combine(s1, s2)
    out = deepcopy(s2)
    out.t .= s1.t .* s2.t
    update_headers!(out)
end
Base.:*(s1::SACtr, s2::SACtr) = multiply(s1, s2)

"""
    divide(s1::SACtr, s2::SACtr) -> quotient::SACtr
    s1 / s2

Return a new SACtr `quotient` which is the result of dividing the first trace `s1`
by the second `s2`.
"""
function divide(s1::SACtr, s2::SACtr)
    check_traces_combine(s1, s2)
    out = deepcopy(s1)
    out.t .= s1.t ./ s2.t
    update_headers!(out)
end
Base.:/(s1::SACtr, s2::SACtr) = divide(s1, s2)


# Temporal operations
"""
    cut!(s::Union{SACtr,AbstractArray{SACtr}}, b, e)
    cut!(s::Union{SACtr,AbstractArray{SACtr}}, beg_head::Symbol, b, end_head::Symbol, e)

In the first form, cut a trace or array of traces `s` in memory between times
`b` and `e`, relative to the O marker.

In the second form, cut based on headers `beg_head` and `end_head`, from
`beg_head+b` to `end_head+e`.

    cut!(s::AbstractArray{SACtr}, a::Array, b::Array)

Cut the array of traces `s` between the times in arrays `b` and `e`, which must be
the same length as `s`.
"""
function cut!(s::SACtr, b::Real, e::Real)
    if b < s.b
        info("SAC.cut!: beginning cut is before start of trace.  Setting to $(s.b).")
        b = s.b
    end
    b > s.e && error("SAC.cut!: end cut time is later than end of trace.")
    if e > s.e
        info("SAC.cut!: end cut is after end of trace.  Setting to $(s.e).")
        e = s.e
    end
    e < s.b && error("SAC.cut!: end time is earlier than start of trace.")
    ib = round(Int, (b - s.b)/s.delta) + 1
    ie = s.npts - round(Int, (s.e - e)/s.delta)
    s.t = s.t[ib:ie]
    s.b, s.e = s.b + (ib - 1)*s.delta, s.b + (ie - 1)*s.delta
    s.npts = ie - ib + 1
    update_headers!(s)
end
function cut!(s::SACtr, bh::Symbol, b::Real, eh::Symbol, e::Real)
    isundefined(s[bh]) && error("Header $bh is not defined")
    isundefined(s[eh]) && error("Header $eh is not defined")
    b, e = s[bh] + b, s[eh] + e
    cut!(s, b, e)
end

# Array version of cut!
function cut!(a::AbstractArray{SACtr}, b::Number, e::Number)
    for s in a
        SAC.cut!(s, b, e)
    end
    a
end
function cut!(a::AbstractArray{SACtr}, bh::Symbol, b::Real, eh::Symbol, e::Real)
    for s in a cut!(s, bh, b, eh, e) end
    a
end

function cut!(a::AbstractArray{SACtr}, b::Array{B}, e::Array{E}) where {B<:Real,E<:Real}
    @assert length(a) == length(b) == length(e) "Arrays `a`, `b` and `e` must be the same length"
    for (s, beg, en) in zip(a, b, e)
        SAC.cut!(s, beg, en)
    end
    a
end

"""
    differentiate!(s::SACtr, npoints::Integer=2)

Differentiate the SAC trace `s`, replacing it with its time derivative `dsdt`.
Select the mode of numerical differentiation with `npoints`.

### Available algorithms

- `npoints == 2`: Two-point.  `dsdt.t[i] = (s.t[i+1] - s.t[i])/s.delta`.
  Non-central difference, so `s.b` is increased by half `s.delta`.  `npts` is
  reduced by 1.
- `npoints == 3`: Three-point. `dsdt.t[i] = (s.t[i+1] - s.t[i-1])/(2 * s.delta)`.
  Central difference.  `s.b` is increased by `s.delta`; `npts` reduced by 2.
- `npoints == 5`: Five-point. `dsdt.t[i] =
  (2/3)*(s.t[i+1] - s.t[i-1])/s.delta - (1/12)*(s.t[i+2] - s.t[i-2])/s.delta`.
  Central difference.  `s.b` is increased by `2s.delta`; `npts` reduced by 4.
"""
function differentiate!(s::SACtr, npoints::Integer=2)
    npoints in (2, 3, 5) ||
        throw(ArgumentError("`npoints` cannot be $(npoints); must be one of (2, 3, 5)"))
    if npoints == 2
        t = Vector{SACFloat}(undef, s.npts - 1)
        @inbounds for i in 1:(s.npts-1)
            s.t[i] = (s.t[i+1] - s.t[i])/s.delta
        end
        pop!(s.t)
        s.npts -= 1
        s.b += s.delta/2
    elseif npoints == 3
        @inbounds for i in 2:(s.npts-1)
            s.t[i-1] = (s.t[i+1] - s.t[i-1])/(2*s.delta)
        end
        pop!(s.t); pop!(s.t)
        s.npts -= 2
        s.b += s.delta
    elseif npoints == 5
        t1 = (s.t[3] - s.t[1])/(2*s.delta)
        t2 = (s.t[end] - s.t[end-2])/(2*s.delta)
        d1 = 2/(3*s.delta)
        d2 = 1/(12*s.delta)
        t_minus_2 = s.t[1]
        t_minus_1 = s.t[2]
        t = s.t[3]
        t_plus_1 = s.t[4]
        @inbounds for i in 2:(s.npts-3)
            t_plus_2 = s.t[i+3]
            s.t[i] = d1*(t_plus_1 - t_minus_1) - d2*(t_plus_2 - t_minus_2)
            t_minus_2 = t_minus_1
            t_minus_1 = t
            t = t_plus_1
            t_plus_1 = t_plus_2
        end
        s.t[1] = t1
        s.t[end-2] = t2
        pop!(s.t); pop!(s.t)
        s.npts -= 2
        s.b += s.delta
    end
    update_headers!(s)
end

"""
    envelope!(::SACtr)

Find the envelope of a SAC trace
"""
function envelope!(a::AbstractArray{SACtr})
    for s in a
        s.t = abs(DSP.hilbert(s.t))
    end
    update_headers!(a)
end
envelope!(s::SACtr) = envelope!([s])

"""
    fft(s::SACtr) -> f, S

Return the Fourier-transformed trace from `s` as `S`, with the frequencies
which correspond to each point in `f`.
"""
function fft(s::SACtr)
    N = round(Int, s.npts/2) + 1
    fmax = 1.0/(s.npts*s.delta)
    f = collect(1:N)*fmax
    S = Base.fft(s.t)[1:N]
    return f, S
end

function fft(a::AbstractArray{SACtr})
    n = length(a)
    f, S = Array{Array}(n), Array{Array}(n)
    for i = 1:n
        f[i], S[i] = fft(a[i])
    end
    return f, S
end

"""
    integrate!(s::SACtr, method=:trapezium)

Replace `s` with its time-integral.  This is done by default using the trapezium rule.
Use `method=:rectangle` to use the rectangle rule.

If `method==:trapezium` (the default), then `s.npts` is reduced by one and `s.b` is
increased by `s.delta/2`.
"""
function integrate!(s::SACtr, method::Symbol=:trapezium)
    method in (:trapezium, :rectangle) ||
        throw(ArgumentError("`methodod` must by one of `:trapezium` or `:rectangle` " *
                            "(got '$method')"))
    if method == :trapezium
        total = zero(s.t[1])
        h = s.delta/2
        @inbounds for i in 1:(s.npts-1)
            total += h*(s.t[i] + s.t[i+1])
            s.t[i] = total
        end
        s.npts -= 1
        pop!(s.t)
        s.b += s.delta/2
    elseif method == :rectangle
        h = s.delta
        @inbounds for i in 2:s.npts
            s.t[i] = h*s.t[i] + s.t[i-1]
        end
    end
    update_headers!(s)
end
const int! = integrate!

"""
    interpolate!(::SACtr, npts=npts)
    interpolate!(::SACtr, delta=delta)
    interpolate!(::SACtr, n=n)

Resample a SAC trace by supplying one of three things:

* A new number of samples (`npts`)
* A new sampling interval (`delta` in seconds)
* A multiple by which to increase the sampling (`n`)

Interpolation is performed using quadratic splines using the `Dierckx` package.
"""
function interpolate!(s::SACtr; npts::Integer=0, delta::Real=0.0, n::Integer=0)
    # Calculate new points at which to evaluate time series
    interp_t = if npts != 0
        npts >= 0 || error("`npts` cannot be negative")
        delta = (s.e - s.b)/(npts - 1)
        s.b + (0:(npts-1))*delta
    elseif delta != 0.0
        delta >= 0.0 || error("`delta` cannot be negative")
        delta < (s.e - s.b) || error("`delta`")
        times = s.b:delta:s.e
        npts = length(times)
        times
    elseif n != 0
        n > 0 || error("`n` cannot be negative")
        npts = (s.npts - 1)*n + 1
        delta = (s.e - s.b)/(npts - 1)
        s.b + (0:(npts-1))*delta
    else
        error("Must supply one keyword argument of `npts`, `n` or `delta`")
    end
    @assert npts == length(interp_t)
    # Create fit using degree-2 Bsplines
    spl = Dierckx.Spline1D(collect(SAC.time(s)), s.t, k=2)
    s.t = Dierckx.evaluate(spl, interp_t)
    s.npts = npts
    s.delta = delta
    update_headers!(s)
end

"""
    rmean!(::SACtr)

Remove the mean in-place for a SAC trace.
"""
function rmean!(s::SACtr)
    s.t .= s.t .- mean(s.t)
    update_headers!(s)
end
function rmean!(a::AbstractArray{SACtr})
    for s in a
        rmean!(s)
    end
    a
end

"""
    rotate_through!(::SACtr, ::SACtr, phi)
    rotate_through!(::AbstractArray{SACtr}, phi)

In the first form, with two SAC traces which are horizontal and orthgonal, rotate
them *clockwise* by `phi`° about the vertical axis.

In the second form, rotate each sequential pair of traces (i.e., indices 1 and 2,
3 and 4, ..., end-1 and end).

This is a reference frame transformation (passive rotation) and hence particle motion
will appear to rotate anti-clockwise.
"""
function rotate_through!(s1::SACtr, s2::SACtr, phi)
    if !traces_are_orthogonal(s1, s2)
        error("SAC.rotate_through!: traces must be orthogonal")
    elseif s1.npts != s2.npts
        error("SAC.rotate_through!: traces must be same length")
    elseif s1.delta != s2.delta
        error("SAC.rotate_through!: traces must have same delta")
    end
    phir = deg2rad(phi)
    cosp, sinp = cos(phir), sin(phir)
    @inbounds for i = 1:s1.npts
        s1.t[i], s2.t[i] = cosp*s1.t[i] - sinp*s2.t[i], sinp*s1.t[i] + cosp*s2.t[i]
    end
    for t in (s1, s2)
        setfield!(t, :cmpaz, SAC.SACFloat(mod(getfield(t, :cmpaz) + phi, 360.)))
        setfield!(t, :kcmpnm, SAC.sacstring(round(getfield(t, :cmpaz), 2)))
        update_headers!(t)
    end
    s1, s2
end
function rotate_through!(a::AbstractArray{SACtr}, phi::Real)
    length(a)%2 != 0 && error("SAC.rotate_through!: Array of traces must be a multiple of two long")
    for i = 1:length(a)÷2
        rotate_through!(a[2*i - 1], a[2*i], phi)
    end
    a
end
function rotate_through!(a::AbstractArray{SACtr}, phi::AbstractArray{<:Real})
    length(a)%2 != 0 &&
        throw(ArgumentError("SAC.rotate_through!: Array of traces must be a multiple of two long"))
    if length(phi) == length(a)
        for i in 1:length(a)÷2
            rotate_through!(a[2i-1], a[2i], phi[2i])
        end
    elseif length(phi) != length(a)÷2 ||
        for i in 1:length(a)÷2
            rotate_through!(a[2i-1], a[2i], phi[i])
        end
    else
        throw(ArgumentError("Length of `phi` must be the same as the number of traces " *
            "or half the number of traces"))
    end
    a
end

"""
    rotate_through(s1::SACtr, s2::SACtr, phi) -> new_s1, new_s2

Copying version of `rotate_through` which returns modified versions of the traces
in `s1` and `s2`, leaving the originals unaltered.  See docs of `rotate_through!` for details.
"""
function rotate_through(s1::T, s2::T, phi) where {T<:Union{SACtr,AbstractArray{SACtr}}}
    s1_new, s2_new = deepcopy(s1), deepcopy(s2)
    rotate_through!.(s1_new, s2_new, phi)
    s1_new, s2_new
end
rotate_through(a::AbstractArray{SACtr}, phi) = rotate_through!(deepcopy(a), phi)

"""
    rotate_to_gcp!(s1::SACtr, s2::SACtr, reverse=false) -> s1, s2

Rotate a pair of SAC traces in place, and return them, so that s1 points along
the radial direction (the backazimuth plus 180°), and s2 is 90° clockwise from
that.

If `reverse` is true, then s2 is rotated to be 90° anticlockwise, so that the
polarity is reversed.

The component names of the radial and transverse traces are updated to be
'R', and either 'T' or '-T' respectively for normal and reverse polarity.

    rotate_to_gcp!(a::AbstractArray{SACtr}, reverse=false) -> a

Rotate traces where each subsuequent pair of traces in the array are considered
as the two horizontal components
"""
function rotate_to_gcp!(s1::SACtr, s2::SACtr, reverse::Bool=false)
    s2_is_clockwise_of_s1 = angle_difference(s1[:cmpaz], s2[:cmpaz]) > 0
    s2_is_clockwise_of_s1 || ((s1, s2) = (s2, s1))
    !isundefined(s1[:baz]) && !isundefined(s2[:baz]) ||
        throw(ArgumentError("Backazimuth is not defined for both traces"))
    s1[:baz] ≈ s2[:baz] || throw(ArgumentError("Backazimuth not the same for both traces"))
    phi = mod(s1[:baz] + 180 - s1[:cmpaz], 360)
    rotate_through!(s2, s1, phi) # Checks for orthogonality
    reverse && flip_component!(s2)
    s1[:kcmpnm] = sacstring("Radial")
    s2[:kcmpnm] = sacstring((reverse ? "-" : "")*"Trans")
    s1, s2
end
function rotate_to_gcp!(a::AbstractArray{SACtr}, reverse::Bool=false)
    length(a)%2 == 0 || throw(ArgumentError("Array of traces must be a multiple of two long"))
    for i in 1:length(a)÷2
        rotate_to_gcp!(a[2i-1], a[2i], reverse)
        # FIXME: This shouldn't be necessary, but is probably due to swapping
        # of the components in rotate_to_gcp!(s1, s2) above.
        if strip(a[2i-1][:kcmpnm]) == "Trans" a[2i-1], a[2i] = a[2i], a[2i-1] end
    end
    a
end

rotate_to_gcp(s1::SACtr, s2::SACtr, args...) = rotate_to_gcp!(deepcopy(s1), deepcopy(s2), args...)
rotate_to_gcp(a::AbstractArray{SACtr}, args...) = rotate_to_gcp!(deepcopy(a), args...)

"""
    rtrend!(::SACtr)

Remove the trend from a SAC trace in place.
"""
function rtrend!(s::SACtr)
    t = time(s)
    x0, x1 = linreg(t, s.t)
    s.t .= s.t .- (x0 .+ x1.*t)
    update_headers!(s)
end
function rtrend!(a::AbstractArray{SACtr})
    for s in a
        rtrend!(s)
    end
    a
end

"""
    taper!(s::SACtr, width=0.05, form=:hanning)
    taper!(S::AbstractArray{SACtr}, width=0.05, form=:hanning)

Apply a symmetric taper to each end of the data in SAC trace `s` or traces `S`.

`form` may be one of `:hanning`, `:hamming` or `:cosine`.

`width` represents the fraction (at both ends) of the trace tapered, up to 0.5.
"""
function taper!(s::SACtr, width=0.05, form=:hanning::Symbol)
    form in [:hamming, :hanning, :cosine] ||
        error("SAC.taper!: `form` must be one of `:hamming`, `:hanning` or `:cosine`")
    0 < width <= 0.5 || error("SAC.taper!: width must be between 0 and 0.5")
    n = max(2, floor(Int, (s.npts + 1)*width))

    if form in [:hamming, :hanning]
        omega = SAC.SACFloat(pi/n)
        if form == :hanning
            f0 = f1 = SAC.SACFloat(0.50)
        elseif form == :hamming
            f0 = SAC.SACFloat(0.54)
            f1 = SAC.SACFloat(0.46)
        end

        @inbounds for i in 0:n-1
            amp = f0 - f1*cos(omega*SAC.SACFloat(i))
            j = s.npts - i
            s.t[i+1] *= amp
            s.t[j] *= amp
        end
    end

    if form == :cosine
        omega = SAC.SACFloat(pi/(2*n))
        @inbounds for i in 0:n-1
            amp = sin(omega*i)
            j = s.npts - i
            s.t[i+1] *= amp
            s.t[j] *= amp
        end
    end

    SAC.update_headers!(s)
end
taper!(S::AbstractArray{SACtr}, width=0.05, form::Symbol=:hamming) =
    (for s in S taper!(s, width, form) end; S)

"""
    tshift!(::SACtr, tshift, wrap=true)

Shift a SAC trace backward in time by `t` seconds.

If `wrap` true (default), then points which move out the back of the trace
are added to the front (and vice versa).  Setting it to false instead pads the
trace with zeroes.
"""
function tshift!(s::SACtr, tshift::Number, wrap=true)
    n = round(Int, tshift/s.delta)
    if n == 0
        sac_verbose && info("SAC.tshift!: t ($tshift) is less than delta ($(s.delta)) so no shift applied")
        return
    end
    s.t = circshift(s.t, n)
    if !wrap
        n > 0 ? s.t[1:n] = 0. : s.t[end+n+1:end] = 0.
    end
    update_headers!(s)
end

"""
    update_headers!(s::SACtr)

Ensure that header values which are based on the trace or other header values
are consistent, such as `depmax`.  Should be called after any operation on the trace
`s.t`.
"""
function update_headers!(s::SACtr)
    s.depmax = maximum(s.t)
    s.depmin = minimum(s.t)
    s.depmen = mean(s.t)
    s.e = s.b + s.delta*(s.npts - 1)
    s
end

function update_headers!(a::AbstractArray{SACtr})
    for s in a
        update_headers!(s)
    end
    a
end
