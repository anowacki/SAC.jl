# Definition, constructors and other operations only to do with the SACtr type

# Composite type for SAC evenly-spaced time series data
@eval mutable struct SACtr
    $([:($(s)::SACFloat) for s in sac_float_hdr]...)
    $([:($(s)::SACInt) for s in sac_int_hdr]...)
    $([:($(s)::SACBool) for s in sac_bool_hdr]...)
    $([:($(s)::SACChar) for s in sac_char_hdr]...)
    # The time series, accessed with .t
    t::Array{SACFloat,1}
    function SACtr(delta_in::Number, npts_in::Integer, b_in=0.)
        delta_in > 0 || error("SACtr: delta must be positive")
        npts_in >= 0 || error("SACtr: npts must be 0 or larger")
        # Variables are by default undefined, or false for bools
        $([:($(s) = sac_rnull) for s in sac_float_hdr]...)
        $([:($(s) = sac_inull) for s in sac_int_hdr]...)
        $([:($(s) = false) for s in sac_bool_hdr]...)
        $([:($(s) = sac_cnull) for s in sac_char_hdr]...)
        # Variables which must be present
        npts = convert(SACInt, npts_in)
        delta = convert(SACFloat, delta_in)
        b = b_in
        e = b + (npts - 1)*delta
        t = zeros(npts)
        depmin = 0.
        depmax = 0.
        nvhdr = sac_ver_num
        iftype = 1
        idep = 5
        iztype = 9
        ievtyp = 5
        leven = true
        lovrok = true
        lcalda = true
        new($([:($(s)) for s in [sac_float_hdr; sac_int_hdr; sac_bool_hdr; sac_char_hdr]]...),
            t)
    end
end
@doc """
    SACtr(delta, npts, b=0.) -> ::SACtr

Construct a composite type holding an evenly-spaced SAC time-series trace, where the trace
is accessed through the field name `t`.  Supply the constant sampling interval `delta`
in seconds, and the number of points in the trace `t`.  Optionally, specify the trace
start time `b` in seconds.

    SACtr(v::AbstractVector, delta, b=0.) -> ::SACtr

Construct a `SACtr` by supplying an array `v`, sampling interval `delta` and optionally
the starting time.

    SACtr(d::Vector{UInt8}, file=""; swap=true, terse=false, check_npts=true) -> ::SACtr

Construct a SACtr from a raw array of bytes representing some data in SAC format.
If `swap` is false, then non-native-endian files are not converted.  If `terse` is
true, then warnings about swapping are not written.  If `check_npts` is false, then
parts of files are read without error.
""" SACtr


function SACtr(data::Vector{UInt8}, file=""; swap::Bool=true, terse::Bool=false,
        check_npts::Bool=true)
    len = sac_byte_len
    clen = 2*sac_byte_len
    # Determine endianness and act accordingly
    nvhdr = reinterpret(SACInt, data[sac_nvhdr_pos*len+1:(sac_nvhdr_pos+1)*len])[1]
    native = if nvhdr == sac_ver_num
        true
    elseif bswap(nvhdr) == sac_ver_num
        false
    else
        error("Array does not appear to be SAC data")
    end
    native && machine_is_little_endian && !swap &&
        error("Data are little-endian but `swap` is `false`.  Not attempting to swap bytes" *
            (file!="" ? " for file '$file'." : "."))
    native && machine_is_little_endian && !terse &&
        @info("Data are little-endian; byteswapping")
    byteswap(x) = native ? x : bswap(x)
    # Create an empty object
    npts_in_file = (length(data) - sac_header_len)÷len
    trace = SACtr(1, npts_in_file)

    ## Read header
    # Float part
    off = 0
    for (i, field) in enumerate(sac_float_hdr)
        setfield!(trace, field,
            byteswap(reinterpret(SACFloat, data[off+1:off+len])[1]))
        off += len
    end
    # Int part
    for (i, field) in enumerate(sac_int_hdr)
        setfield!(trace, field,
            byteswap(reinterpret(SACInt, data[off+1:off+len])[1]))
        off += len
    end
    # Boolean part
    for (i, field) in enumerate(sac_bool_hdr)
        setfield!(trace, field,
            1 == byteswap(reinterpret(SACInt, data[off+1:off+len])[1]))
        off += len
    end
    # Character part
    # kevnm header is double length, so treat separately
    trace.kstnm = String(strip(String(data[off+1:off+clen])))
    off += clen
    trace.kevnm = String(strip(String(data[off+1:off+2clen])))
    trace.kevnm == strip(rpad(sac_cnull, clen)^2) && (trace.kevnm = sac_cnull)
    off += 2clen
    for (i, field) in enumerate(sac_char_hdr)
        i <= 2 && continue
        setfield!(trace, field,
            String(strip(String(data[1+off:clen+off]))))
        off += clen
    end

    # Check length
    @assert off == sac_header_len
    check_npts && npts_in_file < trace.npts &&
        error("Number of points is not as expected: have $npts_in_file versus npts = " *
            "$(trace.npts) in header" * (file!="" ? " for file '$file'." : "."))

    # Now read in the trace
    trace.t .= reinterpret(SACFloat, data[(sac_header_len+1):end])
    native || (trace.t .= bswap.(trace.t))
    update_headers!(trace)
    any(isundefined.([trace.gcarc, trace.az, trace.baz])) && update_great_circle!(trace)
    trace
end

function SACtr(v::AbstractVector, delta::Real, b=0.0)
    s = SACtr(delta, length(v), b)
    s.t .= v
    update_headers!(s)
    s
end

"""
    getindex(A::AbstractArray{SACtr}, s::Symbol) -> Array{typeof(A[:].s)}
    A[:s] -> Array{typeof(A[:].s)}

Return an array of values containing the header with name `s` for the SACtr
traces.  This allows one to get all the headers values by doing A[:kstnm],
for example.
"""
getindex(A::AbstractArray{SACtr}, s::Symbol) = Array{typeof(getfield(A[1], s))}([getfield(a, s) for a in A])
getindex(t::SACtr, s::Symbol) = getfield(t, s) # Also define for single trace for consistency

"""
    setindex!(s::SACtr, value, s::Symbol)
    s[s] = value

Set the header with name `s` (a `Symbol`) to `value`.  E.g.:

    A[:kevnm] = "Blast 1"

    setindex!(A::AbstractArray{SACtr}, value, s::Symbol)
    A[s] = value

Set the header with name `s` for all the SACtr traces in the array `A`.  This
allows one to set all the headers at once for a set of traces by doing e.g.:

    A[:user0] = 1:length(A)

#### Setting trace values

When setting the values of the trace, use `s[:t] = k`.  This will automatically
update headers to reflect the trace and use broadcasting internally so that
the trace cannot be set to the wrong length.  Note however that a
broadcast assignment (`s[:t] .= k`) will not update the headers; this can be
done manually with `SAC.update_headers!(s)`.

    s[:t] = 1:s[:npts] # Sets to trace to be a line with constant slope
    s[:t] = sin.(2π*time(s)) # Sets the trace to a sine function
"""
function setindex!(t::SACtr, v, s::Symbol)
    if s == :t
        t.t .= v
        update_headers!(t)
    else
        setfield!(t, s, convert(typeof(getfield(t, s)), v))
        s in (:evlo, :evla, :stlo, :stla) && update_great_circle!(t)
    end
    t
end
setindex!(A::AbstractArray{SACtr}, V, s::Symbol) = setindex!.(A, V, s)

"""
    (==)(a::SACtr, b::SACtr) -> ::Bool

Return `true` if the traces `a` and `b` are equal (that is, have all fields the same),
and `false` otherwise.
"""
function (==)(a::SACtr, b::SACtr)
    for f in fieldnames(SACtr)
        if getfield(a, f) != getfield(b, f)
            return false
        end
    end
    true
end

"""
    copy(s::SACtr) -> t::SACtr

Return a copy of SAC trace `s`.
"""
function copy(s::SACtr)
    return SACtr(s.delta, s.depmin, s.depmax, s.scale, s.odelta, s.b, s.e, s.o, s.a, s.internal0,
        s.t0, s.t1, s.t2, s.t3, s.t4, s.t5, s.t6, s.t7, s.t8, s.t9, s.f,
        s.resp0, s.resp1, s.resp2, s.resp3, s.resp4, s.resp5, s.resp6, s.resp7, s.resp8, s.resp9,
        s.stla, s.stlo, s.stel, s.stdp, s.evla, s.evlo, s.evel, s.evdp, s.mag,
        s.user0, s.user1, s.user2, s.user3, s.user4, s.user5, s.user6, s.user7, s.user8, s.user9,
        s.dist, s.az, s.baz, s.gcarc, s.internal1, s.internal2, s.depmen, s.cmpaz, s.cmpinc,
        s.xminimum, s.xmaximum, s.yminimum, s.ymaximum,
        s.unused1, s.unused2, s.unused3, s.unused4, s.unused5, s.unused6, s.unused7,
        s.nzyear, s.nzjday, s.nzhour, s.nzmin, s.nzsec, s.nzmsec,
        s.nvhdr, s.norid, s.nevid, s.npts, s.internal3, s.nwfid, s.nxsize, s.nysize, s.unused8,
        s.iftype, s.idep, s.iztype, s.unused9, s.iinst, s.istreg, s.ievreg, s.ievtyp, s.iqual,
        s.isynth, s.imagtyp, s.imagsrc, s.unused10, s.unused11, s.unused12, s.unused13,
        s.unused14, s.unused15, s.unused16, s.unused17,
        s.leven, s.lpspol, s.lovrok, s.lcalda, s.unused18,
        s.kstnm, s.kevnm, s.khole, s.ko, s.ka, s.kt0, s.kt1, s.kt2, s.kt3, s.kt4, s.kt5, s.kt6, s.kt7,
        s.kt8, s.kt9, s.kf, s.kuser0, s.kuser1, s.kuser2, s.kcmpnm, s.knetwk, s.kdatrd, s.kinst,
        Base.copy(s.t))
end

@static if VERSION >= v"0.7-"
    Base.broadcastable(x::SACtr) = Ref(x)
end
