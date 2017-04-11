# Definition, constructors and other operations only to do with the SACtr type

# Composite type for SAC evenly-spaced time series data
@eval type SACtr
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

    SACtr(d::Vector{UInt8}, file=""; swap=true, terse=false, check_npts=true) -> ::SACtr

Construct a SACtr from a raw array of bytes representing some data in SAC format.
If `swap` is false, then non-native-endian files are not converted.  If `terse` is
true, then warnings about swapping are not written.  If `check_npts` is false, then
parts of files are read without error.
""" SACtr


@eval function SACtr(data::Vector{UInt8}, file=""; swap::Bool=true, terse::Bool=false,
        check_npts::Bool=true)
    const len = sac_byte_len
    const clen = 2*sac_byte_len
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
        info("Data are little-endian; byteswapping")
    byteswap(x) = native ? x : bswap(x)

    ## Read header
    # Float part
    $([:($s = byteswap(reinterpret(SACFloat, data[(($i-1)*len)+1:$i*len])[1])) for (s, i) in zip(sac_float_hdr, 1:length(sac_float_hdr))]...)
    off = length(sac_float_hdr)*len
    # Int part
    $([:($s = byteswap(reinterpret(SACInt, data[(($i-1)*len)+1+off:$i*len+off])[1])) for (s, i) in zip(sac_int_hdr, 1:length(sac_int_hdr))]...)
    off += length(sac_int_hdr)*len
    # Boolean part
    $([:($s = 0 != byteswap(reinterpret(SACInt, data[(($i-1)*len)+1+off:$i*len+off])[1])) for (s, i) in zip(sac_bool_hdr, 1:length(sac_bool_hdr))]...)
    off += length(sac_bool_hdr)*len
    # Character part
    # kevnm header is double length, so split into two then recombine
    char_sym_list = [sac_char_hdr[1]; :kevnm1; :kevnm2; sac_char_hdr[3:end]]
    $([:($s = ascii(String(reinterpret(UInt8, data[(($i-1)*clen)+1+off:$i*clen+off])))) for (s, i) in zip([sac_char_hdr[1]; :kevnm1; :kevnm2; sac_char_hdr[3:end]], 1:length(sac_char_hdr)+1)]...)
    kevnm = kevnm1 * kevnm2
    off += (length(sac_char_hdr) + 1)*clen

    # Check length
    @assert off == sac_header_len
    npts_in_file = (length(data) - sac_header_len)÷len
    check_npts && npts_in_file < npts &&
        error("Number of points is not as expected: have $npts_in_file versus npts = " *
            "$npts in header" * (file!="" ? " for file '$file'." : "."))

    # Create an empty object...
    trace = SACtr(delta, npts_in_file)
    # ...and fill the headers...
    $([:(trace.$s = $s) for s in sac_all_hdr]...)
    # ...then read in the trace
    trace.t .= reinterpret(SACFloat, data[(sac_header_len+1):end])
    native || (trace.t .= bswap.(trace.t))
    update_headers!(trace)
    any(isundefined.([trace.gcarc, trace.az, trace.baz])) && update_great_circle!(trace)
    trace
end

"""
    getindex(A::Array{SACtr}, s::Symbol) -> Array{typeof(A[:].s)}
    A[:s] -> Array{typeof(A[:].s)}

Return an array of values containing the header with name `s` for the SACtr
traces.  This allows one to get all the headers values by doing A[:kstnm],
for example.
"""
getindex(A::Array{SACtr}, s::Symbol) = Array{typeof(getfield(A[1], s))}([getfield(a, s) for a in A])
getindex(t::SACtr, s::Symbol) = getfield(t, s) # Also define for single trace for consistency

"""
    setindex!(A::Array{SACtr}, value, s::Symbol)
    A[:s] = value

Set the header with name `s` for all the SACtr traces in the array `A`.  This
allows one to set all the headers at once for a set of traces by doing e.g.:

    A[:kevnm] = "Blast 1"

or

    A[:user0] = 1:length(A)
"""
function setindex!(A::Array{SACtr}, V, s::Symbol)
    fieldtype = typeof(getfield(A[1], s))
    if length(A) == length(V)
        for (a, v) in zip(A, V)
            setfield!(a, s, convert(fieldtype, v))
        end
    elseif length(V) == 1
        for a in A
            setfield!(a, s, convert(fieldtype, V))
        end
    else
        error("Number of header values must be one or the number of traces")
    end
    s in (:evlo, :evla, :stlo, :stla) && update_great_circle!.(t)
    A
end
function setindex!(t::SACtr, v, s::Symbol)
    setfield!(t, s, convert(typeof(getfield(t, s)), v))
    s in (:evlo, :evla, :stlo, :stla) && update_great_circle!(t)
    t
end

"""
    (==)(a::SACtr, b::SACtr) -> ::Bool

Return `true` if the traces `a` and `b` are equal (that is, have all fields the same),
and `false` otherwise.
"""
function (==)(a::SACtr, b::SACtr)
    for f in fieldnames(a)
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
