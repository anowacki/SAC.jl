# Input/output routines

"""
    read(file; swap=true, terse=false) -> s::SACtr

Return the SAC trace as read from file `file` as a `SACtr` object.  If `swap` is false,
then auto-byteswapping is not performed and an error is returned if the file is not
of the assumed endianness.  Autoswapping is reported unless `terse` is `true`.
"""
function read(file; swap::Bool=true, terse::Bool=false, check_npts::Bool=true)
    data = open(file, "r") do f
        Base.read(f)
    end
    SACtr(data, file, swap=swap, terse=terse, check_npts=check_npts)
end

"""
    read_cut(file, b, e; swap=true, terse=false) -> s::SACtr
    read_cut(file, b_header, b_time, e_header, e_time; swap=true, terse=false) -> s::SACtr

Return the trace `s`, cut between beginning time `b` and end time `e` in seconds
relative to the zero time.

Optionally specify headers to which `b_time` and `e_time` are relative as
`b_header` and `e_header` respectively.  For example:

    s = read_cut(file, :a, -1, :f + 20)

If `b` is before the start of the trace, then the trace is read from the beginning;
similarly, if `e` is after the end of the trace, the trace is read up until the end.
In either case, a warning is issued unless `terse` is `true`.
"""
function read_cut(file, hb::Symbol, b::Real, he::Symbol, e::Real; swap=true, terse=false)
    b < e || throw(ArgumentError("Start cut must be before stop cut"))
    len = sac_byte_len
    clen = 2*sac_byte_len
    open(file, "r") do f
        header = Base.read(f, sac_header_len)
        # Determine endianness and act accordingly
        nvhdr = reinterpret(SACInt, header[sac_nvhdr_pos*len+1:(sac_nvhdr_pos+1)*len])[1]
        native = if nvhdr == sac_ver_num
            true
        elseif bswap(nvhdr) == sac_ver_num
            false
        else
            error("File does not appear to be SAC data (nvhdr = $nvhdr)")
        end
        native || swap || error("File $file is non-native-endian but `swap` is `false`")
        # Create new SACtr object, setting the trace to have one zero value
        append!(header, zeros(UInt8, 4))
        s = SACtr(header, check_npts=false)
        # Check args
        isundefined(s, hb) && error("Header $hb is undefined")
        isundefined(s, he) && error("Header $he is undefined")
        b < s.b && !terse &&
            @warn("Start cut $b is before trace start $(s.b); reading from trace start")
        e > s.e && !terse &&
            @warn("End cut $e is after trace end $(s.e); reading to trace end")
        # Now do cutting
        ib = round(Int, (b - getfield(s, hb))/s.delta) + 1
        ie = round(Int, (e - getfield(s, he))/s.delta) + 1
        npts = ie - ib + 1
        @assert ib >= 1 && ie <= s.npts
        s.b += (ib - 1)*s.delta
        s.npts = npts
        seek(f, sac_header_len + len*(ib - 1))
        s.t = reinterpret(SACFloat, Base.read(f, len*npts))
        native || (s.t .= bswap.(s.t))
        update_headers!(s)
    end
end
read_cut(file, b::Real, e::Real; kwargs...) = read_cut(file, :b, b, :b, e; kwargs...)

"""
    file_is_native_endian(file)

Return `true` if `file` is a native-endian (defined by a constant in the module to be
little-endian on this machine) SAC file, and `false` if not.

The heuristic is thus: native-endian files have bytes 305:308 which are
a representation of the value `6`.  `6` is the current magic SAC file version number,
hard-coded into the routine.
"""
function file_is_native_endian(file::String)
    nvhdr = try
        d = open(file, "r") do f
            seek(f, sac_nvhdr_pos*sac_byte_len)
            reinterpret(SACInt, Base.read(f, sac_byte_len))[1]
	    end
    catch err
        error("SAC.file_is_native_endian: Cannot open file '$file' for reading " *
            "or file is not the correct type (error $err)")
    end
    if nvhdr == sac_ver_num
        return true
    elseif bswap(nvhdr) == sac_ver_num
        return false
    else
        error("SAC.file_is_native_endian: File '$file' does not appear to be a " *
            "valid SAC file (nvhdr is $nvhdr)")
    end
end

"Write floats or integers either swapped or native-endian"
_write_swap(swap::Bool, F::IOStream, x) = Base.write(F, swap ? Base.bswap.(x) : x)
"Write a String as a SAC string of the correct length, padded with ' 's"
_write_string(F::IOStream, x::String, maxlen::Integer) =
    Base.write(F, x[1:min(length(x),maxlen)]*" "^(max(0,maxlen-length(x))))

@eval function write(s::SACtr, file; byteswap=sac_force_swap)
    open(file, "w") do f
        # Write header
        $([:(_write_swap(byteswap, f, s.$s)) for s in [sac_float_hdr; sac_int_hdr]]...)
        $([:(_write_swap(byteswap, f, SACInt(s.$s))) for s in sac_bool_hdr]...)
        # No byte-swapping needed for characters, but pad them to the correct length
        _write_string(f, s.kstnm, saccharlen)
        # Handle special case of double-length kevnm  header
        kevnm = isundefined(s.kevnm) ? "-12345  -12345  " : s.kevnm
        _write_string(f, kevnm, 2*saccharlen)
        $([:(_write_string(f, s.$s, saccharlen)) for s in sac_char_hdr[3:end]]...)
        # Trace
        _write_swap(byteswap, f, s.t)
    end
end
@doc """
    write(s::SACtr, file; byteswap)
    write(S::AbstractArray{SACtr}, files; byteswap)

Write a SAC trace `s` to `file`, or a set of traces `S` to a set of files `files`.
Set `byteswap` to `false` to force writing in native-endian format; set to `true`
to write bigendian files (MacSAC type).  The default is to write bigendian format.
""" write

function write(s::AbstractArray{SACtr}, file::Array{String}; args...)
    length(s) == length(file) || error("SAC.write: Arrays must be same length")
    for i = 1:length(s)
        write(s[i], file[i]; args...)
    end
    return
end

"""
    read_wild(pat, dir="."; echo=true) -> A, files

Read files matching globbing pattern `pat` from directory `dir`.
If `echo` is false, do not show which files are being read.

Returns an array of SACtr types `A`, and an array of file names `files`.
**NB:** If no files are matched, then empty `SACtr` and `String` arrays
are returned and a warning message printed.
"""
function read_wild(pat::String, dir::String="."; echo::Bool=true)
    files = Glob.glob(pat, dir)
    n = size(files, 1)
    if n == 0
        @info("SAC.read_wild: No files matching '$pat' in directory '$dir'")
        return SACtr[], String[]
    end
    A = Array{SACtr}(undef, n)
    for i = 1:n
        echo && @info("SAC.read: '$(files[i])'")
        A[i] = SAC.read(files[i]; terse=!echo)
    end
    return A, files
end

# List printing
function Base.show(io::IO, s::SACtr)
    out = "SAC.SACtr(delta=$(s[:delta])"
    for hdr in (:b, :npts, :kstnm, :gcarc, :az, :baz)
        if !isundefined(s[hdr])
            out *= ", " * string(hdr) * "=" * strip(string(s[hdr]))
        end
    end
    print(io, out * ")")
end

# Printing to the REPL, for instance
function Base.show(io::IO, ::MIME"text/plain", s::SACtr)
    hdr_string_len = maximum(length.(string.(sac_all_hdr)))
    print(io, "SAC.SACtr:")
    for hdr in sac_all_hdr
        !isundefined(s[hdr]) &&
            print(io, "\n", lpad(string(hdr), hdr_string_len, ' ') * ": ", s[hdr])
    end
end
