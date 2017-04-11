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
        _write_string(f, s.kevnm, 2*saccharlen)
        $([:(_write_string(f, s.$s, saccharlen)) for s in sac_char_hdr[3:end]]...)
        # Trace
        _write_swap(byteswap, f, s.t)
    end
end
@doc """
    write(s::SACtr, file; byteswap)
    write(S::Array{SACtr}, files; byteswap)

Write a SAC trace `s` to `file`, or a set of traces `S` to a set of files `files`.
Set `byteswap` to `false` to force writing in native-endian format; set to `true`
to write bigendian files (MacSAC type).  The default is to write bigendian format.
""" write

function write(s::Array{SACtr}, file::Array{String}; args...)
    length(s) == length(file) || error("SAC.write: Arrays must be same length")
    for i = 1:length(s)
        write(s[i], file[i]; args...)
    end
    return
end

"""
    read_wild(pat, dir=\"./\"; echo=true) -> A, files

Read files matching globbing pattern `pat` from directory `dir`.
If `echo` is false, do not show which files are being read.

Returns an array of SACtr types `A`, and an array of file names `files`.
"""
function read_wild(pat::String, dir::String="."; echo::Bool=true)
    # Return an array of SACtr types, and an array which gives the file path
    # for each trace.  Return nothing if there are no files.
    # Defaults to current directory.
    if !isdir(dir)
        info("SAC.read_wild: No directory '$dir'")
        return
    end
    files = Glob.glob(pat, dir)
    n = size(files, 1)
    if n == 0
        info("SAC.read_wild: No files matching '$pat' in directory '$dir'")
        return
    end
    A = Array(SACtr, n)
    for i = 1:n
        echo && info("SAC.read: '$(files[i])'")
        A[i] = SAC.read(files[i]; terse=!echo)
    end
    return A, files
end
