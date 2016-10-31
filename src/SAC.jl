"""
SAC.jl provides routines for dealing with SAC-formatted time series files,
including reading, writing, filtering, mean removal, rotating, and so on.
Sister library SACPlot.jl can be used for plotting.
"""
module SAC

__precompile__()

import DSP
import Glob
import Base: ==, copy, getindex, fft, setindex!, time, write
import Compat

const String = Compat.String

export
    SACtr,
    bandpass!,
    bp!,
    copy,
    cut!,
    envelope!,
    fft,
    highpass!,
    hp!,
    lp!,
    lowpass!,
    read_wild,
    rmean!,
    rotate_through!,
    rtrend!,
    taper!,
    time,
    tshift!,
    write


# SAC types
const SACFloat = Float32
const SACChar = String
const SACInt = Int32
const SACBool = Bool
# Constructors
sacfloat(x) = map(Float32, x)
sacint(x) = map(Int32, x)
const sacchar = ascii
sacbool(x) = x != 0
# Length of SAC floats and ints
const sac_byte_len = 4
# Length of SAC character headers (except kevnm, which is twice the length)
const saccharlen = 8
# SAC file version number
const sac_ver_num = SACInt(6)
# Whether this machine is big- or little endian.  SAC files are meant to be big-endian,
# so this determines whether a file is native-endian or not.
const machine_is_little_endian = true

# Convert a number into a SACChar
sacstring(x, maxlen=saccharlen) = sacchar(string(x)[1:minimum((length(string(x)),maxlen))]*" "^(maximum((0,maxlen-length(string(x))))))

# SAC unset values
const sac_rnull = -12345.
const sac_inull = -12345
const sac_cnull = "-12345"

# Default values for filtering
const sac_npoles = 2
const sac_passes = 1

# For SAC/BRIS, files are always big-endian, so set this to true to swap by default
const sac_force_swap = true

# Flag for verbosity
sac_verbose = true

# Lists of SAC headers as symbols
const sac_float_hdr = [:delta, :depmin, :depmax, :scale, :odelta, :b,
    :e, :o, :a, :internal0, :t0, :t1,
    :t2, :t3, :t4, :t5, :t6, :t7,
    :t8, :t9, :f, :resp0, :resp1, :resp2,
    :resp3, :resp4, :resp5, :resp6, :resp7, :resp8,
    :resp9, :stla, :stlo, :stel, :stdp, :evla,
    :evlo, :evel, :evdp, :mag, :user0, :user1,
    :user2, :user3, :user4, :user5, :user6, :user7,
    :user8, :user9, :dist, :az, :baz, :gcarc,
    :internal1, :internal2, :depmen, :cmpaz, :cmpinc, :xminimum,
    :xmaximum, :yminimum, :ymaximum, :unused1, :unused2, :unused3,
    :unused4, :unused5, :unused6, :unused7]
const sac_int_hdr = [:nzyear, :nzjday, :nzhour, :nzmin, :nzsec, :nzmsec,
    :nvhdr, :norid, :nevid, :npts, :internal3, :nwfid,
    :nxsize, :nysize, :unused8, :iftype, :idep, :iztype,
    :unused9, :iinst, :istreg, :ievreg, :ievtyp, :iqual,
    :isynth, :imagtyp, :imagsrc, :unused10, :unused11, :unused12,
    :unused13, :unused14, :unused15, :unused16, :unused17]
const sac_bool_hdr = [:leven, :lpspol, :lovrok, :lcalda, :unused18]
const sac_char_hdr = [:kstnm, :kevnm, :khole, :ko, :ka, :kt0,
    :kt1, :kt2, :kt3, :kt4, :kt5, :kt6,
    :kt7, :kt8, :kt9, :kf, :kuser0, :kuser1,
    :kuser2, :kcmpnm, :knetwk, :kdatrd, :kinst]
const sac_all_hdr = [sac_float_hdr; sac_int_hdr; sac_bool_hdr; sac_char_hdr]

# Where in the file the NVHDR value is
const sac_nvhdr_pos = length(sac_float_hdr) + find(sac_int_hdr .== :nvhdr)[1] - 1

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

    SACtr(d::Vector{UInt8}) -> ::SACtr

Construct a SACtr from a raw array of bytes representing some data in SAC format.
""" SACtr


@eval function SACtr(data::Vector{UInt8}; terse::Bool=false)
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

    # Create an empty object...
    trace = SACtr(delta, npts)
    # ...and fill the headers...
    $([:(trace.$s = $s) for s in sac_all_hdr]...)
    # ...then read in the trace
    for i = 1:npts
        trace.t[i] = byteswap(reinterpret(SACFloat, data[(i-1)*len+1+off:i*len+off])[1])
    end
    update_headers!(trace)
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
end
setindex!(t::SACtr, v, s::Symbol) = setfield(t, s, v)

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

@eval function read(file; swap::Bool=true, terse::Bool=false)
    const len = sac_byte_len
    const clen = 2*sac_byte_len
    # Determine endianness and act accordingly
    native = file_is_native_endian(file)
    native && machine_is_little_endian && !swap &&
        error("File '$file' is little-endian, but `swap` is `false`; use set to `true` to auto-byteswap")
    native && machine_is_little_endian && !terse &&
        info("File '$file' is little-endian; byteswapping")
    byteswap(x) = native ? x : bswap(x)

    ## Read data
    data = open(file, "r") do f
        Base.read(f)
    end

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

    # Create an empty object...
    trace = SACtr(delta, npts)
    # ...and fill the headers...
    $([:(trace.$s = $s) for s in sac_all_hdr]...)
    # ...then read in the trace
    for i = 1:npts
        trace.t[i] = byteswap(reinterpret(SACFloat, data[(i-1)*len+1+off:i*len+off])[1])
    end
    update_headers!(trace)
    trace
end
@doc """
    read(file; swap=true, terse=false) -> s::SACtr

Return the SAC trace as read from file `file` as a `SACtr` object.  If `swap` is false,
then auto-byteswapping is not performed and an error is returned if the file is not
of the assumed endianness.  Autoswapping is reported unless `terse` is `true`.
""" read

"""
    file_is_native_endian(file)

Return `true` if `file` is a native-endian (defined by a constant in the module to be
little-endian on this machine) SAC file, and `false` if not.

The heuristic is thus: native-endian files have bytes 305:308 which are
a representation of the value `6`.  `6` is the current magic SAC file version number,
hard-coded into the routine.
"""
function file_is_native_endian(file::Compat.String)
    nvhdr = try
        d = open(file, "r") do f
            Compat.readbytes(f, (sac_nvhdr_pos + 1)*sac_byte_len)
	    end
        reinterpret(SACInt, d[end-sac_byte_len+1:end])[1]
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


@eval function write(s::SACtr, file::String; byteswap=sac_force_swap)
    # Write a SAC composed type to file
    # Call with byteswap=true to write non-native-endian files
    f = open(file, "w")
    if byteswap
        w(F::IOStream, x) = Base.write(F, Base.bswap(x))
    else
        w(F::IOStream, x) = Base.write(F, x)
    end
    # Define routine to shorten/pad character headers
    w(F::IOStream, x::String, maxlen::Integer) =
        Base.write(F, x[1:minimum((length(x),maxlen))]*" "^(maximum((0,maxlen-length(x)))))
    # Write header
    $([:(w(f, s.$s)) for s in [sac_float_hdr; sac_int_hdr]]...)
    $([:(w(f, SACInt(s.$s))) for s in sac_bool_hdr]...)
    # No byte-swapping needed for characters, but pad them to the correct length
    w(f, s.kstnm, saccharlen)
    w(f, s.kevnm, 2*saccharlen)
    $([:(w(f, s.$s, saccharlen)) for s in sac_char_hdr[3:end]]...)
    # Trace
    for i = 1:s.npts
        w(f, s.t[i])
    end
    close(f)
end
@doc """
    write(s::SACtr, file; byteswap)
    write(S::Array{SACtr}, files; byteswap)

Write a SAC trace `s` to `file`, or a set of traces `S` to a set of files `files`.
Set `byteswap` to `true` to force writing in non-native-endian format; set to `false`
to write native-endian files.  The default is to write bigendian format.
""" -> write

function write(s::Array{SACtr}, file::Array{String}; args...)
    length(s) == length(file) || error("SAC.write: Arrays must be same length")
    for i = 1:length(s)
        write(s[i], file[i]; args...)
    end
    return
end


@doc """
`copy(s::SACtr) -> t::SACtr`

Return a copy of SAC trace `s`.
""" ->
function copy(s::SACtr)
    # Return a copy of a SAC trace
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

@doc """
    read_wild(pat, dir=\"./\"; echo=true) -> A, files

Read files matching globbing pattern `pat` from directory `dir`.
If `echo` is false, do not show which files are being read.

Returns an array of SACtr types `A`, and an array of file names `files`.
""" ->
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

@doc """
    sample() -> ::SACtr

Return some sample SAC data
""" ->
function sample()
    # Return some sample data, which is what you get when calling `fg seis' in SAC
    file = joinpath(dirname(@__FILE__()), "../data/seis.sac")
    return SAC.read(file)
end


@doc """
    cut!(s::SACtr, b::Number, e::Number)
    cut!(s::Array{SACtr}, b::Number, e::Number)

Cut a trace or array of traces `s` in memory between times `b` and `e`, relative
to the O marker.

    cut!(s::Array{SACtr}, a::Array, b::Array)

Cut the array of traces `s` between the times in arrays `b` and `e`, which must be
the same length as `s`.
""" ->
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
    s.b, s.e = b, e
    s.npts = ie - ib + 1
    update_headers!(s)
    return
end

# Array version of cut!
function cut!(a::Array{SACtr}, b::Number, e::Number)
    for s in a
        SAC.cut!(s, b, e)
    end
end

function cut!{B<:Real,E<:Real}(a::Array{SACtr}, b::Array{B}, e::Array{E})
    @assert length(a) == length(b) == length(e) "Arrays `a`, `b` and `e` must be the same length"
    for (s, beg, en) in zip(a, b, e)
        SAC.cut!(s, beg, en)
    end
end

@doc """
    fft(s::SACtr) -> f, S

Return the Fourier-transformed trace from `s` as `S`, with the frequencies
which correspond to each point in `f`.
""" ->
function fft(s::SACtr)
    # Return the fourier-transformed trace and the frequencies to go along with it
    N = round(Int, s.npts/2) + 1
    fmax = 1./(s.npts*s.delta)
    f = collect(1:N)*fmax
    S = Base.fft(s.t)[1:N]
    return f, S
end

function fft(a::Array{SACtr})
    # Return arrays containing f and S for an array of SACtr objects
    n = length(a)
    f, S = Array(Array, n), Array(Array, n)
    for i = 1:n
        f[i], S[i] = fft(a[i])
    end
    return f, S
end

@doc """
    rmean!(::SACtr)

Remove the mean in-place for a SAC trace.
""" ->
function rmean!(s::SACtr)
    # Remove the mean in-place
    s.t = s.t - mean(s.t)
    update_headers!(s)
    return
end

function rmean!(a::Array{SACtr})
    for s in a
        rmean!(s)
    end
end

@doc """
    rtrend!(::SACtr)

Remove the trend from a SAC trace in place.
""" ->
function rtrend!(s::SACtr)
    # Remove the trend in-place
    t = time(s)
    x0, x1 = linreg(t, s.t)
    s.t = s.t - (x0 + x1*t)
    update_headers!(s)
    return
end

function rtrend!(a::Array{SACtr})
    for s in a
        rtrend!(s)
    end
end

"""
    taper!(s::SACtr, width=0.05, form=:hanning)
    taper!(S::Array{SACtr}, width=0.05, form=:hanning)

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
    return
end
taper!(S::Array{SACtr}, width=0.05, form=:hamming::Symbol) = for s in S taper!(s, width, form) end

function update_headers!(s::SACtr)
    # Update headers which are automatically calculated from the trace
    s.depmax = maximum(s.t)
    s.depmin = minimum(s.t)
    s.depmen = mean(s.t)
    return
end

function update_headers!(a::Array{SACtr})
    for s in a
        update_headers!(s)
    end
end

@doc """
    time(::SACtr) -> t

Return a FloatRange `t` which contains the times for each sample of the SAC trace.
""" ->
time(s::SACtr) = s.b + (0:s.npts-1)*s.delta

@doc """
    bandpass!(::SACtr, c1, c2; ftype=\"butterworth\", npoles=2, passes=1)

Perform a bandpass filter on the SAC trace, between frequency corners `c1`
and `c2`.\n
Select type of filter with `ftype`: current options are: `butterworth`.
Set number of poles with `npoles`.\n
`passes` may be 1 (forward) or 2 (forward and reverse).
""" ->
function bandpass!(s::SACtr, c1::Number, c2::Number;
        ftype::String="butterworth", npoles::Integer=sac_npoles,
        passes::Integer=sac_passes)
                  # tranbw::Number=0.3, atten::Number=30)
    # Perform a bandpass on the trace, using either a Butterworth, Bessel or
    # Chebyshev (type 1 or 2) filter.
    # INPUT:
    #    s::SACtr     : SACtr composite type
    #    c1::Number   : Low corner / Hz
    #    c2::Number   : High corner / Hz
    # INPUT (OPTIONAL):
    #    type::String : Name of type.  Unambiguous short forms for the
    #                        following are acceptable:
    #                        [bu]tterworth [Default]
    #   npoles::Int       : Number of poles (1-10) [Default 2]
    #    passes::Int       : Number of passes (1-2) [Default 1]

    # Check arguments
    c1 >= c2 &&    error("SAC.bandpass: Upper corner must be larger than lower corner")
    response = DSP.Bandpass(c1, c2; fs=1./s.delta)
    prototype = get_filter_prototype(ftype, npoles)
    # Create apply the filter
    f = DSP.digitalfilter(response, prototype)
    apply_filter!(s, f, passes)
    return
end

function bandpass!(a::Array{SACtr}, c1, c2; ftype="butterworth", npoles=sac_npoles,
        passes=sac_passes)
    for s in a
        bandpass!(s, c1, c2; ftype=ftype, npoles=npoles, passes=passes)
    end
end

bp! = bandpass!

@doc """
    highpass!(::SACtr, c; ftype=\"butterworth\", npoles=2, passes=1)

Perform a highpass filter on the SAC trace, above frequency corner `c`.\n
Select type of filter with `ftype`: current options are: `butterworth`.
Set number of poles with `npoles`.\n
`passes` may be 1 (forward) or 2 (forward and reverse).
""" ->
function highpass!(s::SACtr, c::Number;
        ftype::String="butterworth", npoles::Integer=sac_npoles,
        passes::Integer=sac_passes)
    # Perform a highpass on the trace, in-place.
    response = DSP.Highpass(c; fs=1./s.delta)
    prototype = get_filter_prototype(ftype, npoles)
    f = DSP.digitalfilter(response, prototype)
    apply_filter!(s, f, passes)
    return
end

function highpass!(a::Array{SACtr}, c;
        ftype="butterworth", npoles=sac_npoles, passes=sac_passes)
    for s in a
        highpass!(s, c; ftype=ftype, npoles=npoles, passes=passes)
    end
end

hp! = highpass!

@doc """
    lowpass!(::SACtr, c; ftype=\"butterworth\", npoles=2, passes=1)

Perform a lowpass filter on the SAC trace, above frequency corner `c`.\n
Select type of filter with `ftype`: current options are: `butterworth`.
Set number of poles with `npoles`.\n
`passes` may be 1 (forward) or 2 (forward and reverse).
""" ->
function lowpass!(s::SACtr, c::Number;
        ftype::String="butterworth", npoles::Integer=sac_npoles,
        passes::Integer=sac_passes)
    # Perform a lowpass on the trace, in-place.
    response = DSP.Lowpass(c; fs=1./s.delta)
    prototype = get_filter_prototype(ftype, npoles)
    f = DSP.digitalfilter(response, prototype)
    apply_filter!(s, f, passes)
    return
end

function lowpass!(a::Array{SACtr}, c;
        ftype="butterworth", npoles=sac_npoles, passes=sac_passes)
    for s in a
        lowpass!(s, c; ftype=ftype, npoles=npoles, passes=passes)
    end
end

lp! = lowpass!

@doc """
    rotate_through!(::SACtr, ::SACtr, phi)

Given two SAC traces which are horizontal and orthgonal, rotate them clockwise
by `phi`° about the vertical axis.  This is a reference frame transformation
(passive rotation) and hence particle motion will appear to rotate
anti-clockwise.
""" ->
function rotate_through!(s1::SACtr, s2::SACtr, phi)
    # Rotate two orthogonal horizontal traces clockwise by ('through') phi (degrees).
    # This has the effect of changing the reference frame (passive rotation),
    # and hence the particle motion appears to rotate anti-clockwise.
    if mod(s2.cmpaz - s1.cmpaz, 180.) != 90.
        error("SAC.rotate_through!: traces must be orthogonal")
    elseif s1.npts != s2.npts
        error("SAC.rotate_through!: traces must be same length")
    elseif s1.delta != s2.delta
        error("SAC.rotate_through!: traces must have same delta")
    end
    phir = deg2rad(phi)
    R = [cos(phir) sin(phir);
        -sin(phir) cos(phir)]
    for i = 1:s1.npts
        (s1.t[i], s2.t[i]) = R*[s1.t[i]; s2.t[i]]
    end
    update_headers!(s1)
    s1.cmpaz = mod(s1.cmpaz + phi, 360.)
    s1.kcmpnm = sacstring(s1.cmpaz)
    update_headers!(s2)
    s2.cmpaz = mod(s2.cmpaz + phi, 360.)
    s2.kcmpnm = sacstring(s2.cmpaz)
    return
end

function rotate_through!(a::Array{SACtr}, phi)
    length(a)%2 != 0 && error("SAC.rotate_through!: Array of traces must be a multiple of two long")
    for i = 1:length(a)/2
        rotate_through!(a[2*i - 1], a[2*i], phi)
    end
end

@doc """
    tshift!(::SACtr, tshift; wrap=true)

Shift a SAC trace backward in time by `t` seconds.

If `wrap` true (default), then points which move out the back of the trace
are added to the front (and vice versa).  Setting it to false instead pads the
trace with zeroes.
""" ->
function tshift!(s::SACtr, tshift::Number; wrap=true)
    # Shift a trace backward in time by t seconds, wrapping around by default,
    # or optionally zeroing the front/endmost samples if pad=false
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
    return
end

@doc """
    envelope!(::SACtr)

Find the envelope of a SAC trace
""" ->
function envelope!(a::Array{SACtr})
    for s in a
        s.t = abs(DSP.hilbert(s.t))
    end
    update_headers!(a)
end
envelope!(s::SACtr) = envelope!([s])

@doc """
    multiply!(::SACtr, value)

Multiply the values in a SAC trace by `value`
""" ->
function multiply!(a::Array{SACtr}, val)
    for s in a s.t[:] = s.t[:]*val end
    update_headers!(a)
end
multiply!(s::SACtr, val) = multiply!([s], val)
mul! = multiply!

@doc """
    add!(::SACtr, value)

Add a constant value to a SAC trace
""" ->
function add!(a::Array{SACtr}, val)
    for s in a s.t[:] = s.t[:] + val end
    update_headers!(a)
end
add!(s::SACtr, val) = add!([s], val)

@doc """
    divide!(::SACtr, value)

Divide the values in a SAC trace by `value`
""" ->
function divide!(a::Array{SACtr}, value)
    value != 0. || error("SAC.divide!: Cannot divide by 0")
    multiply!(a, 1./value)
end
divide!(s::SACtr, value) = divide!([s], value)
div! = divide!


function apply_filter!(s::SACtr, f, passes::Integer)
        passes < 1 || passes > 2 && error("SAC.apply_filter!: Number of passes must be 1 or 2")
    if passes == 1
        DSP.filt!(s.t, f, s.t)
    elseif passes == 2
        s.t = DSP.filtfilt(f, s.t)
    else
        error("SAC.apply_filter!: passes must be 1 or 2")
    end
    update_headers!(s)
    return
end

function get_filter_prototype(ftype::String, npoles::Integer)
    # Return a filter prototype for use with filtering
    # INPUT:
    #    type::String : Name of type.  Unambiguous short forms for the
    #                        following are acceptable:
    #                        [bu]tterworth [Default]
    #                       [be]ssel
    #                       chebyshev1 [c1]
    #                       chebyshev2 [c2]
    #   npoles::Int       : Number of poles (1-10) [Default 2]
    npoles < 1 || npoles > 10 &&
        error("SAC.get_filter_prototype: npoles must be in range 1 - 10")
    length(ftype) < 2 && error("SAC.get_filter_prototype: ftype must " *
        "be at least two characters long")
    if lowercase(ftype[1:2]) == "bu"
        prototype = DSP.Butterworth(npoles)
    elseif lowercase(ftype[1:2]) == "be"
        error("SAC.get_filter_prototype: Bessel filter type not implemented yet")
    elseif lowercase(ftype[1:2]) == "c1" || lowercase(ftype) == "chebyshev1"
        error("SAC.get_filter_prototype: Chebyshev1 filter type not implemented yet")
        prototype = DSP.Chebyshev1(npoles)
    elseif lowercase(ftype[1:2]) == "c2" || lowercase(ftype) == "chebyshev2"
        error("SAC.get_filter_prototype: Chebyshev2 filter type not implemented yet")
        prototype = DSP.Chebyshev2(npoles)
    else
        error("SAC.get_filter_prototype: unrecognised filter type '$ftype'")
    end
    return prototype
end

end # module SAC
