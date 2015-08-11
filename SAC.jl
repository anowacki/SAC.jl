module SAC
# Contains routines for dealing with SAC-formatted files

import DSP
import Glob
import Base.read, Base.write, Base.copy, Base.fft, Base.time
if VERSION < v"0.4"
	import Docile.@doc
end

export SACtr,
	bandpass!,
	bp!,
	copy,
	cut!,
	fft,
	highpass!,
	hp!,
	lp!,
	lowpass!,
	rmean!,
	rotate_through!,
	rtrend!,
	time,
	tshift!,
	write


# SAC types
const SACFloat = Float32
const SACChar = ASCIIString
const SACInt = Int32
const SACBool = Bool
# Constructors
sacfloat = float32
sacint = int32
sacchar = ascii
sacbool = bool
# Length of SAC character headers (except kevnm, which is twice the length)
const saccharlen = 8
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

# Composite type for SAC evenly-spaced time series data
type SACtr
	# Header, floating point.  These are Float32, but we use
	# FloatingPoint internally for ease of use.
	delta::SACFloat
	depmin::SACFloat
	depmax::SACFloat
	scale::SACFloat
	odelta::SACFloat
	b::SACFloat
	e::SACFloat
	o::SACFloat
	a::SACFloat
	internal0::SACFloat
	t0::SACFloat
	t1::SACFloat
	t2::SACFloat
	t3::SACFloat
	t4::SACFloat
	t5::SACFloat
	t6::SACFloat
	t7::SACFloat
	t8::SACFloat
	t9::SACFloat
	f::SACFloat
	resp0::SACFloat
	resp1::SACFloat
	resp2::SACFloat
	resp3::SACFloat
	resp4::SACFloat
	resp5::SACFloat
	resp6::SACFloat
	resp7::SACFloat
	resp8::SACFloat
	resp9::SACFloat
	stla::SACFloat
	stlo::SACFloat
	stel::SACFloat
	stdp::SACFloat
	evla::SACFloat
	evlo::SACFloat
	evel::SACFloat
	evdp::SACFloat
	mag::SACFloat
	user0::SACFloat
	user1::SACFloat
	user2::SACFloat
	user3::SACFloat
	user4::SACFloat
	user5::SACFloat
	user6::SACFloat
	user7::SACFloat
	user8::SACFloat
	user9::SACFloat
	dist::SACFloat
	az::SACFloat
	baz::SACFloat
	gcarc::SACFloat
	internal1::SACFloat
	internal2::SACFloat
	depmen::SACFloat
	cmpaz::SACFloat
	cmpinc::SACFloat
	xminimum::SACFloat
	xmaximum::SACFloat
	yminimum::SACFloat
	ymaximum::SACFloat
	unused1::SACFloat
	unused2::SACFloat
	unused3::SACFloat
	unused4::SACFloat
	unused5::SACFloat
	unused6::SACFloat
	unused7::SACFloat
	# Integer parts.  These are Int32; again, we use Integer
    nzyear::SACInt
	nzjday::SACInt
	nzhour::SACInt
	nzmin::SACInt
	nzsec::SACInt
	nzmsec::SACInt
	nvhdr::SACInt
	norid::SACInt
	nevid::SACInt
	npts::SACInt
	internal3::SACInt
	nwfid::SACInt
	nxsize::SACInt
	nysize::SACInt
	unused8::SACInt
	iftype::SACInt
	idep::SACInt
	iztype::SACInt
	unused9::SACInt
	iinst::SACInt
	istreg::SACInt
	ievreg::SACInt
	ievtyp::SACInt
	iqual::SACInt
	isynth::SACInt
	imagtyp::SACInt
	imagsrc::SACInt
	unused10::SACInt
	unused11::SACInt
	unused12::SACInt
	unused13::SACInt
	unused14::SACInt
	unused15::SACInt
	unused16::SACInt
	unused17::SACInt
	# Logical part: boolean
	leven::SACBool
	lpspol::SACBool
	lovrok::SACBool
	lcalda::SACBool
	unused18::SACBool
	# Character parts
	kstnm::SACChar
	kevnm::SACChar
	khole::SACChar
	ko::SACChar
	ka::SACChar
	kt0::SACChar
	kt1::SACChar
	kt2::SACChar
	kt3::SACChar
	kt4::SACChar
	kt5::SACChar
	kt6::SACChar
	kt7::SACChar
	kt8::SACChar
	kt9::SACChar
	kf::SACChar
	kuser0::SACChar
	kuser1::SACChar
	kuser2::SACChar
	kcmpnm::SACChar
	knetwk::SACChar
	kdatrd::SACChar
	kinst::SACChar
	# The time series
	t::Array{SACFloat,1}
end

function SACtr(delta, npts)
	# Return a new SAC type with field filled in
	# Define variables which must be present
	b = 0.
	e = (npts - 1)*delta
	t = zeros(npts)
	depmin = 0.
	depmax = 0.
	nvhdr = 6
	iftype = 1
	idep = 5
	iztype = 9
	ievtyp = 5
	leven = true
	lpspol = false
	lovrok = true
	lcalda = true
	unused18 = false
	# Other variables are by default undefined
	scale = odelta = o = a = internal0 =
	        t0 = t1 = t2 = t3 = t4 = t5 = t6 = t7 = t8 = t9 = f =
	        resp0 = resp1 = resp2 = resp3 = resp4 = resp5 = resp6 = resp7 = resp8 = resp9 =
	        stla = stlo = stel = stdp = evla = evlo = evel = evdp = mag =
	        user0 = user1 = user2 = user3 = user4 = user5 = user6 = user7 = user8 = user9 =
	        dist = az = baz = gcarc = internal1 = internal2 = depmen = cmpaz = cmpinc =
	        xminimum = xmaximum = yminimum = ymaximum =
	        unused1 = unused2 = unused3 = unused4 = unused5 = unused6 = unused7 =
		sac_rnull
	nzyear = nzjday = nzhour = nzmin = nzsec = nzmsec =
			nvhdr = norid = nevid = npts = internal3 = nwfid = nxsize = nysize = unused8 =
			iftype = idep = iztype = unused9 = iinst = istreg = ievreg = ievtyp = iqual =
			isynth = imagtyp = imagsrc = unused10 = unused11 = unused12 = unused13 =
			unused14 = unused15 = unused16 = unused17 =
		sac_inull
	kevnm = kstnm = khole = ko = ka = kt0 = kt1 = kt2 = kt3 = kt4 = kt5 = kt6 = kt7 =
			kt8 = kt9 = kf = kuser0 = kuser1 = kuser2 = kcmpnm = knetwk = kdatrd = kinst =
		sac_cnull


	return SACtr(delta, depmin, depmax, scale, odelta, b, e, o, a, internal0,
        t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, f,
		resp0, resp1, resp2, resp3, resp4, resp5, resp6, resp7, resp8, resp9,
		stla, stlo, stel, stdp, evla, evlo, evel, evdp, mag,
		user0, user1, user2, user3, user4, user5, user6, user7, user8, user9,
		dist, az, baz, gcarc, internal1, internal2, depmen, cmpaz, cmpinc,
		xminimum, xmaximum, yminimum, ymaximum,
		unused1, unused2, unused3, unused4, unused5, unused6, unused7,
		# Integers
		nzyear, nzjday, nzhour, nzmin, nzsec, nzmsec,
		nvhdr, norid, nevid, npts, internal3, nwfid, nxsize, nysize, unused8,
		iftype, idep, iztype, unused9, iinst, istreg, ievreg, ievtyp, iqual,
		isynth, imagtyp, imagsrc, unused10, unused11, unused12, unused13,
		unused14, unused15, unused16, unused17,
		# Logical
		leven, lpspol, lovrok, lcalda, unused18,
		# Character
		kstnm, kevnm, khole, ko, ka, kt0, kt1, kt2, kt3, kt4, kt5, kt6, kt7,
		kt8, kt9, kf, kuser0, kuser1, kuser2, kcmpnm, knetwk, kdatrd, kinst,
		# Trace
		t)
end

@doc """
read(file; byteswap=true, terse=false) -> s::SACtr

   Read a SAC trace from \"file\".  If \"byteswap\" is true,
   enforce swapping between native-endian and non-native-endian.
   If \"terse\" is true, do not print warnings.
   Returns a SACtr type \"s\".
""" ->
function read(file; byteswap="auto", terse::Bool=false)
	# Read a binary SAC evenly-spaced time series file from disk.  Try to byteswap
	# if possible; prevent with bswap="no", force with bswap="yes"
	const len = 4
	const clen = 8
	data = readbytes(open(file, "r"))
	# Check the header version is as expected given our byte-swapping choice
	# (nvhdr is the 77th len-byte record)
	nvhdr = reinterpret(Int32, data[76*len+1:77*len])[1]
	if lowercase(byteswap[1:1]) == "a"
		if nvhdr != 6
			swap(x) = bswap(x)
			if bswap(nvhdr) != 6
				error("SAC.read: Header version '$nvhdr' for file '$file' " *
				      "is out of date for big- and little-endianness.")
			end
			!terse && info("SAC.read: Auto-byteswapping file '$file'...")
		else
			swap(x) = x
		end
	elseif lowercase(byteswap[1:1]) == "n"
		swap(x) = x
		if nvhdr != 6
			error("SAC.read: Header version '$nvhdr' for file '$file' " *
			      "is out of date (file byte-swapped?)")
		end
	elseif lowercase(byteswap[1:1]) == "y"
		swap(x) = bswap(x)
		if bswap(nvhdr) != 6
			error("SAC.read: Header version '$nvhdr' for file '$file' " *
			      "is out of date for non-native endianness (file native endian?)")
		end
	else
		error("SAC.read: byteswap must be [y]es, [n]o or [a]uto.")
	end
	# Read in header
	i = 1
	delta = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	depmin = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	depmax = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	scale = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	odelta = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	b = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	e = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	o = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	a = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	internal0 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	t0 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	t1 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	t2 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	t3 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	t4 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	t5 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	t6 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	t7 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	t8 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	t9 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	f = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	resp0 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	resp1 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	resp2 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	resp3 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	resp4 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	resp5 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	resp6 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	resp7 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	resp8 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	resp9 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	stla = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	stlo = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	stel = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	stdp = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	evla = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	evlo = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	evel = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	evdp = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	mag = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	user0 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	user1 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	user2 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	user3 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	user4 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	user5 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	user6 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	user7 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	user8 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	user9 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	dist = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	az = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	baz = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	gcarc = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	internal1 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	internal2 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	depmen = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	cmpaz = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	cmpinc = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	xminimum = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	xmaximum = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	yminimum = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	ymaximum = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	unused1 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	unused2 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	unused3 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	unused4 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	unused5 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	unused6 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
	unused7 = swap(reinterpret(Float32, data[i:i+len-1])[1]);  i+=len
    nzyear = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	nzjday = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	nzhour = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	nzmin = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	nzsec = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	nzmsec = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	nvhdr = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	norid = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	nevid = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	npts = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	internal3 = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	nwfid = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	nxsize = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	nysize = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	unused8 = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	iftype = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	idep = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	iztype = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	unused9 = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	iinst = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	istreg = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	ievreg = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	ievtyp = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	iqual = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	isynth = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	imagtyp = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	imagsrc = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	unused10 = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	unused11 = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	unused12 = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	unused13 = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	unused14 = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	unused15 = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	unused16 = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	unused17 = swap(reinterpret(Int32, data[i:i+len-1])[1]);  i+=len
	leven = bool(swap(reinterpret(Int32, data[i:i+len-1])[1]));  i+=len
	lpspol = bool(swap(reinterpret(Int32, data[i:i+len-1])[1]));  i+=len
	lovrok = bool(swap(reinterpret(Int32, data[i:i+len-1])[1]));  i+=len
	lcalda = bool(swap(reinterpret(Int32, data[i:i+len-1])[1]));  i+=len
	unused18 = bool(swap(reinterpret(Int32, data[i:i+len-1])[1]));  i+=len
	# 8 characters long
	kstnm = ascii(char(data[i:i+clen-1]));  i+=clen
	# 16 characters long
	kevnm = ascii(char(data[i:i+2*clen-1]));  i+=2*clen
	# 8 characters long again
	khole = ascii(char(data[i:i+clen-1]));  i+=clen
	ko = ascii(char(data[i:i+clen-1]));  i+=clen
	ka = ascii(char(data[i:i+clen-1]));  i+=clen
	kt0 = ascii(char(data[i:i+clen-1]));  i+=clen
	kt1 = ascii(char(data[i:i+clen-1]));  i+=clen
	kt2 = ascii(char(data[i:i+clen-1]));  i+=clen
	kt3 = ascii(char(data[i:i+clen-1]));  i+=clen
	kt4 = ascii(char(data[i:i+clen-1]));  i+=clen
	kt5 = ascii(char(data[i:i+clen-1]));  i+=clen
	kt6 = ascii(char(data[i:i+clen-1]));  i+=clen
	kt7 = ascii(char(data[i:i+clen-1]));  i+=clen
	kt8 = ascii(char(data[i:i+clen-1]));  i+=clen
	kt9 = ascii(char(data[i:i+clen-1]));  i+=clen
	kf = ascii(char(data[i:i+clen-1]));  i+=clen
	kuser0 = ascii(char(data[i:i+clen-1]));  i+=clen
	kuser1 = ascii(char(data[i:i+clen-1]));  i+=clen
	kuser2 = ascii(char(data[i:i+clen-1]));  i+=clen
	kcmpnm = ascii(char(data[i:i+clen-1]));  i+=clen
	knetwk = ascii(char(data[i:i+clen-1]));  i+=clen
	kdatrd = ascii(char(data[i:i+clen-1]));  i+=clen
	kinst = ascii(char(data[i:i+clen-1]));  i+=clen
	# Trace
	t = zeros(npts)
	for i = 1:npts
		j = int(633 + (i-1)*len)
		t[i] = swap(reinterpret(Float32, data[j:j+len-1])[1])
	end

	return SACtr(delta, depmin, depmax, scale, odelta, b, e, o, a, internal0,
        t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, f,
		resp0, resp1, resp2, resp3, resp4, resp5, resp6, resp7, resp8, resp9,
		stla, stlo, stel, stdp, evla, evlo, evel, evdp, mag,
		user0, user1, user2, user3, user4, user5, user6, user7, user8, user9,
		dist, az, baz, gcarc, internal1, internal2, depmen, cmpaz, cmpinc,
		xminimum, xmaximum, yminimum, ymaximum,
		unused1, unused2, unused3, unused4, unused5, unused6, unused7,
		# Integers
		nzyear, nzjday, nzhour, nzmin, nzsec, nzmsec,
		nvhdr, norid, nevid, npts, internal3, nwfid, nxsize, nysize, unused8,
		iftype, idep, iztype, unused9, iinst, istreg, ievreg, ievtyp, iqual,
		isynth, imagtyp, imagsrc, unused10, unused11, unused12, unused13,
		unused14, unused15, unused16, unused17,
		# Logical
		leven, lpspol, lovrok, lcalda, unused18,
		# Character
		kstnm, kevnm, khole, ko, ka, kt0, kt1, kt2, kt3, kt4, kt5, kt6, kt7,
		kt8, kt9, kf, kuser0, kuser1, kuser2, kcmpnm, knetwk, kdatrd, kinst,
		# Trace
		t)
end

@doc """
write(s::SACtr, file; byteswap)

   Write a SAC trace \"s\" to \"file\".  If \"byteswap\" is true, write
   in non-native-endian format.
""" ->
function write(s::SACtr, file::ASCIIString; byteswap=sac_force_swap)
	# Write a SAC composed type to file
	# Call with byteswap=true to write non-native-endian files
	f = open(file, "w")
	if byteswap
		w(F::IOStream, x) = Base.write(F, Base.bswap(x))
	else
		w(F::IOStream, x) = Base.write(F, x)
	end
	# Define routine to shorten/pad character headers
	w(F::IOStream, x::ASCIIString, maxlen::Integer) =
		Base.write(F, x[1:minimum((length(x),maxlen))]*" "^(maximum((0,maxlen-length(x)))))
	# Write header
	w(f, s.delta)
	w(f, s.depmin)
	w(f, s.depmax)
	w(f, s.scale)
	w(f, s.odelta)
	w(f, s.b)
	w(f, s.e)
	w(f, s.o)
	w(f, s.a)
	w(f, s.internal0)
	w(f, s.t0)
	w(f, s.t1)
	w(f, s.t2)
	w(f, s.t3)
	w(f, s.t4)
	w(f, s.t5)
	w(f, s.t6)
	w(f, s.t7)
	w(f, s.t8)
	w(f, s.t9)
	w(f, s.f)
	w(f, s.resp0)
	w(f, s.resp1)
	w(f, s.resp2)
	w(f, s.resp3)
	w(f, s.resp4)
	w(f, s.resp5)
	w(f, s.resp6)
	w(f, s.resp7)
	w(f, s.resp8)
	w(f, s.resp9)
	w(f, s.stla)
	w(f, s.stlo)
	w(f, s.stel)
	w(f, s.stdp)
	w(f, s.evla)
	w(f, s.evlo)
	w(f, s.evel)
	w(f, s.evdp)
	w(f, s.mag)
	w(f, s.user0)
	w(f, s.user1)
	w(f, s.user2)
	w(f, s.user3)
	w(f, s.user4)
	w(f, s.user5)
	w(f, s.user6)
	w(f, s.user7)
	w(f, s.user8)
	w(f, s.user9)
	w(f, s.dist)
	w(f, s.az)
	w(f, s.baz)
	w(f, s.gcarc)
	w(f, s.internal1)
	w(f, s.internal2)
	w(f, s.depmen)
	w(f, s.cmpaz)
	w(f, s.cmpinc)
	w(f, s.xminimum)
	w(f, s.xmaximum)
	w(f, s.yminimum)
	w(f, s.ymaximum)
	w(f, s.unused1)
	w(f, s.unused2)
	w(f, s.unused3)
	w(f, s.unused4)
	w(f, s.unused5)
	w(f, s.unused6)
	w(f, s.unused7)
    w(f, s.nzyear)
	w(f, s.nzjday)
	w(f, s.nzhour)
	w(f, s.nzmin)
	w(f, s.nzsec)
	w(f, s.nzmsec)
	w(f, s.nvhdr)
	w(f, s.norid)
	w(f, s.nevid)
	w(f, s.npts)
	w(f, s.internal3)
	w(f, s.nwfid)
	w(f, s.nxsize)
	w(f, s.nysize)
	w(f, s.unused8)
	w(f, s.iftype)
	w(f, s.idep)
	w(f, s.iztype)
	w(f, s.unused9)
	w(f, s.iinst)
	w(f, s.istreg)
	w(f, s.ievreg)
	w(f, s.ievtyp)
	w(f, s.iqual)
	w(f, s.isynth)
	w(f, s.imagtyp)
	w(f, s.imagsrc)
	w(f, s.unused10)
	w(f, s.unused11)
	w(f, s.unused12)
	w(f, s.unused13)
	w(f, s.unused14)
	w(f, s.unused15)
	w(f, s.unused16)
	w(f, s.unused17)
	w(f, int32(s.leven))
	w(f, int32(s.lpspol))
	w(f, int32(s.lovrok))
	w(f, int32(s.lcalda))
	w(f, int32(s.unused18))
	# No byte-swapping needed for characters, but pad them to the correct length
	w(f, s.kstnm, saccharlen)
	w(f, s.kevnm, 2*saccharlen)
	w(f, s.khole, saccharlen)
	w(f, s.ko, saccharlen)
	w(f, s.ka, saccharlen)
	w(f, s.kt0, saccharlen)
	w(f, s.kt1, saccharlen)
	w(f, s.kt2, saccharlen)
	w(f, s.kt3, saccharlen)
	w(f, s.kt4, saccharlen)
	w(f, s.kt5, saccharlen)
	w(f, s.kt6, saccharlen)
	w(f, s.kt7, saccharlen)
	w(f, s.kt8, saccharlen)
	w(f, s.kt9, saccharlen)
	w(f, s.kf, saccharlen)
	w(f, s.kuser0, saccharlen)
	w(f, s.kuser1, saccharlen)
	w(f, s.kuser2, saccharlen)
	w(f, s.kcmpnm, saccharlen)
	w(f, s.knetwk, saccharlen)
	w(f, s.kdatrd, saccharlen)
	w(f, s.kinst, saccharlen)
	# Trace
	for i = 1:s.npts
		w(f, s.t[i])
	end
	close(f)
end

@doc """
copy(s::SACtr) -> t::SACtr

   Return a copy of SAC trace <s>.
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

   Read files matching globbing pattern \"pat\" from directory \"dir\".
   If \"echo\" is false, do not show which files are being read.

   Returns an array of SACtr types \"A\", and an
""" ->
function read_wild(pat::ASCIIString, dir::ASCIIString="."; echo::Bool=true)
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
	file = dirname(@__FILE__()) * "/data/seis.sac"
	return read(file)
end


@doc """
cut!(s::SACtr, b::Number, e::Number)

Cut a trace \"s\" in memory between times \"b\" and \"e\", relative to O marker
""" ->
function cut!(s::SACtr, b::Number, e::Number)
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
	ib = int((b - s.b)/s.delta) + 1
	ie = s.npts - int((s.e - e)/s.delta)
	s.t = s.t[ib:ie]
	s.b, s.e = b, e
	s.npts = ie - ib + 1
	update_headers!(s)
	return
end

# Array version of cut!
function cut!(a::Array{SACtr}, b, e)
	for s in a
		SAC.cut!(s, b, e)
	end
end

@doc """
fft(s::SACtr) -> f, S

   Return the Fourier-transformed trace from \"s\" as \"S\", with the frequencies
   which correspond to each point in \"f\".
""" ->
function fft(s::SACtr)
	# Return the fourier-transformed trace and the frequencies to go along with it
	N = int(s.npts/2) + 1
	fmax = 1./(s.npts*s.delta)
	f = [1:N]*fmax
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

   Return an array \"t\" which contains the time for each sample of the SAC trace.
""" ->
function time(s::SACtr)
	# Return an array containing the times for each sample
	return s.b + [0:s.npts-1;]*s.delta
end

@doc """
bandpass!(::SACtr, c1, c2; ftype=\"butterworth\", npoles=2, passes=1)

   Perform a bandpass filter on the SAC trace, between frequency corners \"c1\"
   and \"c2\".\n
   Select type of filter with \"ftype\": current options are: \"butterworth\".
   Set number of poles with \"npoles\".\n
   \"passes\" may be 1 (forward) or 2 (forward and reverse).
""" ->
function bandpass!(s::SACtr, c1::Number, c2::Number;
		ftype::ASCIIString="butterworth", npoles::Integer=sac_npoles,
		passes::Integer=sac_passes)
				  # tranbw::Number=0.3, atten::Number=30)
	# Perform a bandpass on the trace, using either a Butterworth, Bessel or
	# Chebyshev (type 1 or 2) filter.
	# INPUT:
	#	s::SACtr     : SACtr composite type
	#	c1::Number   : Low corner / Hz
	#	c2::Number   : High corner / Hz
	# INPUT (OPTIONAL):
	#	type::ASCIIString : Name of type.  Unambiguous short forms for the
	#	                    following are acceptable:
	#	                    [bu]tterworth [Default]
	#   npoles::Int       : Number of poles (1-10) [Default 2]
	#	passes::Int       : Number of passes (1-2) [Default 1]

	# Check arguments
	c1 >= c2 &&	error("SAC.bandpass: Upper corner must be larger than lower corner")
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

   Perform a highpass filter on the SAC trace, above frequency corner \"c\".\n
   Select type of filter with \"ftype\": current options are: \"butterworth\".
   Set number of poles with \"npoles\".\n
   \"passes\" may be 1 (forward) or 2 (forward and reverse).
""" ->
function highpass!(s::SACtr, c::Number;
		ftype::ASCIIString="butterworth", npoles::Integer=sac_npoles,
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

   Perform a lowpass filter on the SAC trace, above frequency corner \"c\".\n
   Select type of filter with \"ftype\": current options are: \"butterworth\".
   Set number of poles with \"npoles\".\n
   \"passes\" may be 1 (forward) or 2 (forward and reverse).
""" ->
function lowpass!(s::SACtr, c::Number;
		ftype::ASCIIString="butterworth", npoles::Integer=sac_npoles,
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
   by \"phi\"° about the vertical axis.  This is a reference frame transformation
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

@doc """
tshift!(::SACtr, tshift; warp=true)

   Shift a SAC trace backward in time by \"t\" seconds.

   If \"wrap\" true (default), then points which move out the back of the trace
   are added to the front (and vice versa).  Setting it to false instead pads the
   trace with zeroes.
""" ->
function tshift!(s::SACtr, tshift::Number; wrap=true)
	# Shift a trace backward in time by t seconds, wrapping around by default,
	# or optionally zeroing the front/endmost samples if pad=false
	n = int(tshift/s.delta)
	if n == 0
		sac_verbose && info("SAC.tshift!: t ($tshift) is less than delta ($(s.delta)) so no shift applied")
		return
	end
	s.t = circshift(s.t, n)
	if !wrap
		n > 0 ? s.t[1:n] = 0. : s.t[end+n+1:end] = 0.
	end
	return
end

function apply_filter!(s::SACtr, f, passes::Integer)
		passes < 1 || passes > 2 && error("SAC.apply_filter!: Number of passes must be 1 or 2")
	if passes == 1
		DSP.filt!(s.t, f, s.t)
	elseif passes == 2
		s.t = DSP.filtfilt(f, s.t)
	else
		error("SAC.apply_filter!: passes must be 1 or 2")
	end
	return
end

function get_filter_prototype(ftype::ASCIIString, npoles::Integer)
	# Return a filter prototype for use with filtering
	# INPUT:
	#	type::ASCIIString : Name of type.  Unambiguous short forms for the
	#	                    following are acceptable:
	#	                    [bu]tterworth [Default]
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