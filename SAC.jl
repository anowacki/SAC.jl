module SAC
# Contains routines for dealing with SAC-formatted files

import DSP

export SACtr, fft, rmean!, rtrend!, time, bandpass!, highpass!, write

# SAC unset values
const sac_rnull = -12345.
const sac_inull = -12345
const sac_cnull = "-12345"

# Default values for filtering
const sac_npoles = 2
const sac_passes = 1

# For SAC/BRIS, files are always big-endian, so set this to true to swap by default
const sac_force_swap = true

# Composite type for SAC evenly-spaced time series data
type SACtr
	# Header, floating point.  These are Float32, but we use
	# FloatingPoint internally for ease of use.
	delta::Float32
	depmin::Float32
	depmax::Float32
	scale::Float32
	odelta::Float32
	b::Float32
	e::Float32
	o::Float32
	a::Float32
	internal0::Float32
	t0::Float32
	t1::Float32
	t2::Float32
	t3::Float32
	t4::Float32
	t5::Float32
	t6::Float32
	t7::Float32
	t8::Float32
	t9::Float32
	f::Float32
	resp0::Float32
	resp1::Float32
	resp2::Float32
	resp3::Float32
	resp4::Float32
	resp5::Float32
	resp6::Float32
	resp7::Float32
	resp8::Float32
	resp9::Float32
	stla::Float32
	stlo::Float32
	stel::Float32
	stdp::Float32
	evla::Float32
	evlo::Float32
	evel::Float32
	evdp::Float32
	mag::Float32
	user0::Float32
	user1::Float32
	user2::Float32
	user3::Float32
	user4::Float32
	user5::Float32
	user6::Float32
	user7::Float32
	user8::Float32
	user9::Float32
	dist::Float32
	az::Float32
	baz::Float32
	gcarc::Float32
	internal1::Float32
	internal2::Float32
	depmen::Float32
	cmpaz::Float32
	cmpinc::Float32
	xminimum::Float32
	xmaximum::Float32
	yminimum::Float32
	ymaximum::Float32
	unused1::Float32
	unused2::Float32
	unused3::Float32
	unused4::Float32
	unused5::Float32
	unused6::Float32
	unused7::Float32
	# Integer parts.  These are Int32; again, we use Integer
    nzyear::Int32
	nzjday::Int32
	nzhour::Int32
	nzmin::Int32
	nzsec::Int32
	nzmsec::Int32
	nvhdr::Int32
	norid::Int32
	nevid::Int32
	npts::Int32
	internal3::Int32
	nwfid::Int32
	nxsize::Int32
	nysize::Int32
	unused8::Int32
	iftype::Int32
	idep::Int32
	iztype::Int32
	unused9::Int32
	iinst::Int32
	istreg::Int32
	ievreg::Int32
	ievtyp::Int32
	iqual::Int32
	isynth::Int32
	imagtyp::Int32
	imagsrc::Int32
	unused10::Int32
	unused11::Int32
	unused12::Int32
	unused13::Int32
	unused14::Int32
	unused15::Int32
	unused16::Int32
	unused17::Int32
	# Logical part: boolean
	leven::Bool
	lpspol::Bool
	lovrok::Bool
	lcalda::Bool
	unused18::Bool
	# Character parts
	kstnm::ASCIIString
	kevnm::ASCIIString
	khole::ASCIIString
	ko::ASCIIString
	ka::ASCIIString
	kt0::ASCIIString
	kt1::ASCIIString
	kt2::ASCIIString
	kt3::ASCIIString
	kt4::ASCIIString
	kt5::ASCIIString
	kt6::ASCIIString
	kt7::ASCIIString
	kt8::ASCIIString
	kt9::ASCIIString
	kf::ASCIIString
	kuser0::ASCIIString
	kuser1::ASCIIString
	kuser2::ASCIIString
	kcmpnm::ASCIIString
	knetwk::ASCIIString
	kdatrd::ASCIIString
	kinst::ASCIIString
	# The time series
	t::Array{Float32,1}
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

function read(file; byteswap="auto")
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
			info("SAC.read: Auto-byteswapping file '$file'...")
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

function write(s::SACtr, file::ASCIIString, byteswap=sac_force_swap)
	# Write a SAC composed type to file
	# Call with byteswap=true to write non-native-endian files
	f = open(file, "w")
	if byteswap
		w(F::IOStream, x) = Base.write(F, Base.bswap(x))
	else
		w(F::IOStream, x) = Base.write(F, x)
	end
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
	# No byte-swapping needed for characters
	Base.write(f, s.kstnm)
	Base.write(f, s.kevnm)
	Base.write(f, s.khole)
	Base.write(f, s.ko)
	Base.write(f, s.ka)
	Base.write(f, s.kt0)
	Base.write(f, s.kt1)
	Base.write(f, s.kt2)
	Base.write(f, s.kt3)
	Base.write(f, s.kt4)
	Base.write(f, s.kt5)
	Base.write(f, s.kt6)
	Base.write(f, s.kt7)
	Base.write(f, s.kt8)
	Base.write(f, s.kt9)
	Base.write(f, s.kf)
	Base.write(f, s.kuser0)
	Base.write(f, s.kuser1)
	Base.write(f, s.kuser2)
	Base.write(f, s.kcmpnm)
	Base.write(f, s.knetwk)
	Base.write(f, s.kdatrd)
	Base.write(f, s.kinst)
	# Trace
	for i = 1:s.npts
		w(f, s.t[i])
	end
	close(f)
end	

function sample()
	file = dirname(@__FILE__()) * "/data/seis.sac"
	return read(file)
end

function fft(s::SACtr)
	# Return the fourier-transformed trace and the frequencies to go along with it
	N = int(s.npts/2) + 1
	fmax = 1./(s.npts*s.delta)
	f = [1:N]*fmax
	S = Base.fft(s.t)[1:N]
	return f, S
end

function rmean!(s::SACtr)
	# Remove the mean
	s.t = s.t - mean(s.t)
	update_headers!(s)
end

function rtrend!(s::SACtr)
	# Remove the trend
	t = time(s)
	x0, x1 = linreg(t, s.t)
	s.t = s.t - (x0 + x1*t)
	update_headers!(s)
end
	
function time(s::SACtr)
	# Return an array containing the times for each sample
	return [s.b:s.delta:s.e]
end

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
end

function highpass!(s::SACtr, c::Number;
		ftype::ASCIIString="butterworth", npoles::Integer=sac_npoles,
		passes::Integer=sac_passes)
	# Perform a highpass on the trace.
	response = DSP.Highpass(c; fs=1./s.delta)
	prototype = get_filter_prototype(ftype, npoles)
	f = DSP.digitalfilter(response, prototype)
	apply_filter!(s, f, passes)
end

function lowpass!(s::SACtr, c::Number;
		ftype::ASCIIString="butterworth", npoles::Integer=sac_npoles,
		passes::Integer=sac_passes)
	# Perform a lowpass on the trace.
	response = DSP.Lowpass(c; fs=1./s.delta)
	prototype = get_filter_prototype(ftype, npoles)
	f = DSP.digitalfilter(response, prototype)
	apply_filter!(s, f, passes)
end

function apply_filter!(s::SACtr, f::DSP.ZPKFilter{Complex{Float64},Complex{Float64},Float64},
		passes::Integer)
		passes < 1 || passes > 2 && error("SAC.apply_filter!: Number of passes must be 1 or 2")
	if passes == 1
		DSP.filt!(s.t, f, s.t)
	elseif passes == 2
		s.t = DSP.filtfilt(f, s.t)
	else
		error("SAC.apply_filter!: passes must be 1 or 2")
	end
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