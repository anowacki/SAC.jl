# Contains constants and constructors for SAC file precision, and defaults
# for routines and endianness

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
# so this determines whether a file is 'native-endian' or not.
const machine_is_little_endian = bits(1)[end] == '1'

# Convert a number into a SACChar
sacstring(x, maxlen=saccharlen) = sacchar(string(x)[1:minimum((length(string(x)),maxlen))]*" "^(maximum((0,maxlen-length(string(x))))))

# SAC unset values
const sac_rnull = -12345.
const sac_inull = -12345
const sac_cnull = "-12345"

# Default values for filtering
const sac_npoles = 2
const sac_passes = 1

# For SAC/BRIS (MacSAC), files are always big-endian, so set this appropriately
const sac_force_swap = machine_is_little_endian

# Flag for verbosity
const sac_verbose = Ref(true)

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
