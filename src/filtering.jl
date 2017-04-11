# Trace filtering operations

const available_filters = [:butterworth]

"""
    bandpass!(s::SACtr, c1, c2; ftype=:butterworth, npoles=2, passes=1) -> s

Perform a bandpass filter on the SAC trace `s`, between frequency corners `c1`
and `c2`, returning the modified trace.

Select type of filter with `ftype`: current options are: `$(available_filters)`.
Set number of poles with `npoles`.

`passes` may be 1 (forward) or 2 (forward and reverse).
"""
function bandpass!(s::SACtr, c1::Number, c2::Number;
        ftype=:butterworth, npoles::Integer=sac_npoles, passes::Integer=sac_passes)
    c1 >= c2 && error("SAC.bandpass: Upper corner must be larger than lower corner")
    ftype in available_filters || throw(ArgumentError("'$ftype' is not a valid filter type"))
    response = DSP.Bandpass(c1, c2; fs=1./s.delta)
    prototype = get_filter_prototype(ftype, npoles)
    # Create apply the filter
    f = DSP.digitalfilter(response, prototype)
    apply_filter!(s, f, passes)
end

function bandpass!(a::Array{SACtr}, c1, c2; ftype="butterworth", npoles=sac_npoles,
        passes=sac_passes)
    for s in a
        bandpass!(s, c1, c2; ftype=ftype, npoles=npoles, passes=passes)
    end
    a
end
const bp! = bandpass!

"""
    highpass!(s::SACtr, c; ftype=:butterworth, npoles=2, passes=1) -> s

Perform a highpass filter on the SAC trace `s`, above frequency corner `c`,
returning the modified trace.

Select type of filter with `ftype`: current options are: `$(available_filters)`.
Set number of poles with `npoles`.

`passes` may be 1 (forward) or 2 (forward and reverse).
"""
function highpass!(s::SACtr, c::Number;
        ftype=:butterworth, npoles::Integer=sac_npoles,
        passes::Integer=sac_passes)
    ftype in available_filters || throw(ArgumentError("'$ftype' is not a valid filter type"))
    response = DSP.Highpass(c; fs=1./s.delta)
    prototype = get_filter_prototype(ftype, npoles)
    f = DSP.digitalfilter(response, prototype)
    apply_filter!(s, f, passes)
end

function highpass!(a::Array{SACtr}, c;
        ftype=:butterworth, npoles=sac_npoles, passes=sac_passes)
    for s in a
        highpass!(s, c; ftype=ftype, npoles=npoles, passes=passes)
    end
    a
end
const hp! = highpass!

"""
    lowpass!(s::SACtr, c; ftype=:butterworth, npoles=2, passes=1) -> s

Perform a lowpass filter on the SAC trace `s`, above frequency corner `c`,
returning the modified trace.

Select type of filter with `ftype`: current options are: `$(available_filters)`.
Set number of poles with `npoles`.

`passes` may be 1 (forward) or 2 (forward and reverse).
"""
function lowpass!(s::SACtr, c::Number;
        ftype=:butterworth, npoles::Integer=sac_npoles, passes::Integer=sac_passes)
    ftype in available_filters || throw(ArgumentError("'$ftype' is not a valid filter type"))
    response = DSP.Lowpass(c; fs=1./s.delta)
    prototype = get_filter_prototype(ftype, npoles)
    f = DSP.digitalfilter(response, prototype)
    apply_filter!(s, f, passes)
end

function lowpass!(a::Array{SACtr}, c;
        ftype=:butterworth, npoles=sac_npoles, passes=sac_passes)
    for s in a
        lowpass!(s, c; ftype=ftype, npoles=npoles, passes=passes)
    end
    a
end
const lp! = lowpass!

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
end

"""
    get_filter_prototype(ftype, npoles::Integer) -> prototype

Return a filter prototype for use with filtering via `DSP.filt`.

Filter type `ftype` must be one of:
    $(available_filters)

The number of poles of the filter, `npoles`, must be in the range 1–10.
"""
function get_filter_prototype(ftype, npoles::Integer)
    npoles < 1 || npoles > 10 &&
        error("SAC.get_filter_prototype: npoles must be in range 1 - 10")
    ftype = string(ftype)
    if lowercase(ftype[1:2]) == "bu"
        prototype = DSP.Butterworth(npoles)
    else
        error("SAC.get_filter_prototype: unrecognised filter type '$ftype'")
    end
    prototype
end
