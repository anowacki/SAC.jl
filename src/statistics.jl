"""
    stalta(s::SACtr, t, ws, wl) -> ::SACFloat

Return the ratio of the short-term average (window length `ws` s) to long-term
average (length `wl` s) at time `t` s.

    stalta(s::SACtr, t1, t2, ws, wl) -> ::Vector{SACFloat}

Return STA/LTA for each sample between times `t1` and `t2` s.

    stalta(s::SACtr, ws, wl) -> ::SACtr

Return STA/LTA for whole trace as a `SACtr` object.

#### Note on window positions

In this implemetation, the ratio refers to that between a long time window *before*
time `t` to a short window *after* `t`, and here the short window includes the sample
*at* `t` (or the nearest thereto).
"""
stalta(s::SACtr, t, ws, wl) = window_average(s, t, ws)/window_average(s, t-wl-s.delta, wl)

function stalta(s::SACtr, t1, t2, ws, wl)
   n = round(Int, (t2 - t1)/s.delta) + 1
   ratio = Array{eltype(s.t)}(undef, n)
   for i in 1:n
       t = t1 + (i-1)*s.delta
       ratio[i] = stalta(s, t, ws, wl)
   end
   ratio
end

function stalta(s::SACtr, ws, wl)
    t1 = s.b + wl + 2s.delta
    t2 = s.e - ws - s.delta
    t = stalta(s, t1, t2, ws, wl)
    s2 = deepcopy(s)
    s2.t = t
    s2.npts = length(t)
    s2.b = t1
    update_headers!(s2)
end

"""
    window_average(s::SACtr, ref_time, window_length, centred=false)

Return the mean value of the absolute amplitude within a window
`window_length` s long.  `ref_time` in s may be the start of the window (default)
or the centre if `centred` is `true`.

Note that window times down to 0 are allowed, in which case the amplitude of
the nearest point is returned.
"""
function window_average(s::SACtr, time, window_length, centred=false)
    window_length >= 0 || error("`window_length` ($window_length) must be positive")
    n = round(Int, window_length/s.delta) + 1
    it1 = round(Int, (time - s.b - (centred ? window_length/2 : 0.0))/s.delta) + 1
    it2 = it1 + n - 1
    it1 > 0 && it2 <= s.npts ||
        error("Window of times $time - $(time+window_length) s lies outside " *
              "trace (trace indices $it1 - $it2; start time $(s.b) s; end time $(s.e)) s")
    ave = zero(SACFloat)
    @inbounds for i in it1:it2
        ave += abs(s.t[i])
    end
    ave/n
end
