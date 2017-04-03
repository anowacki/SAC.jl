# Routines to calculate the azimuth, backazimuth and great circle distance
# from headers in traces.

"Earth elliposoid semi-axes in WGS84"
const earth_r_major_WGS84 = 6378137.0000
const earth_r_minor_WGS84 = 6356752.3142
"Flattening of the Earth in WGS84"
const f_WGS84 = (earth_r_major_WGS84 - earth_r_minor_WGS84)/earth_r_major_WGS84

"""
    _great_circle(lon0, lat0, lon1, lat1, f=f_WGS84) -> gcarc, az, baz

Return the great-circle distance, `gcarc`, forward azimuth `az` and backazimuth `baz`
between two points, all specified in degrees.
"""
function _great_circle(lon0, lat0, lon1, lat1, f=f_WGS84)
    lon0, lat0, lon1, lat1 = Float64(lon0), Float64(lat0), Float64(lon1), Float64(lat1)
    gcarc, az, baz = GreatCircle.vincentydist(f, 1.0, deg2rad(lat0), deg2rad(lon0),
                                              deg2rad(lat1), deg2rad(lon1))
    rad2deg(gcarc), rad2deg(az), rad2deg(baz)
end

"""
    update_great_circle!(s::SACtr)

If all headers `evlo`, `evla`, `stlo` and `stla` are set, update the values of
`az`, `baz` and `gcarc`.
"""
function update_great_circle!(s::SACtr)
    any([s.evlo, s.evla, s.stlo, s.stla] .== sac_rnull) && return
    s.gcarc, s.az, s.baz = _great_circle(s.evlo, s.evla, s.stlo, s.stla)
end
