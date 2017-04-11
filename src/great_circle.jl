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
    gcarc, az, baz = vincentydist(f, 1.0, deg2rad(lat0), deg2rad(lon0),
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

# The following functions have been taken from
#   https://github.com/acrosby/GreatCircle.jl
# commit
#   c4860071103c1895a86561c0a0a772bf89705fba
#
# See the file LICENSE.md for details of the licence.
#
# -----------------------------------------------------------------------
# | Algorithms from Geocentric Datum of Australia Technical Manual |
# | |
# | http://www.anzlic.org.au/icsm/gdatum/chapter4.html |
# | |
# | This page last updated 11 May 1999 |
# | |
# | Computations on the Ellipsoid |
# | |
# | There are a number of formulae that are available |
# | to calculate accurate geodetic positions, |
# | azimuths and distances on the ellipsoid. |
# | |
# | Vincenty's formulae (Vincenty, 1975) may be used |
# | for lines ranging from a few cm to nearly 20,000 km, |
# | with millimetre accuracy. |
# | The formulae have been extensively tested |
# | for the Australian region, by comparison with results |
# | from other formulae (Rainsford, 1955 & Sodano, 1965). |
# | |
# | * Inverse problem: azimuth and distance from known |
# | latitudes and longitudes |
# | * Direct problem: Latitude and longitude from known |
# | position, azimuth and distance. |
# | * Sample data |
# | * Excel spreadsheet |
# | |
# | Vincenty's Inverse formulae |
# | Given: latitude and longitude of two points |
# | (phi1, lembda1 and phi2, lembda2), |
# | Calculate: the ellipsoidal distance (s) and |
# | forward and reverse azimuths between the points (alpha12, alpha21). |
# | |
# -----------------------------------------------------------------------
"""
    vincentypt(flattening, semimajor_radius, lat, lon, azimuth, distance) -> lat, lon, backazimuth

Return lat and long of projected point and backazimuth,
given a reference point and a distance and azimuth to project.

Coordinates and azimuth are in radians.
"""
function vincentypt{T<:AbstractFloat}(f::T, a::T, phi1::T, lembda1::T, alpha12::T, s::T)
    const two_pi = 2π
    if alpha12 < 0.0
        alpha12 = alpha12 + two_pi
    end
    if alpha12 > two_pi
        alpha12 = alpha12 - two_pi
    end

    b = a .* (1.0 - f)

    TanU1 = (1 - f) .* tan(phi1)
    U1 = atan( TanU1 )
    sigma1 = atan2( TanU1, cos(alpha12) )
    Sinalpha = cos(U1) .* sin(alpha12)
    cosalpha_sq = 1.0 - Sinalpha .* Sinalpha

    u2 = cosalpha_sq .* (a .* a - b .* b ) ./ (b .* b)
    A = 1.0 + (u2 ./ 16384) .* (4096 + u2 .* (-768 + u2 .* (320 - 175 .* u2) ) )
    B = (u2 ./ 1024) .* (256 + u2 .* (-128 + u2 .* (74 - 47 .* u2) ) )

    # Starting with the approximation
    sigma = (s ./ (b .* A))

    # Not moving anywhere. We can return the location that was passed in.
    if sigma == 0
        return phi1, lembda1, alpha12
    end

    last_sigma = 2.0 .* sigma + 2.0 # something impossible

    # Iterate the following three equations
    # until there is no significant change in sigma
    # two_sigma_m , delta_sigma
    while ( abs( (last_sigma - sigma) ./ sigma) > 1.0e-9 )
        global two_sigma_m = 2 .* sigma1 + sigma
        delta_sigma = B .* sin(sigma) .* ( cos(two_sigma_m) + (B./4) .* (cos(sigma) .* (-1 + 2 .* cos(two_sigma_m).^2 - (B./6) .* cos(two_sigma_m) .* (-3 + 4 .* sin(sigma).^2) .* (-3 + 4 .* cos(two_sigma_m).^2 ))))
        last_sigma = sigma
        sigma = (s ./ (b .* A)) + delta_sigma
    end

    phi2 = atan2( (sin(U1) .* cos(sigma) + cos(U1) .* sin(sigma) .* cos(alpha12) ), ((1-f) .* sqrt( Sinalpha.^2 + (sin(U1) .* sin(sigma) - cos(U1) .* cos(sigma) .* cos(alpha12)).^2)))

    lembda = atan2( (sin(sigma) .* sin(alpha12 )), (cos(U1) .* cos(sigma) - sin(U1) .* sin(sigma) .* cos(alpha12)))

    C = (f./16) .* cosalpha_sq .* (4 + f .* (4 - 3 .* cosalpha_sq ))

    omega = lembda - (1-C) .* f .* Sinalpha .* (sigma + C .* sin(sigma) .* (cos(two_sigma_m) + C .* cos(sigma) .* (-1 + 2 .* cos(two_sigma_m).^2 )))

    lembda2 = lembda1 + omega

    alpha21 = atan2( Sinalpha, (-sin(U1) .* sin(sigma) + cos(U1) .* cos(sigma) .* cos(alpha12)))
    alpha21 = alpha21 + pi

    if alpha21 < 0.0
        alpha21 = alpha21 + two_pi
    end
    if alpha21 > two_pi
        alpha21 = alpha21 - two_pi
    end
    return phi2, lembda2, alpha21

end

function vincentypt{T<:AbstractFloat}(f::T, a::T, phi1::Array{T,1}, lembda1::Array{T,1}, alpha12::Array{T,1}, s::Array{T,1})
    lenphi = length(phi1)
    @assert lenphi == length(lembda1)
    @assert lenphi == length(alpha12)
    @assert lenphi == length(s)
    phi2 = zeros(phi1)
    lembda2 = zeros(phi1)
    alpha21 = zeros(phi1)
    for i=1:lenphi
        @inbounds phi2[i], lembda2[i], alpha21[i] = vincentypt(f, a, phi1[i], lembda1[i], alpha12[i], s[i])
    end
    return phi2, lembda2, alpha21
end

vincentypt{T<:AbstractFloat}(f::T, a::T, phi1::Array{T,1}, lembda1::Array{T,1}, alpha12::Array{T,1}, s::T) = vincentypt(f, a, phi1, lembda1, alpha12, ones(length(phi1)) * s)
vincentypt{T<:AbstractFloat}(f::T, a::T, phi1::Array{T,1}, lembda1::Array{T,1}, alpha12::T, s::Array{T,1}) = vincentypt(f, a, phi1, lembda1, ones(length(phi1)) * alpha12, s)
vincentypt{T<:AbstractFloat}(f::T, a::T, phi1::Array{T,1}, lembda1::Array{T,1}, alpha12::T, s::T) = vincentypt(f, a, phi1, lembda1, ones(length(phi1)) * alpha12, ones(length(phi1)) * s)


"""
    vincentydist(flattening, semimajor_radius, lat1, lon1, lat2, lon2) -> dist, azimuth, backazimuth

Returns the distance between two geographic points on the ellipsoid,
and the forward and reverse azimuths between these points.

Latitude, longitude and (back)azimuth are in radians; distances are in the units
provided for the semimajor axis length.
"""
function vincentydist{T<:AbstractFloat}(f::T, a::T, phi1::T, lembda1::T, phi2::T, lembda2::T)
    const two_pi = 2π
    if (abs( phi2 - phi1 ) < 1e-8) & ( abs( lembda2 - lembda1) < 1e-8 )
        return 0.0, 0.0, 0.0
    end

    b = a .* (1.0 - f)

    TanU1 = (1 - f) .* tan( phi1 )
    TanU2 = (1 - f) .* tan( phi2 )

    U1 = atan(TanU1)
    U2 = atan(TanU2)

    lembda = lembda2 - lembda1
    last_lembda = -4000000.0 # an impossibe value
    omega = lembda

    # Iterate the following equations until there is no significant change in lembda
    alpha, sigma, Sin_sigma, Cos2sigma_m, Cos_sigma, sqr_sin_sigma = -999999., -999999., -999999., -999999., -999999., -999999.
    while ( (last_lembda < -3000000.0) | (lembda != 0) ) & ( abs( (last_lembda - lembda)./lembda) > 1.0e-9 )
        sqr_sin_sigma = (cos(U2) .* sin(lembda)).^2 + ( (cos(U1) .* sin(U2) - sin(U1) .* cos(U2) .* cos(lembda) )).^2
        Sin_sigma = sqrt( sqr_sin_sigma )
        Cos_sigma = sin(U1) .* sin(U2) + cos(U1) .* cos(U2) .* cos(lembda)
        sigma = atan2( Sin_sigma, Cos_sigma )

        Sin_alpha = cos(U1) .* cos(U2) .* sin(lembda) ./ sin(sigma)

        if (Sin_alpha >= 1)# & (Sin_alpha .== 1.0)
            Sin_alpha = 1.0
        elseif (Sin_alpha <= -1)# & (Sin_alpha .== -1.0)
            Sin_alpha = -1.0
        end

        alpha = asin( Sin_alpha )
        Cos2sigma_m = cos(sigma) - (2 .* sin(U1) .* sin(U2) ./ cos(alpha).^2 )
        C = (f./16) .* cos(alpha).^2 .* (4 + f .* (4 - 3 .* cos(alpha).^2))
        last_lembda = lembda
        lembda = omega + (1-C) .* f .* sin(alpha) .* (sigma + C .* sin(sigma) .* (Cos2sigma_m + C .* cos(sigma) .* (-1 + 2 .* Cos2sigma_m.^2 )))
    end

    u2 = cos(alpha).^2 .* (a.*a-b.*b) ./ (b.*b)

    A = 1 + (u2./16384) .* (4096 + u2 .* (-768 + u2 .* (320 - 175 .* u2)))

    B = (u2./1024) .* (256 + u2 .* (-128+ u2 .* (74 - 47 .* u2)))

    delta_sigma = B .* Sin_sigma .* (Cos2sigma_m + (B./4) .* (Cos_sigma .* (-1 + 2 .* Cos2sigma_m.^2 ) - (B./6) .* Cos2sigma_m .* (-3 + 4 .* sqr_sin_sigma) .* (-3 + 4 .* Cos2sigma_m.^2 )))

    s = b .* A .* (sigma - delta_sigma)

    alpha12 = atan2( (cos(U2) .* sin(lembda)), (cos(U1) .* sin(U2) - sin(U1) .* cos(U2) .* cos(lembda)))
    alpha21 = atan2( (cos(U1) .* sin(lembda)), (-sin(U1) .* cos(U2) + cos(U1) .* sin(U2) .* cos(lembda)))

    if alpha12 < 0.0
        alpha12 = alpha12 + two_pi
    end
    if alpha12 > two_pi
        alpha12 = alpha12 - two_pi
    end
    alpha21 = alpha21 + two_pi ./ 2.0

    if alpha21 < 0.0
        alpha21 = alpha21 + two_pi
    end
    if alpha21 > two_pi
        alpha21 = alpha21 - two_pi
    end
    return s, alpha12, alpha21
end

function vincentydist{T<:AbstractFloat}(f::T, a::T, phi1::Array{T,1}, lembda1::Array{T,1}, phi2::Array{T,1}, lembda2::Array{T,1})
    lenphi = length(phi1)
    @assert lenphi == length(lembda1)
    @assert lenphi == length(phi2)
    @assert lenphi == length(lembda2)
    s = zeros(phi1)
    alpha12 = zeros(phi1)
    alpha21 = zeros(phi1)
    for i=1:lenphi
        @inbounds s[i], alpha12[i], alpha21[i] = vincentydist(f, a, phi1[i], lembda1[i], phi2[i], lembda2[i])
    end
    return s, alpha12, alpha21
end

vincentydist{T<:AbstractFloat}(f::T, a::T, phi1::T, lembda1::Array{T,1}, phi2::T, lembda2::Array{T,1}) = vincentydist(f, a, phi1 * ones(lembda1), lembda1, phi2 * ones(lembda1), lembda2)
vincentydist{T<:AbstractFloat}(f::T, a::T, phi1::Array{T,1}, lembda1::T, phi2::Array{T,1}, lembda2::T) = vincentydist(f, a, phi1, lembda1 * ones(phi1), phi2, lembda2 * ones(phi1))
