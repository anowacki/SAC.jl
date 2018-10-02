# Test rotation to great circle path
using Test
using SAC

@testset "Great circle path rotation" begin
    let
        # The azimuth and distance of the event
        local az = rand(0:359) # degrees, approximate for now
        local d = 0.1 # degrees

        # Pair of SAC traces representing a spike arriving on the radial component
        local npts = 10
        local delta = 1
        local s = SACtr.([delta,delta], [npts,npts])
        s[:stlo], s[:stla] = 0, 0
        s[:evla] = -d*cosd(az)
        s[:evlo] = -d*sind(az)

        # Make azimuth more accurate as using great circle calculation
        az = mod(s[1][:baz] + 180, 360)
        s[:t] = zeros.(s[:npts])
        # Randomly swap the order of the components, which are in random orientations
        local ne = rand(Bool)
        local az1 = rand(0:359)
        s[:cmpaz] = ne ? [az1, az1+90] : [az1+90, az1]
        s[:kcmpnm] = ne ? ["1", "2"] : ["2", "1"]
        s[:cmpinc] = 90
        local imax = npts÷2
        s[1].t[imax], s[2].t[imax] = ne ? (cosd(az-az1), sind(az-az1)) :
                                          (sind(az-az1), cosd(az-az1))


        # Create all possible combinations of rotations
        local (r, t) = rotate_to_gcp(s[1], s[2])
        local (r′, t′) = rotate_to_gcp(s[2], s[1])
        local rt = rotate_to_gcp(s)
        local rt′ = rotate_to_gcp(s[[2,1]])
        local (R, T) = rotate_to_gcp!(deepcopy(s[1]), deepcopy(s[2]))
        local (R′, T′) = rotate_to_gcp!(deepcopy(s[2]), deepcopy(s[1]))
        local RT = rotate_to_gcp!(deepcopy(s))
        local RT′ = rotate_to_gcp!(deepcopy(s[[2,1]]))

        local list_of_radials = [r, r′, rt[1], rt′[1], R, R′, RT[1], RT′[1]]
        local list_of_transverses = [t, t′, rt[2], rt′[2], T, T′, RT[2], RT′[2]]

        # Test rotation has been done correctly
        local atol = sqrt(eps(SAC.SACFloat))
        for (r1,t1) in zip(list_of_radials, list_of_transverses)
            @test r1.t[imax] ≈ 1
            @test all(isapprox.([r1.t[1:imax-1]; r1.t[imax+1:end]], 0, atol=atol))
            @test all(isapprox.(t1.t, 0, atol=atol))
        end

        # Test order of arguments doesn't stop returns being in order r, t
        for r1 in list_of_radials, r2 in list_of_radials
            @test r1 == r2
        end
        for t1 in list_of_transverses, t2 in list_of_transverses
            @test t1 == t2
        end
    end
end
