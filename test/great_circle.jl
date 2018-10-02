# Test the grat circle calculation functions
using Test
using SAC

@testset "Great circles" begin
    let s = SACtr(1, 2)
        s[:evlo], s[:evla], s[:stlo], s[:stla] = 0, 0, 12, 15
        @test s[:gcarc] ≈ 19.047431084635505
        @test s[:az] ≈ 37.992268575139384
        @test s[:baz] ≈ 219.57789440455088
    end
end
