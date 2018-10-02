using Test
using SAC

@testset "Utility" begin
    # Component flipping
    let s = SACtr(rand(3), 1)
        # Horizontal component pointing at 15°
        s[:cmpaz], s[:cmpinc] = 15, 90
        @test flip_component(s)[:cmpaz] ≈ 195
        @test flip_component(s)[:cmpinc] ≈ 90
        # Arbitrary component
        s[:cmpaz], s[:cmpinc] = 271, 30
        s′ = flip_component!(s)
        @test s′[:cmpaz] ≈ 91
        @test s′[:cmpinc] ≈ 150
        @test flip(s) == flip_component(s)
    end
end