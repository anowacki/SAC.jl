# Test the creation of SACtr's
using Compat.Test
using SAC

@testset "Construction" begin
    let npts = 1000, delta = 1, b = -5
        t = rand(SAC.SACFloat, npts)
        s = SACtr(t, delta)
        @test s.t == t
        @test s.delta == delta
        @test s.npts == npts

        s = SACtr(t, delta, b)
        @test s.b == b
        @test s.e ≈ b + (npts - 1)*delta
    end
end
