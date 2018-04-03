# Test the creation of SACtr's
using Compat.Test
using SAC

@testset "Construction" begin
    # Construction from array
    let npts = 1000, delta = 1, b = -5
        t = rand(SAC.SACFloat, npts)
        s = SACtr(t, delta)
        @test s.t == t
        @test s.delta == delta
        @test s.npts == npts

        s = SACtr(t, delta, b)
        @test s.b == b
        @test s.e ≈ b + (npts - 1)*delta
        @test s.depmax == maximum(t)
        @test s.depmin == minimum(t)
        @test s.depmen == mean(t)
    end

    # Construction without array
    let npts = rand(1:1000), b = 2, delta = 0.1
        @test SACtr(delta, npts).npts == npts
        @test SACtr(delta, npts).b ≈ 0.0
        @test SACtr(delta, npts, b).b ≈ b
    end

    ## Field modification and access
    # Single traces
    let a = rand(SAC.SACFloat, 100), delta = rand(), s = SACtr(a, delta)
        s[:kevnm] = "ABCD"
        @test s.kevnm == "ABCD"
        @test s[:kevnm] == "ABCD"
        s[:user0] = 1
        @test s.user0 ≈ 1
        @test s[:user0] ≈ 1
        s[:t] = 1:100
        @test s.t ≈ 1:100
        @test s.t ≈ s[:t]
        s[:t] = 2:101
        @test s.t ≈ 2:101
        @test s[:depmin] ≈ 2
        @test s[:depmax] ≈ 101
        @test s[:depmen] ≈ mean(2:101)
        k = 1
        s[:t] = k
        @test s[:depmin] ≈ k
        @test s[:depmax] ≈ k
        @test all(s[:t] .== SAC.SACFloat(k))
        k = 1:100
        s[:t] = k
        @test all(s[:t] .== SAC.SACFloat.(k))
    end

    # Arrays of traces
    let N = 5, b = -10, delta = 2, npts = 100
        a = [SACtr(delta, npts, b) for _ in 1:N]
        @test all(a[:b] .≈ b)
        @test all(a[:delta] .≈ delta)
        @test all(a[:npts] .== npts)
        a[:a] = 1:N
        @test all(a[:a] .== SAC.SACFloat.(1:N))
        a[:t] = [1:npts for _ in 1:N]
        @test all(a[:t] .== [SAC.SACFloat.(1:npts) for _ in 1:N])
        k = -1
        a[:t] = k
        @test all(a[:t] .== [[SAC.SACFloat(k) for j in 1:npts] for i in 1:N])
        a[:kcmpnm] = "ABC"
        @test all(a[:kcmpnm] .== "ABC")
        a[:kcmpnm] = ["a", "b", "c", "d", "e"]
        @test all(a[:kcmpnm] .== ["a", "b", "c", "d", "e"])
    end
end
