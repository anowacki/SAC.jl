# Test statistical measures of traces
using SAC, Base.Test

@testset "Statistics" begin
    @testset "Window average" begin
        # Trace of constant absolute value of 1 with random polarity
        let delta = 0.1, npts = 1001
            s = SACtr(rand([-1, 1], npts), delta)
            @test window_average(s, s.b, s.delta*(s.npts-1)) ≈ 1.0
            @test_throws ErrorException window_average(s, s.b-s.delta, s.delta*(s.npts))
            @test_throws ErrorException window_average(s, s.b, -1)
        end
        # Trace of amplitude 1 at 1 sample after zero
        let delta = 1.0, npts = 101, b = -50.0
            s = SACtr(zeros(npts), delta, b)
            s.t[52] = 1
            @test window_average(s, s.b, 40delta) ≈
                window_average(s, 2delta, 40delta) atol=eps(SAC.SACFloat)
            @test window_average(s, 0, 2delta) ≈ 1/3
            @test window_average(s, delta, 2delta, true) ≈ 1/3
        end
    end

    @testset "STA/LTA" begin
        # Step function from 1 to 2 one sample after 0 time
        let delta = 1.0, npts = 101, b = -50.0
            s = SACtr(ones(npts), delta, b)
            s.t[52:end] = 2
            @test stalta(s, delta, 10, 10) ≈ 2.0
            @test stalta(s, -delta, 2delta, 0.5, 0.5) ≈ [1, 1, 2, 2]
            @test typeof(stalta(s, 0.1, 1)) == SACtr
            t1, t2 = 2, 4
            s′ = stalta(s, t1, t2)
            n = npts - t1 - t2 - 3
            @test s′.npts == n
            @test time(s′)[indmax(s′.t)] ≈ 1.0
        end
    end
end
