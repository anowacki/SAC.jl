# Test stacking of SAC traces
using SAC, Base.Test

normal_dist(x, μ, σ) = exp(-(x - μ)/2σ^2)/sqrt(2σ^2*π)

@testset "Stacking" begin
    # Create five traces which are just a normal distribution and stack the traces
    # with shifts of the means, meaning the stack should be just a normal distribution
    # itself.
    @test begin
        N = 5
        μ = 1:5
        σ = 0.1

        delta = 0.01
        b = 0.
        e = 10
        npts = round(Int, (e - b)/delta) + 1

        S = Array{SACtr}(N)
        for i in 1:N
            S[i] = SACtr(delta, npts, b)
            S[i].t[:] = normal_dist.(time(S[i]), μ[i], σ)
            S[i].a = μ[i]
            S[i].ka = "mu"
            SAC.update_headers!(S[i])
        end

        # Create the stack
        t1 = -0.5
        t2 = 0.5
        st = stack(S, (-0.5,0.5), collect(μ), ones(N))
        st.t ≈ normal_dist.(t1:delta:t2, 0, σ)
    end
end
