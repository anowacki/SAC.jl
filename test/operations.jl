# Test linear operations on SAC traces
using Test
using SAC
if VERSION >= v"0.7-"
    import FFTW: rfft, irfft
end

function inplace_compare(s, f, f!, val)
    s = deepcopy(s)
    s′ = deepcopy(s)
    f!(s, val)
    f(s′, val) == s
end

function inplace_raw_compare(a::AbstractVector, b, delta, f!, val)
    s = SACtr(a, delta)
    f!(s, val)
    all(s.t .== b)
end

function inplace_raw_compare(a::AbstractVector, b, delta, f)
    s = SACtr(a, delta)
    modify!(s, f)
    all(s.t .== b)
end


@testset "Operations" begin
    ## In-place versions
    # Compare copying and inplace operations
    let s = SAC.sample(), val = rand()
        @test inplace_compare(s, add, add!, val)
        @test inplace_compare(s, subtract, subtract!, val)
        @test inplace_compare(s, divide, divide!, val)
        @test inplace_compare(s, multiply, multiply!, val)
    end
    # Compare inplace operations with raw array operations
    let a = rand(SAC.SACFloat, 100), val = rand(SAC.SACFloat)-SAC.SACFloat(2), delta = 0.01
        @test inplace_raw_compare(a, a.+val, delta, add!, val)
        @test inplace_raw_compare(a, a.-val, delta, add!, -val)
        @test inplace_raw_compare(a, a.*val, delta, multiply!, val)
        @test inplace_raw_compare(a, a./val, delta, divide!, val)
    end
    # General modification
    let a = rand(SAC.SACFloat, 100), delta = 2
        @test inplace_raw_compare(a, a.^2, delta, x->x^2)
        @test inplace_raw_compare(a, rem.(a, SAC.SACFloat(0.1)), delta, x->rem(x, SAC.SACFloat(0.1)))
        @test inplace_raw_compare(a, 2 .* (a .- 1)/20, delta, x->2*(x-1)/20)
    end

    ## Copying versions
    # Trace and constant
    let s = SAC.sample(), val = rand(SAC.SACFloat)
        # Addition/subtraction: raw values
        @test typeof(s + val) == SACtr
        @test all(s.t .+ val .== (s + val).t)
        @test all(s.t .+ val .== (val + s).t)
        @test all(val .+ s.t .== (val + s).t)
        @test all(val .+ s.t .== (s + val).t)
        @test s + 0 == s
        @test typeof(s - val) == SACtr
        @test all(s.t .- val .== (s - val).t)
        @test all(s.t .- val .== (-val + s).t)
        @test all(val .- s.t .== (val - s).t)
        @test all(val .- s.t .== (-s + val).t)
        @test s - 0 == s
        # Order
        @test s + val == val + s
        @test s - val == -val + s
        # Multiplication/divison: raw values
        val += SAC.SACFloat(1) # Ensure we're not dividing by zero
        @test all(s.t .* val .== (s*val).t)
        @test all(s.t .* val .== (val*s).t)
        @test s*1 == s
        @test all(s.t ./ val .== (s/val).t)
        @test all(val ./ s.t .== (val/s).t)
        @test s/1 == s
        # Order
        @test s*val == val*s
        @test s/val == s*1/val
        @test val/s == val*1/s
    end
    # General modification
    let s = SACtr(rand(100), 2.0)
        @test all(modify(s, x->x^2).t .== s.t.^2)
        @test isapprox(modify(s, x->rem(x, 0.1)).t, rem.(s.t, 0.1))
    end

    # Array of traces and constant
    let a = sample(:array), vals = rand(SAC.SACFloat, length(a))
        local N = length(a)
        # Array with constant
        val = vals[1]
        @test typeof(a .+ val) == Array{SACtr,1}
        n = rand(1:N)
        @test (a .+ val)[n] == a[n] + val
        @test (a .- val)[n] == a[n] - val
        @test (a .* val)[n] == a[n] * val
        @test (a ./ val)[n] == a[n] / val
        # Array with array
        @test typeof(a .+ vals) == Array{SACtr,1}
        @test a .+ vals == vals .+ a
        @test typeof(a .- vals) == Array{SACtr,1}
        @test a .- vals == -vals .+ a
        @test a .* vals == vals .* a
        @test typeof(a ./ vals) == Array{SACtr,1}
        @test a ./ vals == a .* 1 ./ vals
        @test all(a .+ zeros(N) .== a)
        @test all(a .- zeros(N) .== a)
        @test all(a .* ones(N) .== a)
        @test all(a ./ ones(N) .== a)
    end

    # Multiple trace versions
    let s1 = SAC.sample(), s2 = SAC.sample() + 1
        @test_throws ArgumentError SAC.check_traces_combine(s1, cut(s1, :a, 1, :e, -1))
        s1′ = deepcopy(s1)
        s1′.delta = 2s1.delta
        @test_throws ArgumentError SAC.check_traces_combine(s1, s1′)
        @test SAC.check_traces_combine(s1, s1)
        @test SAC.check_traces_combine(s1, s2)
        @test typeof(s1 + s2) == SACtr
        @test all((s1 + s2).t .== s1.t .+ s2.t)
        @test all((s1 - s2).t .== s1.t .- s2.t)
        @test all((s1 * s2).t .== s1.t .* s2.t)
        @test all((s1 / s2).t .== s1.t ./ s2.t)
        # Order
        @test s1 + s2 == s2 + s1
        @test s1 - s2 == -s2 + s1
        @test s1 * s2 == s2 * s1
        @test s1/s2 == s1*1/s2
    end
    
    ## Trend and mean
    # Single trace
    let s = SACtr([i%2 for i in 1:1000], 0.01)
        @test s[:depmen] ≈ 0.5
        @test rmean(s)[:depmen] ≈ 0.0
        rmean!(s)
        @test s[:depmen] ≈ 0.0
    end
    let s = SACtr(1:100, 0.1)
        @test all(isapprox.(rtrend(s)[:t], 0.0, atol=1e-4))
        rtrend!(s)
        @test all(isapprox.(s[:t], 0.0, atol=1e-4))
    end
    # Arrays
    let a = [SAC.sample() for _ in 1:5], s = SAC.sample()
        @test all(rmean(a) .== rmean(s))
        rmean!(a)
        rmean!(s)
        @test all(a .== s)
    end
    let a = [SAC.sample() for _ in 1:5], s = SAC.sample()
        @test all(rtrend(a) .== rtrend(s))
        rtrend!(a)
        rtrend!(s)
        @test all(s .== s)
    end

    ## Normalisation
    let s = SACtr(rand(SAC.SACFloat, 100), 0.1), s′ = deepcopy(s)
        @test all(normalise(s).t .== s′.t./maximum(abs, s.t))
        @test normalise!(s) == normalise(s′)
    end

    ## FFT
    let s = SACtr(rand(SAC.SACFloat, 10000), 0.01)
        f, S = fft(s)
        @test s.t ≈ irfft(S, s[:npts])
        @test length(f) == length(S) == s[:npts]÷2 + 1
        @test minimum(f) == zero(typeof(f[1]))
        @test maximum(f) ≈ 1/s[:delta]/2
        # Array version
        ss = [s for _ in 1:3]
        f, S = fft(ss)
        @test length(f) == length(S) == length(ss)
        @test all((f[1],) .== f)
        @test all((S[1],) .== S)
    end
end
