# Test the integrate function
using Compat.Test
using SAC

@testset "Integration" begin
    let sample_data = SAC.sample()

        # Don't accept wrong arguments
        @test_throws ArgumentError SAC.integrate(sample_data, :weird_method)

        # Compare against output from SAC
        for (method, name) in zip((:trapezium, :rectangle), ("trapezoidal", "rectangular"))
            local file = joinpath(dirname(@__FILE__), "data", "test_integrate_data_method_$(name).sac")
            isfile(file) || error("Test data file '$(file)' for differentiate does not exist")
            local reference_data = SAC.read(file)
            local test_data = SAC.integrate(SAC.rmean(sample_data), method)
            @test test_data.t ≈ reference_data.t
            @test test_data.b ≈ reference_data.b
        end
    end
end
