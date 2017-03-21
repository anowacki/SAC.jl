# Test the integrate function
using SAC, Base.Test

sample_data = SAC.sample()

# Don't accept wrong arguments
@test_throws ArgumentError SAC.integrate(sample_data, :weird_method)

# Compare against output from SAC
for (method, name) in zip((:trapezium, :rectangle), ("trapezoidal", "rectangular"))
    file = joinpath(dirname(@__FILE__), "data", "test_integrate_data_method_$(name).sac")
    isfile(file) || error("Test data file '$(file)' for differentiate does not exist")
    reference_data = SAC.read(file)
    test_data = SAC.integrate(SAC.rmean(sample_data), method)
    @test test_data.t ≈ reference_data.t
    @test test_data.b ≈ reference_data.b
end
