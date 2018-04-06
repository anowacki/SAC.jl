# Test the differentiate function
using Compat.Test
using SAC

@testset "Differentiation" begin
    let sample_data = SAC.sample()

        # Don't accept wrong arguments
        for i in (0, 1, 4, 6)
            @test_throws ArgumentError SAC.differentiate(sample_data, i)
        end

        # Compare against output from SAC
        for (i, str) in zip((2, 3, 5), ("two", "three", "five"))
            local file = joinpath(dirname(@__FILE__), "data", "test_diff_data_points_$(str).sac")
            isfile(file) || error("Test data file '$(file)' for differentiate does not exist")
            local reference_data = SAC.read(file)
            @test SAC.differentiate(sample_data, i).t ≈ reference_data.t
        end
    end
end
