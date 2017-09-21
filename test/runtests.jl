using Base.Test

@testset "SAC tests" begin
    include("construction.jl")
    include("io.jl")
    include("diff.jl")
    include("integrate.jl")
    include("stack.jl")
    include("great_circle.jl")
    include("rotate_to_gcp.jl")
end
