# Test the grat circle calculation functions
using SAC, Base.Test

@test begin
    s = SACtr(1, 2)
    s[:evlo], s[:evla], s[:stlo], s[:stla] = 0, 0, 12, 15
    s[:gcarc] ≈ 19.047431084635505 &&
        s[:az] ≈ 37.992268575139384 &&
        s[:baz] ≈ 219.57789440455088
end
