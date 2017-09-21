# Test reading and writing of SAC files
using SAC, Base.Test

@testset "IO" begin
    # Reading
    @test SAC.sample().t[1:100] ≈ Float32[-0.09728,-0.09728,-0.09856,-0.09856,-0.09728,-0.096,-0.09472,-0.09344,-0.09344,-0.09344,-0.09344,-0.09344,-0.09472,-0.09472,-0.09344,-0.09344,-0.09216,-0.09216,-0.09216,-0.09216,-0.09088,-0.09088,-0.09216,-0.09344,-0.09472,-0.09472,-0.09472,-0.09472,-0.09472,-0.09472,-0.09344,-0.09344,-0.09216,-0.09088,-0.09088,-0.09216,-0.09216,-0.09216,-0.09344,-0.09472,-0.096,-0.09856,-0.09856,-0.09856,-0.09728,-0.09728,-0.09856,-0.09984,-0.09984,-0.09984,-0.09984,-0.09984,-0.10112,-0.10112,-0.10112,-0.10112,-0.1024,-0.1024,-0.10368,-0.1024,-0.10496,-0.10496,-0.10624,-0.10368,-0.10368,-0.1024,-0.10368,-0.10368,-0.10368,-0.10368,-0.10496,-0.10624,-0.10624,-0.10496,-0.10368,-0.10368,-0.10496,-0.10624,-0.10496,-0.10496,-0.10496,-0.10624,-0.1088,-0.1088,-0.1088,-0.1088,-0.1088,-0.10624,-0.10496,-0.10496,-0.10496,-0.10624,-0.10624,-0.10624,-0.10368,-0.10368,-0.10368,-0.10368,-0.1024,-0.10112]

    # Writing
    tempfile = tempname()
    @test begin
        s1 = SAC.sample()
        write([s1], [tempfile])
        s2 = SAC.read(tempfile)
        s1 == s2
    end

    # Cutting
    @test begin
        s1 = SAC.read(tempfile)
        cut!(s1, 12, 13)
        s2 = SAC.read_cut(tempfile, 12, 13)
        s1 == s2
    end

    # Printing
    @test begin
        open(tempfile, "w") do f
            show(f, SAC.SACtr(1, 1))
        end
        readstring(tempfile) == "SAC.SACtr(delta=1.0, b=0.0, npts=1)"
    end
    teststring = VERSION < v"0.6-" ?
        "SAC.SACtr[SAC.SACtr(delta=1.0, b=0.0, npts=1),SAC.SACtr(delta=1.0, b=0.0, npts=1)]" :
        "SAC.SACtr[SAC.SACtr(delta=1.0, b=0.0, npts=1), SAC.SACtr(delta=1.0, b=0.0, npts=1)]"
    @test begin
        open(tempfile, "w") do f
            show(f, [SAC.SACtr(1, 1), SAC.SACtr(1, 1)])
        end
        readstring(tempfile) == teststring
    end
end
