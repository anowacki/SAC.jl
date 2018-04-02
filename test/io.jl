# Test reading and writing of SAC files
using Compat.Test
using SAC

@testset "IO" begin
    tempfile = tempname()

    # Reading
    let s = SAC.sample()
        # Reals
        @test s.t[1:100] ≈ Float32[-0.09728,-0.09728,-0.09856,-0.09856,-0.09728,-0.096,-0.09472,-0.09344,-0.09344,-0.09344,-0.09344,-0.09344,-0.09472,-0.09472,-0.09344,-0.09344,-0.09216,-0.09216,-0.09216,-0.09216,-0.09088,-0.09088,-0.09216,-0.09344,-0.09472,-0.09472,-0.09472,-0.09472,-0.09472,-0.09472,-0.09344,-0.09344,-0.09216,-0.09088,-0.09088,-0.09216,-0.09216,-0.09216,-0.09344,-0.09472,-0.096,-0.09856,-0.09856,-0.09856,-0.09728,-0.09728,-0.09856,-0.09984,-0.09984,-0.09984,-0.09984,-0.09984,-0.10112,-0.10112,-0.10112,-0.10112,-0.1024,-0.1024,-0.10368,-0.1024,-0.10496,-0.10496,-0.10624,-0.10368,-0.10368,-0.1024,-0.10368,-0.10368,-0.10368,-0.10368,-0.10496,-0.10624,-0.10624,-0.10496,-0.10368,-0.10368,-0.10496,-0.10624,-0.10496,-0.10496,-0.10496,-0.10624,-0.1088,-0.1088,-0.1088,-0.1088,-0.1088,-0.10624,-0.10496,-0.10496,-0.10496,-0.10624,-0.10624,-0.10624,-0.10368,-0.10368,-0.10368,-0.10368,-0.1024,-0.10112]
        # Integers
        @test s.nzyear == 1981
        # Strings
        @test s.kevnm == "K8108838"
        @test s.kstnm == "CDV"
        @test s.khole == SAC.sac_cnull
        # Logicals
        @test s.leven == true
    end

    # Writing
    let s1 = SAC.sample()
        @test begin
            write([s1], [tempfile])
            s2 = SAC.read(tempfile)
            s1 == s2
        end
        @test begin
            s1 == SAC.sample()
            write(s1, tempfile)
            s2 = SAC.read(tempfile)
            s1 == s2
        end
    end
        
    let s = SAC.sample()
        s.kevnm = SAC.sac_cnull
        write(s, tempfile)
        d = open(tempfile, "r") do f
            read(f)
        end
        @test String(d[449:(449+15)]) == "-12345  -12345  "
    end

    # Handling spaces in kevnm correctly
    let s = SACtr(0.1, 2), kevnm = "1234567 910"
        s.kevnm = kevnm
        write(s, tempfile)
        s2 = SAC.read(tempfile)
        @test s2.kevnm == kevnm
    end

    # Cutting
    @test begin
        write(SAC.sample(), tempfile)
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
        read(tempfile, String) == "SAC.SACtr(delta=1.0, b=0.0, npts=1)"
    end
    teststring = VERSION < v"0.6-" ?
        "SAC.SACtr[SAC.SACtr(delta=1.0, b=0.0, npts=1),SAC.SACtr(delta=1.0, b=0.0, npts=1)]" :
        VERSION >= v"0.7-" ?
        "SACtr[SAC.SACtr(delta=1.0, b=0.0, npts=1), SAC.SACtr(delta=1.0, b=0.0, npts=1)]" :
        "SAC.SACtr[SAC.SACtr(delta=1.0, b=0.0, npts=1), SAC.SACtr(delta=1.0, b=0.0, npts=1)]"
    @test begin
        open(tempfile, "w") do f
            show(f, [SAC.SACtr(1, 1), SAC.SACtr(1, 1)])
        end
        read(tempfile, String) == teststring
    end
end
