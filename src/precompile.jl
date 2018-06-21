"""
    _precompile()

Precompile a selection of slow-to-JIT functions to help speed up the first
call of common SAC operations.
"""
function _precompile()
    float_types = (Float64, Float32)
    number_types = (Float64, Float32, Int)
    # filtering.jl
    for f1 in number_types
        precompile(highpass!, (SACtr, f1))
        precompile(lowpass!, (SACtr, f1))
        precompile(apply_filter!, (SACtr, f1, Int))
        for f2 in number_types
            precompile(bandpass!, (SACtr, f1, f2))
        end
    end
    precompile(get_filter_prototype, (Symbol, Int))
   
    # great_circle.jl
    for f1 in number_types, f2 in number_types, f3 in number_types, f4 in number_types
        precompile(_great_circle, (f1, f2, f3, f4))
        for f5 in number_types
            precompile(_great_circle, (f1, f2, f3, f4, f5))
        end
    end
    precompile(update_great_circle!, (SACtr,))
    for f1 in number_types, f2 in number_types
        precompile(vincentydist, (f1, f2, f2, f2, f2, f2))
    end
    
    # io.jl
    precompile(SAC.read, (String,))
    for f1 in number_types, f2 in number_types
        precompile(read_cut, (String, Symbol, f1, Symbol, f2))
    end
    precompile(file_is_native_endian, (String,))
    for f in number_types
        precompile(_write_swap, (Bool, IOStream, f))
    end
    precompile(_write_swap, (IOStream, String, Int))
    precompile(write, (SACtr, String))
    precompile(write, (Array{SACtr}, Array{String}))
    precompile(read_wild, (String, String))
    precompile(show, (IO, SACtr))
    precompile(show, (IO, MIME"text/plain", SACtr))
    
    # operations.jl
    for f1 in number_types
        precompile(add!, (SACtr, f1))
        precompile(add!, (Array{SACtr}, f1))
        precompile(divide!, (SACtr, f1))
        precompile(divide!, (Array{SACtr}, f1))
        precompile(multiply!, (SACtr, f1))
        precompile(multiply!, (Array{SACtr}, f1))
        precompile(normalise!, (SACtr,))
        for f2 in number_types
            precompile(cut!, (SACtr, f1, f2))
            precompile(cut!, (Array{SACtr}, f1, f2))
            precompile(cut!, (Array{SACtr}, Array{f1}, Array{f2}))
        end
    end
    precompile(differentiate!, (SACtr, Int))
    for t in (SACtr, Array{SACtr})
        precompile(envelope!, (t,))
        precompile(fft, (t,))
        precompile(rmean!, (t,))
        precompile(rtrend!, (t,))
        precompile(taper!, (t,))
        precompile(update_headers!, (t,))
    end
    precompile(integrate!, (SACtr, Symbol))
    precompile(interpolate!, (SACtr,))
    
    # types.jl
    precompile(SACtr, (Vector{UInt8},))
    precompile(SACtr, (Float64, Int, Float64))
    precompile(getindex, (SACtr, Symbol))
    precompile(getindex, (Array{SACtr}, Symbol))
    for f in [float_types... String]
        precompile(setindex!, (SACtr, f, Symbol))
        precompile(setindex!, (Array{SACtr}, f, Symbol))
    end
    precompile(==, (SACtr, SACtr))
    precompile(copy, (SACtr,))
    
    # util.jl
    precompile(sample, ())
    for t in ((SACFloat,), (SACInt,), (SACChar,), (SACBool,), (SACtr, Symbol))
        precompile(isundefined, t)
    end
end

_precompile()
