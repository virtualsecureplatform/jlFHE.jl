#Reference https://math.stackexchange.com/questions/1435448/negacyclic-fft-multiplication

module MulFFT
    using FFTW
    struct TwistFFTPlan
        fftplan::FFTW.Plan
        ifftplan::FFTW.Plan
        twist::Array{ComplexF64,1}
    end
    function TwistFFTPlan(N::Int32)
        TwistFFTPlan(FFTW.plan_fft(Float64.(rand(-2^9:2^9,div(N,2)))+im*rand(-2^9:2^9,div(N,2))),FFTW.plan_ifft(Float64.(rand(-2^9:2^9,div(N,2)))+im*rand(-2^9:2^9,div(N,2))),TwistGen(N))
    end
    function TwistGen(N::Int32)
        return [exp(im*k*pi/N) for k in 0:div(N,2)-1]
    end

    function TwistFFT(a::Array{Int32,1},fftplan::FFTW.Plan,twist::Array{Complex{Float64}})
        return fftplan*(complex.(Float64.(a[1:div(length(a),2)]),Float64.(a[div(length(a),2)+1:length(a)])).* twist)
    end
    function TwistIFFT(a::Array{Complex{Float64},1},ifftplan::FFTW.Plan,twist::Array{Complex{Float64}})
        b = (ifftplan*a).* conj(twist)
        return append!(real(b),imag(b))
    end

    function PolyMul(a::Array{Int32,1},b::Array{Int32,1},TwistFFTPlan::TwistFFTPlan)
        return map(x->UInt32(mod(round(Int64,x),Int64(2)^32)),TwistIFFT(TwistFFT(a,TwistFFTPlan.fftplan,TwistFFTPlan.twist).*TwistFFT(b,TwistFFTPlan.fftplan,TwistFFTPlan.twist),TwistFFTPlan.ifftplan,TwistFFTPlan.twist))
    end
end