using FFTW
struct TwistFFTPlan
    fftplan::FFTW.Plan
    ifftplan::FFTW.Plan
    twist::Array{ComplexF64}
end
function TwistFFTPlan(N)
    TwistFFTPlan(FFTW.plan_fft(Float64.(rand(-2^9:2^9,N/2))+im*rand(-2^9:2^9,N/2)),FFTW.plan_ifft(Float64.(rand(-2^9:2^9,N/2))+im*rand(-2^9:2^9,N/2)),TwistGen(N))
end
function TwistGen(N::Int32)
    return [exp(1im*k*pi/N) for k in 1:N/2]
end

function TwistFFT(a::Array{Int32,1},fftplan::FFTW.Plan,twist::Array{ComplexF64})
    return fftplan((a[1:length(a)/2] + im*a[length(a)/2+1:length(a)]).*twist)
end
function TwistIFFT(a::Array{ComplexF64,1},ifftplan::FFTW.Plan,twist::Array{ComplexF64})
    return ifftplan(a).*conj(twist)
end

function PolyMul(a::Array{UInt32,1},b::Array{UInt32,1},TwistFFTPlan::TwistFFTPlan)
    unsigned(round(Int32,TwistIFFT(TwistFFT(a,TwistFFTPlan.fftplan,TwistFFTPlan.twist).*TwistFFT(b,TwistFFTPlan.fftplan,TwistFFTPlan.twist),TwistFFTPlan.ifftplan,TwistFFTPlan.twist)%2^32))
end
