module TRLWE
    using ..MulFFT
    using ..Utils
    function trlweSymEncrypt(p::Array{Float64,1},α::Float64,key::Array{UInt32,1},plan::MulFFT.TwistFFTPlan)
        a::Array{UInt32,1} = rand(0:1,length(key))
        b::Array{UInt32,1} = Utils.gaussian32(p,α) + MulFFT.PolyMul(signed.(a),signed.(key),plan)
        return [a,b]
    end
    function trlweSymDecrypt(c::Array{Array{UInt32,1},1},key::Array{UInt32,1},plan::MulFFT.TwistFFTPlan)
        return -((c[2] - MulFFT.PolyMul(signed.(c[1]),signed.(key),plan)).>UInt32(2)^31) .+ 1
    end

    function SampleExtractIndex(trlwe::Array{Array{UInt32,1},1},index::Int32)
        N = length(trlwe[1])
        return vcat([trlwe[1][index - i] for i in 0:index-1],[-trlwe[1][N+index-i] for i in index:N-1],trlwe[2][index])
    end
end