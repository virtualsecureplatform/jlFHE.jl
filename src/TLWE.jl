module TLWE
    using ..Key
    using ..Utils
    using LinearAlgebra
    function tlweSymEncrypt(p::Float64,α::Float64,key::Array{UInt32,1})
        a::Array{UInt32,1} = rand(0:1,length(key))
        b::UInt32 = Utils.gaussian32(p,α)[1] + a ⋅ key
        return append!(a,b)
    end
    function tlweSymDecrypt(c::Array{UInt32,1},key::Array{UInt32,1})
        return 1 - ((c[end] - c[1:end-1] ⋅ key) > UInt32(2)^31)
    end

    function bootsSymEncrypt(p::Array{Int32,1},sk::Key.SecretKey)
        return [tlweSymEncrypt(x, sk.params.α, sk.key.tlwe) for x in (2*p .- 1)*Key.μ]
    end
    function bootsSymDecrypt(c::Array{Array{UInt32,1},1},sk::Key.SecretKey)
        return [tlweSymDecrypt(x,sk.key.tlwe) for x in c]
    end
end