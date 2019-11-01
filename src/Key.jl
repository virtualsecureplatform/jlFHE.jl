module Key
    struct lweParams
        n::Int32
        N::Int32
        l::Int32
        Bgbit::Int32
        t::Int32
        basebit::Int32
        α::Float64
        αᵦₖ::Float64
        αₖₛ::Float64
        offset::Int32
        h::Array{Float64}
        decb::Array{Float64}
    end
    function lweParams(n::Int32,N::Int32,α::Float64,l::Int32,Bgbit::Int32,αᵦₖ::Float64,t::Int32,basebit::Int32,αₖₛ::Float64)
        Bg::Float64 = 2.0 ^ Bgbit
        h::Array{Float64} = [Bg^(-(i+1)) for i in 1:l]
        decb::Array{Float64} = [2.0^(-32) * Bg^(i+1) for i in 1:l]
        offset::UInt32 = Bg/2 * sum(2^32 * h)
        return lweParams(n,N,l,Bgbit,t,basebit,α,αᵦₖ,αₖₛ,offset,h,decb)
    end

    struct lweKey
        tlwe::Array{UInt32,1}
        trlwe::Array{UInt32,1}
    end
    function lweKey(n::Int32,N::Int32)
        return lweKey(rand(0:1,n),rand(0:1,N))
    end

    struct SecretKey
        params::lweParams
        key::lweKey
    end
    function SecretKey(n::Int32,N::Int32,alpha::Float64,l::Int32,Bgbit::Int32,bkalpha::Float64,t::Int32,basebit::Int32,ksalpha::Float64)
        return SecretKey(lweParams(n,N,alpha,l,Bgbit,bkalpha,t,basebit,ksalpha),lweKey(n,N))
    end
    const μ = 2.0^-3
end