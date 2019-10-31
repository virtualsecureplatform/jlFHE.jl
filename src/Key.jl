include("lweParams.jl")
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