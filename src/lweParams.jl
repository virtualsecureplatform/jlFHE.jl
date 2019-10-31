include("MulFFT.jl")

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
μ = 2.0^-3