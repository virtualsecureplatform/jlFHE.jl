module TRGSW
    using ..MulFFT
    using ..TRLWE
    using ..Key
    using LinearAlgebra

    function trgswSymEncrypt(p::Array{UInt32,1},α::Float64,h::Array{Float64},key::Array{UInt32,1},plan::MulFFT.TwistFFTPlan)
        l = length(h)
        c = vcat([TRLWE.trlweSymEncrypt(zeros(length(key)),α,key,plan) for i in 1:2*l]'...)
        ph = p * h'
        for i in 1:l
            c[i,1] += ph[i]
            c[i+l,1] += ph[i]
        end
        return  c
    end
    function trgswfftSymEncrypt(p::Array{UInt32,1},α::Float64,h::Array{Float64},key::Array{UInt32,1},plan::MulFFT.TwistFFTPlan)
        trgsw = trgswSymEncrypt(p,α,h,key,plan)
        return hcat([MulFFT.TwistFFT.(trgsw[:,i],plan.fftplan,plan.twist) for i in 1:2*length(h)]...)
    end

    function Decomposition(trlwe::Array{UInt32,2},params::Key.lweParams)
        return vcat(map(x->map(y->unsigned.(round(Int32,(y.%params.Bg)-params.Bg/2)),(x.+params.offset) * params.decb'),trlwe)...)
    end

    function trgseExternalProduct(trgsw::Array{LinearAlgebra.Adjoint{UInt32,Array{UInt32,1}},2},trlwe::Array{UInt32,2},params::Key.lweParams)
        decvec = Decomposition(trlwe,params)
        return [sum([MulFFT.PoluMul(decvec[i],trgsw[i,j],plan) for i in 1:2*params.l]) for j in 1:2]
    end
end