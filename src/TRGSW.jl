module TRGSW
    using ..MulFFT
    using ..TRLWE
    using ..Key
    using LinearAlgebra

    function trgswSymEncrypt(p::Array{UInt32,1},α::Float64,h::Array{Float64},key::Array{UInt32,1},plan::MulFFT.TwistFFTPlan)
        l = length(h)
        c = [[TRLWE.trlweSymEncrypt(zeros(length(key)),α,key,plan) for i in 1:2*l]...]
        ph = map(x->round.(UInt32,x),p * h')
        for i in 1:l
            c[i][1] .+= ph[i]
            c[i+l][1] .+= ph[i]
        end
        return  c
    end
    function trgswfftSymEncrypt(p::Array{UInt32,1},α::Float64,h::Array{Float64},key::Array{UInt32,1},plan::MulFFT.TwistFFTPlan)
        trgsw = trgswSymEncrypt(p,α,h,key,plan)
        return hcat([MulFFT.TwistFFT.(trgsw[i],plan.fftplan,plan.twist) for i in 1:2*length(h)]...)
    end

    function Decomposition(trlwe::Array{Array{UInt32,1},1},params::Key.lweParams)
        temp = map(x->map(y->round(Int32,(y.%(2^params.Bgbit))-2^(params.Bgbit-1)),(x.+params.offset) * params.decb'),trlwe)
        return vcat([[temp[i][j,:] for j in 1:2] for i in 1:params.l]...)
    end

    function trgswExternalProduct(trgsw::Array{Array{Array{UInt32,1},1},1},trlwe::Array{Array{UInt32,1},1},params::Key.lweParams,plan::MulFFT.TwistFFTPlan)
        decvec = Decomposition(trlwe,params)
        return [sum([MulFFT.PolyMul(decvec[i],signed.(trgsw[i][j]),plan) for i in 1:2*params.l]) for j in 1:2]
    end
end