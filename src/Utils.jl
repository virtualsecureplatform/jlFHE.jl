module Utils
    function dtot32(d)
        return map(x->UInt32(mod(round(Int64,mod(x,1.0)*2.0^32),Int64(2)^32)),d)
    end
    function gaussian32(μ,α)
        return dtot32(α*randn(length(μ)) .+ μ)
    end
end