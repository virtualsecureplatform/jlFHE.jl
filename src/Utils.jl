function dtot32(d)
    return map(x->round(UInt32,x),(d.%1)*2.0^32-(sign.(d)- sign.(d.^2))*2.0^31)
end
function gaussian32(μ,α)
    return dtot32(α*randn(length(μ)) .+ μ)
end