using jlFHE
using Polynomials
using Test

n = 500
N = 1024
Bgbit = 10
α= 2.44e-5
plan = jlFHE.MulFFT.TwistFFTPlan(N)
@testset "jlFHE.jl" begin
    # Write your own tests here.
    testnum = 1000
    sk = jlFHE.Key.SecretKey(n,N,α,2,Bgbit,2^(-24),8,2,2^(-24))
    p = rand(0:1,1000)
    c = jlFHE.TLWE.bootsSymEncrypt(p,sk)
    @test jlFHE.TLWE.bootsSymDecrypt(c,sk) == p

    a = [rand(-2^31:2^31-1,N) for i in 1:testnum]
    b = [rand(-2^(Bgbit-1):2^(Bgbit-1)-1,N) for i in 1:testnum]
    y = [jlFHE.MulFFT.PolyMul(a[i],b[i],plan) for i in 1:testnum]
    z = jlFHE.Utils.dtot32.(coeffs.(map(x->rem(x,Poly(append!(append!([1.0],zeros(N-1)),1.0))),(Poly.(a)*2.0^-32) .* Poly.(b))))
    @test y == z

    p = [rand(0:1,N) for i in 1:testnum]
    c = map(x->jlFHE.TRLWE.trlweSymEncrypt(x,α,sk.key.trlwe,plan),[(2*p[i] .- 1)*jlFHE.Key.μ for i in 1:testnum])
    @test p == map(x->jlFHE.TRLWE.trlweSymDecrypt(x,sk.key.trlwe,plan),c)

    @test p == map(x->[jlFHE.TLWE.tlweSymDecrypt(jlFHE.TRLWE.SampleExtractIndex(x,i),sk.key.trlwe) for i in 1:N],c)
end
