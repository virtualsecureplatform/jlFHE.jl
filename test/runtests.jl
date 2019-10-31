using jlFHE
using Test
include("TLWEenc.jl")

@testset "jlFHE.jl" begin
    # Write your own tests here.
    sk = jlFHE.SecretKey(500,1024,2^(-9),2,10,2^(-24),8,2,2^(-24))
    p = rand(0:1,100)
    c = jlFHE.bootsSymEncrypt(p,sk)
    @test jlFHE.bootsSymDecrypt(c,sk) == p
end
