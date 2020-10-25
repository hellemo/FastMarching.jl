using FastMarching
using Test
using LinearAlgebra

@testset "Test FastMarching" begin
    function testI(n::Integer, endpoint=[1.0, 1.0], stepsize=0.1, T=Float64)
        speedmap = zeros(T, (n, n)) + I .* T(1000) .+ T(0.1)
        source = T.([float(size(speedmap, 1)), float(size(speedmap, 2))])
        distancemap = FastMarching.msfm(speedmap, source, true, true)
        endpoint = T.([1.0, 1.0])
        stepsize = T(stepsize)
        return shortestline = FastMarching.shortestpath(
            distancemap, endpoint, source, stepsize
        )
    end

    endpoint = [1.0, 1.0]
    @testset "testeye" begin
        @testset "size: $t stepsize: $stepsize" for t in 2:2:10,
stepsize in [0.1, 0.01],
T in [Float16, Float32, Float64]

            shortestline = testI(t, endpoint, stepsize, T)
            @test sum(shortestline[:, 2] - shortestline[:, 1]) == 0
            @test sqrt((shortestline[1, 1] - 1.0)^2 - (shortestline[1, 2] - 1.0)^2) <=
                  stepsize
            @test eltype(shortestline) == T
        end
    end

    @testset "correctness" begin
        for N in [2, 5, 10]
            speedimage = ones(N, N)
            sourcepoint = [float(N); float(N)]
            distancemap = FastMarching.msfm(speedimage, sourcepoint, true, true)
            @test distancemap[1, end] ≈ (N - 1)
            @test distancemap[end, 1] ≈ (N - 1)
            @test distancemap[end] ≈ 0
            @test isapprox(distancemap[1], sqrt(2 * (N - 1)^2); atol=0.2)
        end
    end

    include(joinpath(@__DIR__, "multiplestarts.jl"))
end
