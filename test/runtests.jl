using AsteroidPipeline
using Test

@testset "AsteroidPipeline.jl" begin
    @test_throws ErrorException detect_sources(zeros(4, 4); threshold=5.0)
    @test_throws ErrorException link_candidates([], [])
    @test_throws ErrorException crossmatch_catalog([], :skybot; radius=5.0)
end
