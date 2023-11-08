using FeatureSetEnrichments
using Test

@testset "FeatureSetEnrichments.jl" begin
    @test enrichment_score(1:10, 1000) == -1
    @test enrichment_score(991:1000, 1000) == 1

end
