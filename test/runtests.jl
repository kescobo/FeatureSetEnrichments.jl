using FeatureSetEnrichments
using Test

@testset "FeatureSetEnrichments.jl" begin
    @test enrichment_scores(1:10, 1000)[1] == -1
    @test enrichment_scores(991:1000, 1000)[1] == 1

end
