module FeatureSetEnrichments

using HypothesisTests: MannWhitneyUTest
using InvertedIndices: Not
using ThreadsX: count

abstract type FSEATest end

"""
    Permutation(nperm::Int=1000)

A permutation-type fsea as described in
[Mootha et. al. (2013)](https://doi.org/10.1038/ng1180).
First, the [`enrichment_score`](@ref) (E.S.) is calculated
for the true feature set.
Second, the random indices of the same length
as the feature set are selected and the E.S.
is recalculated `nperm` times.

The pseudo-pvalue is the fraction of permutations
where the absolute value of the E.S.
was greater than that of the true feature set.
"""
struct Permutation <: FSEATest
    nperm::Int
end
Permutation() = Permutation(1000)

"""
    MWU()

Uses the Mann-Whitney U-test,
which tests the hypothesis that
a random value pulled from the feature set
will be consistently less than or consistently greater than
a random value pulled from the out-of-set values,
against the null hypothesis that feature set values
are equally likely to be less than or greater than
those not in the feature set.
"""
struct MWU <: FSEATest end


struct FSEAResult
    es::Float64
    pvalue::Float64
    nfeatures::Int
    setranks::Vector{Int}
end

"""
    fsea(model::FSEAModel, vals, fset_idx)

Perform a feature set enrichment analysis (FSEA)
on `vals`, where the feature set being tested
has indices `fset_idx`.

There are currently two types of tests that can be performed:

1. [`Permutation`](@ref) which uses the procedure originally published
   in [Mootha et. al. (2013)](https://doi.org/10.1038/ng1180)
2. [`MWU`](@ref) (the default) which uses the non-parametric
   Mann-Whitney U-test.

"""
function fsea(model::FSEATest, vals, fset_index)
    srtd = zeros(Int, length(vals))
    sortperm!(srtd, vals; rev=true)
    invper
end


fsea(cors, pos) = (cors, pos, MannWhitneyUTest(cors[pos], cors[Not(pos)]))

function fsea(cors, allfeatures::AbstractVector, searchset::Set)
    pos = ThreadsX.findall(x-> x in searchset, allfeatures)

    return fsea(cors, pos)
end


function enrichment_score(cors, pos)
    srt = sortperm(cors; rev=true)
    ranks = invperm(srt)
    setranks = ranks[pos]    
    ng = length(pos)

    setscore =  1 / ng
    notscore = -1 / (length(cors) - ng)

    ys = cumsum(i ∈ setranks ? setscore : notscore for i in eachindex(ranks))
    lower, upper = extrema(ys)
    return abs(lower) > abs(upper) ? lower : upper
end

function fsea_permute(cors, pos; nperm=1000)
    es = enrichment_score(cors, pos)
    ng = length(pos)
    lt = ThreadsX.count(1:nperm) do _
        spos = rand(eachindex(cors), ng)
        ses = enrichment_score(cors, spos)
        es < 0 ? ses < es : ses > es
    end
    return (lt / nperm, es)
end

end
