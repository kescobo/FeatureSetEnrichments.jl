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
    pvalue::Float64
    nfeatures::Int
    setranks::Vector{Int}
end

pvalue(result::FSEAResult) = result.pvalue

function Base.show(io::IO, ::MIME"text/plain", fr::FSEAResult)
    println(io, """
    FSEA Result of dataset with
      n features: $(fr.nfeatures)
      n in-set: $(length(fr.setranks))
      p-value: $(fr.pvalue)"""
    )
end

function _run_fsea(::MWU, fset_ranks, nfeatures)
    mwu = MannWhitneyUTest(fset_ranks, 1:nfeatures)
    return FSEAResult(pvalue(mwu), nfeatures, fset_ranks)
end

function _run_fsea(perm::Permutation, fset_ranks, nfeatures)
    es = enrichment_score(fset_ranks, nfeatures)
    ng = length(fset_ranks)
    nperm = perm.nperm

    lt = ThreadsX.count(1:nperm) do _
        rranks = rand(1:nfeatures, ng)
        ses = enrichment_score(rranks, nfeatures)
        es < 0 ? ses < es : ses > es
    end
    return FSEAResult(lt / nperm, nfeatures, fset_ranks)
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
    ranks = invperm(sortperm(vals))
    fset_ranks = ranks[fset_index]
    nfeatures = length(vals)

    _run_fsea(model, fset_ranks, nfeatures)
end


fsea(vals, fset_index) = fsea(MWU(), vals, fset_index)

_es_at_pos(rank, idx, setscore, notscore) = idx * setscore + (rank - idx) * notscore 


function enrichment_score(ranks, nfeatures)
    nr = length(ranks)
    setscore =  -1 / nr
    notscore = 1 / (nfeatures - nr)

    ranks = sort(ranks)
    scores = ThreadsX.map(i-> _es_at_pos(ranks[i], i, setscore, notscore), eachindex(ranks))
    candidates = [scores; scores .- setscore]
    (_, i) = findmax(abs, candidates)
    return candidates[i]
end

"""
    enrichment_score(result::FSEAResult)

Calculates the enrichment score (E.S.) for a feature set
in an [`FSEAResult`](@ref).
See [Mootha et. al. (2013)](https://doi.org/10.1038/ng1180)
for details about how this is calculated.
"""
enrichment_score(result::FSEAResult) = enrichment_score(result.setranks, result.nfeatures)
