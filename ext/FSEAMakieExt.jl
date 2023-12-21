module FSEAMakieExt

using FeatureSetEnrichments:FSEAResult, enrichment_score
using Makie

const _es_at_pos = FeatureSetEnrichments._es_at_pos

# function plot_fsea(setcors, notcors; label="", ylabel="enrichment", kwargs...)
#     fig = Figure()
#     grid, ax1, ax2 = plot_fsea!(fig.layout, setcors, notcors; label, ylabel, kwargs...)
#     fig, ax1, ax2
# end

function _gen_pts(setranks, nfeatures, setscore, notscore)
    pts = Iterators.map(eachindex(setranks)) do i
        rank = setranks[i]
        prevrank = rank - 1
        score = _es_at_pos(rank, i, setscore, notscore)
        prevscore = score - setscore
        return ((prevrank, prevscore), (rank, score))
    end
    return [(0,0); sort(collect(Iterators.flatten(pts)); lt = (x,y)->x[1]<y[1]); (nfeatures+2,0)]
end

function Makie.plot!(grid::Makie.GridLayout, fr::FSEAResult; title="", ylabel="", linecolor=:gray, kwargs...)
    nfeatures = fr.nfeatures
    setranks = sort(fr.setranks)
    nr = length(setranks)

    setscore =  -1 / nr
    notscore = 1 / (nfeatures - nr)
    
    pts = _gen_pts(setranks, nfeatures, setscore, notscore)
    
    t = "ES: $(round(enrichment_score(fr), digits=3))"

    ax1 = Axis(grid[1,1]; title, ylabel=isempty(ylabel) ? t : ylabel, kwargs...)
    hidexdecorations!(ax1)

    ax2 = Axis(grid[2,1])
    hidedecorations!(ax2)
    
    lines!(ax1, pts; color=linecolor)
    vlines!(ax2, setranks; color=:black)
    
    rowsize!(grid, 2, Relative(1/8))
    rowgap!(grid, 1, Fixed(0))

    linkxaxes!(ax1, ax2)
    tightlimits!.((ax1, ax2))
    grid, ax1, ax2
end

function Makie.plot(fr::FSEAResult; title="", ylabel="", linecolor=:gray, kwargs...)
    fig = Figure()
    grid = GridLayout(fig[1,1])
    (gr, ax1, ax2) = plot!(grid, fr; title, ylabel, linecolor, kwargs...)
    return fig
end

# function plot_corrband!(ax, cors; bandres=5000)
#     ax.xlabel="rank"
#     ax.ylabel="correlation"
#     srt = sortperm(cors)
#     nfeatures = length(cors)
#     rn = nfeatures > bandres ? round.(Int, range(1, nfeatures; length=bandres)) : range(1, nfeatures)
#     blow, bup = extrema(cors)

#     xs = 1:nfeatures
#     band!(ax, xs[rn], fill(blow, length(rn)), fill(bup, length(rn)); color=cors[srt[rn]])
    
#     lower = [x < 0 ? x : 0.0 for x in cors[srt]]
#     upper = [x > 0 ? x : 0.0 for x in cors[srt]]
#     band!(ax, rn, lower[rn], upper[rn]; color=:lightgray)
#     tightlimits!(ax)
# end

end