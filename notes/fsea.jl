module NeuroFSEA

export getneuroactive,
       fsea,
       enrichment_score

using MultivariateStats
using CodecZlib
using Dictionaries
using InvertedIndices
using HypothesisTests


function get_neuroactive_kos(neuroactivepath=joinpath("/brewster", "kevin", "scratch", "raw_data", "gbm.txt"); consolidate=true)
    neuroactive = Dictionary{String, Vector{String}}()
    desc = ""
    for line in eachline(neuroactivepath)
       line = split(line, r"[\t,]")
       if startswith(line[1], "MGB")
            (mgb, desc) = line
            if consolidate
                desc = rstrip(replace(desc, r"\b[IV]+\b.*$"=>""))
                desc = replace(desc, r" \([\w\s\-]+\)"=>"")
                desc = replace(desc, r"^.+ \(([\w\-]+)\) (.+)$"=>s"\1 \2")
                desc = replace(desc, " (AA"=>"")
            end
            @info "getting unirefs for $desc"
            !in(desc, keys(neuroactive)) && insert!(neuroactive, desc, String[])
       else
           filter!(l-> occursin(r"^K\d+$", l), line)
           append!(neuroactive[desc], String.(line))
       end
   end
   return neuroactive
end

function getneuroactive(features; neuroactivepath=joinpath("/brewster", "kevin", "scratch", "raw_data", "gbm.txt"),
                                  map_ko_uniref_path=joinpath("/brewster", "kevin", "scratch", "raw_data", "map_ko_uniref90.txt.gz"), 
                                  consolidate=true)
    @warn "this was changed"
    neuroactivekos = get_neuroactive_kos(neuroactivepath; consolidate)

    kos2uniref = Dictionary{String, Vector{String}}()
    for line in eachline(GzipDecompressorStream(open(map_ko_uniref_path)))
        line = split(line, '\t')
        insert!(kos2uniref, line[1], map(x-> String(match(r"UniRef90_(\w+)", x).captures[1]), line[2:end]))
    end

    neuroactive_index = Dictionary{String, Vector{Int}}()
    for na in keys(neuroactivekos)
        searchfor = Iterators.flatten([kos2uniref[ko] for ko in neuroactivekos[na] if ko in keys(kos2uniref)]) |> Set
        pos = findall(u-> u in searchfor, features)
        insert!(neuroactive_index, na, pos)
    end
    for k in keys(neuroactive_index)
        unique!(neuroactive_index[k])
    end
    return neuroactive_index
end

fsea(cors, pos) = (cors, pos, MannWhitneyUTest(cors[pos], cors[Not(pos)]))

function fsea(cors, allfeatures::AbstractVector, searchset::Set)
    pos = findall(x-> x in searchset, allfeatures)

    return fsea(cors, pos)
end

function fsea(occ::AbstractMatrix, metadatum::AbstractVector, pos::AbstractVector{<:Int})
    cors = let notmissing = map(!ismissing, metadatum)
        occ = occ[:, notmissing]
        metadatum = metadatum[notmissing]
        cor(metadatum, occ, dims=2)'
    end

    return fsea(cors, pos)
end

function enrichment_score(setcors, notcors)
    srt = sortperm([setcors; notcors]; rev=true)
    ranks = invperm(srt)
    setranks = Set(ranks[1:length(setcors)])
    
    setscore =  1 / length(setcors)
    notscore = -1 / length(notcors)
    
    ys = cumsum(i ∈ setranks ? setscore : notscore for i in eachindex(ranks))
    
    lower, upper = extrema(ys)
    return abs(lower) > abs(upper) ? lower : upper
end

end #module NeuroFSEA