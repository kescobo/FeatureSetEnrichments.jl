function get_neuroactive_kos(;neuroactivepath=joinpath("/brewster", "kevin", "scratch", "raw_data", "gbm.txt"), consolidate=true)
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

function get_neuroactive_unirefs(; neuroactivepath=joinpath("/brewster", "kevin", "scratch", "raw_data", "gbm.txt"),
                                    map_ko_uniref_path=joinpath("/brewster", "kevin", "scratch", "raw_data", "map_ko_uniref90.txt.gz"),
                                    consolidate=true)
    kos2uniref = Dictionary{String, Vector{String}}()
    for line in eachline(GzipDecompressorStream(open(map_ko_uniref_path)))
        line = split(line, '\t')
        insert!(kos2uniref, line[1], map(x-> String(match(r"UniRef90_(\w+)", x).captures[1]), line[2:end]))
    end

    nakos = get_neuroactive_kos(; neuroactivepath, consolidate)
    naunirefs = Dictionary{String, Set{String}}()
    for (key, kos) in pairs(nakos)
        set!(naunirefs, key, Set(Iterators.filter(!isnothing, Iterators.flatten(get(kos2uniref, ko, [nothing]) for ko in kos))))
    end
    naunirefs
end
