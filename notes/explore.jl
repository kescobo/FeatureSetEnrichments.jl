using VKCComputing
using BiobakeryUtils
using DataFrames
using CSV

eegbase = CSV.read("data/eeg_3mo_baseline.csv", DataFrame; stringtype=String)
rename!(eegbase, "subject_ID"=> "subject_id")
eegvep = CSV.read("data/eeg_3mo_vep.csv", DataFrame; stringtype=String)
metadata = CSV.read("data/Khula_MicrobiomeMetaData_SA_25052023.csv", DataFrame; stringtype=String)
subset!(metadata, "zymo_code_3m"=> ByRow(!ismissing))
metadata.zymo_code_3m = replace.(metadata.zymo_code_3m, r"z3m[o0]"i => "Z3MO")
select!(metadata, "subject_id", "zymo_code_3m"=> "zymo", "age_zymo_3m_wks"=> "age")
# If new environment, requires running `VKCComputing.set_readonly_pat!()` and `VKCComputing.set_airtable_dir!()`
base = LocalBase()

# `Number_Segs_Post-Seg_Rej` = number of retained trials

#-


#-

seqs = let df = DataFrame()
    subs = filter(!ismissing, map(s-> get(base["Subjects"], "khula-$s", missing), metadata.subject_id))
    for sub in subs
       for biosp in get(sub, :Biospecimens, [])
           biosp = base[biosp]
           !haskey(biosp, :visit) && continue
           for seq in get(biosp, :seqprep, [])
                seq = base[seq]
                seq[:keep] != 1 && continue
                push!(df, (; subject_ref = sub[:uid],
                            subject_id  = replace(sub[:uid], "khula-"=>""),
                            biospecimen = biosp[:uid],
                            visit       = first(base[biosp[:visit]])[:uid],
                            seqprep     = seq[:uid]
                            )
               )
           end
       end
   end
   df
end

gfdf = let
    gffiles = filter(readdir("/grace/sequencing/processed/mgx/humann/main", join=true)) do path
        file = basename(path)
        m = match(r"(SEQ\d+)_S\d+_genefamilies.tsv", file)
        isnothing(m) && return false
        return m[1] in seqs.seqprep
    end

    DataFrame(map(gffiles) do f
        m = match(r"(SEQ\d+)_(S\d+)_genefamilies.tsv", f)
        (seqprep, S_well) = m.captures
        dir = dirname(f)
        file = basename(f)
        (; seqprep, S_well, dir, file, path=f)
    end)
end

leftjoin!(seqs, gfdf; on="seqprep")

#-

leftjoin!(eegbase, subset(seqs, "visit"=> ByRow(==("3mo"))); on="subject_id")
subset!(eegbase, "seqprep"=> ByRow(!ismissing))

leftjoin!(eegvep, subset(seqs, "visit"=> ByRow(==("3mo"))); on="subject_id")
subset!(eegvep, "seqprep"=> ByRow(!ismissing))

#-

prep2sub = Dict(p => s for (p,s) in zip(seqs.seqprep, seqs.subject_id))

gfs = mapreduce(vcat, subset(eegvep, "file"=> ByRow(!ismissing)).path) do f
    seqprep = replace(basename(f), r"_S\d+_genefamilies\.tsv"=> "")
    df = CSV.read(f, DataFrame; skipto = 2, header=["feature", "abundance"])
    df.seqprep .= seqprep
    df.subject .= prep2sub[seqprep]
    subset!(df, "feature"=> ByRow(f-> !contains(f, '|')))
end

gfs_wide = unstack(gfs, "feature", "abundance")

foreach(n-> gfs_wide[!,n] = coalesce.(gfs_wide[!,n], 0.), names(gfs_wide)[3:end])
leftjoin!(gfs_wide, 
          select(eegvep, "subject_id"=> "subject", 
                         "child_sex_0female_1male"=>"sex",
                         "age_3m_eeg"=> "age_weeks",
                         "Number_Segs_Post-Seg_Rej" => "trials",
                         "peak_amp_N1",
                         "peak_latency_N1",
                         "peak_amp_P1",
                         "peak_latency_P1",
                         "peak_amp_N2",
                         "peak_latency_N2"
                );
    on="subject"
)

select!(gfs_wide, "subject", "seqprep","sex","age_weeks", "trials",
                  "peak_amp_N1","peak_latency_N1",
                  "peak_amp_P1","peak_latency_P1",
                  "peak_amp_N2","peak_latency_N2", 
                  Cols(Not("UNMAPPED"))
)

#-

using Statistics
using ThreadsX
include("fsea.jl")
using .NeuroFSEA
using MultipleTesting
using GLM
using HypothesisTests

#-

function runlms(indf, outfile, respcol, featurecols;
                formula = term(:func) ~ term(respcol) + term(:age_weeks) + term(:trials),
        )
        @debug "Respcol: $respcol"
    lmresults = DataFrame(ThreadsX.map(featurecols) do feature
        # @debug "Feature: $feature"

        # ab = collect(indf[!, feature] .+ (minimum(indf[over0, feature])) / 2) # add half-minimum non-zerovalue
        df = select(indf, respcol, "age_weeks", "trials")
        over0 = indf[!, feature] .> 0
        
        df.func = over0
        # @debug "DataFrame: $df"

        mod = glm(formula, df, Binomial(), LogitLink(); dropcollinear=false)
        ct = DataFrame(coeftable(mod))
        ct.feature .= feature
        rename!(ct, "Pr(>|z|)"=>"pvalue", "Lower 95%"=> "lower_95", "Upper 95%"=> "upper_95", "Coef."=> "coef", "Std. Error"=>"std_err")
        select!(ct, Cols(:feature, :Name, :))
        return NamedTuple(only(filter(row-> row.Name == respcol, eachrow(ct))))
    end)

    subset!(lmresults, "pvalue"=> ByRow(!isnan))
    DataFrames.transform!(lmresults, :pvalue => (col-> MultipleTesting.adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(lmresults, :qvalue)

    CSV.write(outfile, lmresults)
    lmresults
end

function runcors(indf, outfile, respcol, featurecols)
    corresults = vec(cor(Matrix(indf[!, featurecols]), indf[!, respcol]; dims=1))
    nact = getneuroactive(replace.(featurecols, "UniRef90_"=>""))
    
    fseas = DataFrame(ThreadsX.map(collect(keys(nact))) do gs
        ixs = nact[gs]
        isempty(ixs) && return (; cortest=respcol, geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = corresults[ixs]
        isempty(cs) && return (; cortest=respcol, geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = corresults[Not(ixs)]
        mwu = MannWhitneyUTest(cs, acs)
        es = enrichment_score(cs, acs)

        return (; cortest=respcol, geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)
    subset!(fseas, "pvalue"=> ByRow(!isnan))
    fseas.qvalue = adjust(fseas.pvalue, BenjaminiHochberg())
    sort!(fseas, :qvalue)
    CSV.write(outfile, fseas)
    fseas
end

#-


res = let nact = getneuroactive(replace.(names(gfs_wide, r"UniRef"), "UniRef90_"=>""))
    mapreduce(vcat, ["peak_amp_N1","peak_latency_N1","peak_amp_P1","peak_latency_P1","peak_amp_N2","peak_latency_N2"]) do eeg_feat 
        res = runlms(gfs_wide, "data/outputs/$(eeg_feat)_lms.csv", eeg_feat, names(gfs_wide, r"UniRef"))

        fseas = DataFrame(ThreadsX.map(collect(keys(nact))) do gs
            cortest = first(res.Name)
            zs = res.z
            
            ixs = nact[gs]
            isempty(ixs) && return (; cortest, geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

            cs = zs[ixs]
            isempty(cs) && return (; cortest, geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

            acs = zs[Not(ixs)]
            mwu = MannWhitneyUTest(cs, acs)
            es = enrichment_score(cs, acs)

            return (; cortest, geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
        end)

        subset!(fseas, :pvalue=> ByRow(!isnan))
        fseas.qvalue = adjust(fseas.pvalue, BenjaminiHochberg())
        sort!(fseas, :qvalue)
        CSV.write("data/outputs/$(eeg_feat)_gsea.csv", fseas)
        fseas
    end
end

res_noage = let nact = getneuroactive(replace.(names(gfs_wide, r"UniRef"), "UniRef90_"=>""))
    mapreduce(vcat, ["peak_amp_N1","peak_latency_N1","peak_amp_P1","peak_latency_P1","peak_amp_N2","peak_latency_N2"]) do eeg_feat

        @debug "EEG feat: $eeg_feat"
        res = runlms(gfs_wide, "data/outputs/$(eeg_feat)_noage_lms.csv", eeg_feat, names(gfs_wide, r"UniRef"); formula = term(:func) ~ term(eeg_feat))

        fseas = DataFrame(ThreadsX.map(collect(keys(nact))) do gs
            cortest = first(res.Name)
            zs = res.z
            
            ixs = nact[gs]
            isempty(ixs) && return (; cortest, geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

            cs = zs[ixs]
            isempty(cs) && return (; cortest, geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

            acs = zs[Not(ixs)]
            mwu = MannWhitneyUTest(cs, acs)
            es = enrichment_score(cs, acs)

            return (; cortest, geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
        end)

        subset!(fseas, :pvalue=> ByRow(!isnan))
        fseas.qvalue = adjust(fseas.pvalue, BenjaminiHochberg())
        sort!(fseas, :qvalue)
        CSV.write("data/outputs/$(eeg_feat)_noage_gsea.csv", fseas)
        fseas
    end
end


#-

using CairoMakie

#-

function plot_fsea!(grid, setcors, notcors; label="", ylabel="enrichment", linecolor=:gray, kwargs...)
    fullcors = [setcors; notcors]
    ncors = length(fullcors)

    srt = sortperm(fullcors)
    ranks = invperm(srt)
    setranks = Set(ranks[1:length(setcors)])
    
    setscore =  1 / length(setcors)
    notscore = -1 / length(notcors)
    
    xs = 1:ncors
    ys = cumsum(i ∈ setranks ? setscore : notscore for i in eachindex(ranks)) .* -1
    
    t = "ES: $(round(enrichment_score(setcors, notcors), digits=3))"
    !isempty(label) && (t = string(label, " ", t))

    ax1 = Axis(grid[1,1]; title=t, ylabel, kwargs...)
    hidexdecorations!(ax1)

    ax2 = Axis(grid[2,1])
    hidedecorations!(ax2)
    
    lines!(ax1, xs, ys; color=linecolor)
    vlines!(ax2, ranks[1:length(setcors)]; color=:black)
    
    rowsize!(grid, 2, Relative(1/8))
    rowgap!(grid, 1, Fixed(0))

    linkxaxes!(ax1, ax2)
    tightlimits!.((ax1, ax2))
    grid, ax1, ax2
end

#-

amp_N2 = CSV.read("data/outputs/peak_amp_N2_lms.csv", DataFrame)
lat_N1 = CSV.read("data/outputs/peak_latency_N1_lms.csv", DataFrame)

@assert names(gfs_wide, r"UniRef") == amp_N2.feature
@assert names(gfs_wide, r"UniRef") == lat_N1.feature

#-

fig = Figure()
grid1=GridLayout(fig[1,1])
grid2=GridLayout(fig[2,1])
grid3=GridLayout(fig[1,2])
grid4=GridLayout(fig[2,2])

plot_fsea!(grid1, amp_N2[nact["Tryptophan synthesis"], "z"], amp_N2[Not(nact["Tryptophan synthesis"]), "z"]; label="Tryptophan synthesis - amp N2")
plot_fsea!(grid2, lat_N1[nact["DOPAC synthesis"], "z"], lat_N1[Not(nact["DOPAC synthesis"]), "z"]; label="DOPAC synthesis - lat N1")
plot_fsea!(grid3, lat_N1[nact["Isovaleric acid synthesis"], "z"], lat_N1[Not(nact["Isovaleric acid synthesis"]), "z"]; label="Isovaleric acid synthesis - lat N1")
plot_fsea!(grid4, lat_N1[nact["S-Adenosylmethionine synthesis"], "z"], lat_N1[Not(nact["S-Adenosylmethionine synthesis"]), "z"]; label="SAM synthesis - lat N1")

fig

#-

aoff_cors = ThreadsX.map(enumerate(gfs.feature)) do (i, feature)
    abunds = collect(gfs[i, 2:end])
    cor(eegbase.visual_Average_aper_offset, abunds)
end

offset_fseas = DataFrame(ThreadsX.map(collect(keys(nact))) do gs
    ixs = nact[gs]
    isempty(ixs) && return (; cortest = "visual_Average_aper_offset", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

    cs = aoff_cors[ixs]
    isempty(cs) && return (; cortest = "visual_Average_aper_offset", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

    acs = aoff_cors[Not(ixs)]
    mwu = MannWhitneyUTest(cs, acs)
    es = enrichment_score(cs, acs)

    return (; cortest = "visual_Average_aper_offset", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
end)

subset!(offset_fseas, :pvalue=> ByRow(!isnan))
offset_fseas.qvalue = adjust(offset_fseas.pvalue, BenjaminiHochberg())
sort!(offset_fseas, :qvalue)

#- 

fig = Figure()
grid1=GridLayout(fig[1,1])

plot_fsea!(grid1, aoff_cors[nact["p-Cresol synthesis"]], aoff_cors[Not(nact["p-Cresol synthesis"])]; label="p-Cresol synthesis")

fig

#-

#-

gfs_vep = outerjoin(map(eegvep.file) do f
    df = CSV.read(f, DataFrame; skipto = 2, header=["feature", replace(basename(f), r"_S\d+_genefamilies\.tsv"=> "")])
    subset!(df, "feature"=> ByRow(f-> !contains(f, '|')))
end...; on="feature")

foreach(n-> gfs_vep[!,n] = coalesce.(gfs_vep[!,n], 0.), names(gfs_vep))

#-

fig = Figure(;resolution = (600, 1800))

vep_fseas = DataFrame()
nact = getneuroactive(replace.(gfs_vep.feature, "UniRef90_"=>""))

i = 0

for eegfeat in ("peak_amp_N1", "peak_latency_N1", "peak_amp_P1", "peak_latency_P1", "peak_amp_N2", "peak_latency_N2")
    @info eegfeat
    vep_cors = ThreadsX.map(enumerate(gfs_vep.feature)) do (i, feature)
        abunds = collect(gfs_vep[i, 2:end])
        cor(eegvep[!, eegfeat], abunds)
    end

    df = DataFrame(map(collect(keys(nact))) do gs
        ixs = nact[gs]
        isempty(ixs) && return (; cortest = eegfeat, geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = vep_cors[ixs]
        isempty(cs) && return (; cortest = eegfeat, geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = vep_cors[Not(ixs)]
        mwu = MannWhitneyUTest(cs, acs)
        es = enrichment_score(cs, acs)

        if any([(eegfeat == "peak_amp_N1"      && gs == "Menaquinone synthesis"),
                 (eegfeat == "peak_latency_P1"  && gs == "GABA synthesis"),
                 (eegfeat == "peak_latency_P1"  && gs == "Glutamate degradation"),
                 (eegfeat == "peak_amp_N2"      && gs == "Quinolinic acid degradation"),
                 (eegfeat == "peak_amp_N2"      && gs == "Butyrate synthesis")])
            global i += 1
            
            grid=GridLayout(fig[i,1])
            plot_fsea!(grid, vep_cors[nact[gs]], aoff_cors[Not(nact[gs])]; label=gs, ylabel="$eegfeat enrichment")
        end


        return (; cortest = eegfeat, geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)
    append!(vep_fseas, df)
end

subset!(vep_fseas, :pvalue=> ByRow(!isnan))
transform!(groupby(vep_fseas, :cortest), "pvalue" => (p-> adjust(collect(p), BenjaminiHochberg())) => "qvalue")
sort!(vep_fseas, :qvalue)

fig

#-

indf = gfs_wide
outfile = "data/outputs/$(eeg_feat)_noage_lms.csv"
respcol = eeg_feat
featurecols = names(gfs_wide, r"UniRef")
formula = term(:func) ~ term(respcol)

#-

newidx = let newseqs = Set(readlines("new_samples.txt"))
    idx = findall(s-> s ∈ newseqs, gfs_wide.seqprep)
end

describe(gfs_wide[newidx, "age_weeks"])
describe(gfs_wide[Not(newidx), "age_weeks"])


#- 
