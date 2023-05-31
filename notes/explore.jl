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

#-

eegbase."subject_ID"

#-

subjects = let all3 = uids(base["Subjects"]) ∩ eegbase."subject_ID" ∩ metadata."subject_id"
    [base["Subjects", sub] for sub in all3]
end

seqs = DataFrame(mapreduce(vcat, base[mapreduce(r-> [b.id for b in base[r[:Biospecimens]]], vcat, subjects)]) do rec
    sample_id = rec[:uid]
    subject_id = base[first(rec[:subject])][:uid]
    seqpreps = base[rec[:seqprep]]

    [(; sample_id, subject_id, seq_id = seq[:uid]) for seq in seqpreps]
end)

leftjoin!(seqs, metadata; on="subject_id")
genefamilies_files = filter(readdir("/grace/sequencing/processed/mgx/humann/main", join=true)) do path
    file = basename(path)
    m = match(r"(SEQ\d+)_S\d+_genefamilies.tsv", file)
    isnothing(m) && return false
    return m[1] in seqs.seq_id
end

leftjoin!(seqs, DataFrame(seq_id = map(f-> match(r"(SEQ\d+)_S\d+", f)[1], genefamilies_files), file = genefamilies_files); on= "seq_id")
subset!(seqs, "file"=> ByRow(!ismissing))
unique!(seqs, "subject_id")
#-

leftjoin!(eegbase, seqs; on="subject_id")
subset!(eegbase, "file"=> ByRow(!ismissing))
leftjoin!(eegvep, seqs; on="subject_id")
subset!(eegvep, "file"=> ByRow(!ismissing))

#-

gfs = outerjoin(map(eegbase.file) do f
    df = CSV.read(f, DataFrame; skipto = 2, header=["feature", replace(basename(f), r"_S\d+_genefamilies\.tsv"=> "")])
    subset!(df, "feature"=> ByRow(f-> !contains(f, '|')))
end...; on="feature")

foreach(n-> gfs[!,n] = coalesce.(gfs[!,n], 0.), names(gfs))

#-

using Statistics
using ThreadsX
include("fsea.jl")
using .NeuroFSEA
using MultipleTesting

aexp_cors = ThreadsX.map(enumerate(gfs.feature)) do (i, feature)
    abunds = collect(gfs[i, 2:end])
    cor(eegbase.visual_Average_aper_exponent, abunds)
end

nact = getneuroactive(replace.(gfs.feature, "UniRef90_"=>""))

fseas = DataFrame(ThreadsX.map(collect(keys(nact))) do gs
    ixs = nact[gs]
    isempty(ixs) && return (; cortest = "visual_Average_aper_exponent", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

    cs = aexp_cors[ixs]
    isempty(cs) && return (; cortest = "visual_Average_aper_exponent", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

    acs = aexp_cors[Not(ixs)]
    mwu = MannWhitneyUTest(cs, acs)
    es = enrichment_score(cs, acs)

    return (; cortest = "visual_Average_aper_exponent", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
end)

subset!(fseas, :pvalue=> ByRow(!isnan))
fseas.qvalue = adjust(fseas.pvalue, BenjaminiHochberg())
sort!(fseas, :qvalue)
CSV.write(outfile, fseas)


#-

using CairoMakie

hist(aoff_cors)

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

fig = Figure()
grid1=GridLayout(fig[1,1])
grid2=GridLayout(fig[2,1])
plot_fsea!(grid1, aexp_cors[nact["Tryptophan synthesis"]], aexp_cors[Not(nact["Tryptophan synthesis"])]; label="Tryptophan synthesis")
plot_fsea!(grid2, aexp_cors[nact["p-Cresol synthesis"]], aexp_cors[Not(nact["p-Cresol synthesis"])]; label="p-Cresol synthesis")

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

