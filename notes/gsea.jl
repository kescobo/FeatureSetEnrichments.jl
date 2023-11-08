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
age_model = CSV.read("data/khula_age_residuals.csv", DataFrame)
# If new environment, requires running `VKCComputing.set_readonly_pat!()` and `VKCComputing.set_airtable_dir!()`
base = LocalBase()

# `Number_Segs_Post-Seg_Rej` = number of retained trials

## corrections

# N1 latency and amplitude stays the same
# then corrected P1 latency is calculated as: P1 latency - N1 latency
# and corrected N2 latency is : N2 latency - P1 latency
# corrected P1 amplitude is: P1 amplitude - N1 amplitude
# and corrected N2 amplitude is: N2 amplitude - P1 amplitude

eegvep.peak_latency_P1_corrected = eegvep.peak_latency_P1 .- eegvep.peak_latency_N1
eegvep.peak_latency_N2_corrected = eegvep.peak_latency_N2 .- eegvep.peak_latency_P1
eegvep.peak_amp_P1_corrected = eegvep.peak_amp_P1 .- eegvep.peak_amp_N1
eegvep.peak_amp_N2_corrected = eegvep.peak_amp_N2 .- eegvep.peak_amp_P1

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
                         "peak_latency_N2",
                         "peak_latency_P1_corrected",
                         "peak_latency_N2_corrected",
                         "peak_amp_P1_corrected",
                         "peak_amp_N2_corrected",
                );
    on="subject"
)

select!(gfs_wide, "subject", "seqprep","sex","age_weeks", "trials",
                  "peak_amp_N1","peak_latency_N1",
                  "peak_amp_P1","peak_latency_P1",
                  "peak_amp_N2","peak_latency_N2",
                  "peak_latency_P1_corrected",
                  "peak_latency_N2_corrected",
                  "peak_amp_P1_corrected",
                  "peak_amp_N2_corrected",
                  Cols(Not("UNMAPPED"))
)

#-

subset!(gfs_wide, "subject"=> ByRow(!=("191-34303377"))) # premie
leftjoin!(gfs_wide, select(age_model, "sample"=> "seqprep", "test_prediction"=> "predicted_age"); on="seqprep")

select!(gfs_wide, Not(r"UniRef90"), sort(names(gfs_wide, r"UniRef")))


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
                age_col = "age_weeks",
                formula = term(:func) ~ term(respcol) + term(age_col) + term(:trials),
        )
        @debug "Respcol: $respcol"
    lmresults = DataFrame(ThreadsX.map(featurecols) do feature
        # @debug "Feature: $feature"

        # ab = collect(indf[!, feature] .+ (minimum(indf[over0, feature])) / 2) # add half-minimum non-zerovalue
        df = select(indf, respcol, age_col, "trials")
        over0 = indf[!, feature] .> 0
        
        df.func = over0
        # @debug "DataFrame: $df"

        try
            mod = glm(formula, df, Binomial(), LogitLink(); dropcollinear=false)
            ct = DataFrame(coeftable(mod))
            ct.feature .= feature
            rename!(ct, "Pr(>|z|)"=>"pvalue", "Lower 95%"=> "lower_95", "Upper 95%"=> "upper_95", "Coef."=> "coef", "Std. Error"=>"std_err")
            select!(ct, Cols(:feature, :Name, :))
            return NamedTuple(only(filter(row-> row.Name == respcol, eachrow(ct))))    
        catch e
            @warn "hit $e for $feature"
            return (; feature, Name = respcol, coef = NaN, std_err = NaN, z = NaN, pvalue = NaN, lower_95 = NaN, upper_95 = NaN, qvalue = NaN)
        end
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
    mapreduce(vcat, ["peak_amp_N1","peak_latency_N1","peak_amp_P1","peak_latency_P1","peak_amp_N2","peak_latency_N2",
                     "peak_latency_P1_corrected", "peak_latency_N2_corrected", "peak_amp_P1_corrected", "peak_amp_N2_corrected"]) do eeg_feat 
        
        @info "EEG feat: $eeg_feat"
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
    mapreduce(vcat, ["peak_amp_N1","peak_latency_N1","peak_amp_P1","peak_latency_P1","peak_amp_N2","peak_latency_N2",
                     "peak_latency_P1_corrected", "peak_latency_N2_corrected", "peak_amp_P1_corrected", "peak_amp_N2_corrected"]) do eeg_feat 

        @info "EEG feat: $eeg_feat"
        res = runlms(gfs_wide, "data/outputs/$(eeg_feat)_noage_lms.csv", eeg_feat, names(gfs_wide, r"UniRef"); formula = term(:func) ~ term(eeg_feat) + term(:trials))

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

res_predage = let nact = getneuroactive(replace.(names(gfs_wide, r"UniRef"), "UniRef90_"=>""))
    mapreduce(vcat, ["peak_amp_N1","peak_latency_N1","peak_amp_P1","peak_latency_P1","peak_amp_N2","peak_latency_N2",
                     "peak_latency_P1_corrected", "peak_latency_N2_corrected", "peak_amp_P1_corrected", "peak_amp_N2_corrected"]) do eeg_feat 

        @info "EEG feat: $eeg_feat"
        res = runlms(gfs_wide, "data/outputs/$(eeg_feat)_pred_lms.csv", eeg_feat, names(gfs_wide, r"UniRef"); 
                    age_col = "predicted_age",
                    formula = term(:func) ~ term(eeg_feat) + term(:predicted_age) + term(:trials))
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
        CSV.write("data/outputs/$(eeg_feat)_pred_gsea.csv", fseas)
        fseas
    end
end

#-

res.model .= "base"
res_noage.model .= "no age"
res_predage.model .= "pred. age"

res_comb = let keep = filter(gs-> length(nact[gs]) > 5, unique(res.geneset))
    comb = vcat(res, res_noage, res_predage)
    subset!(comb, "geneset"=> ByRow(g-> g ∈ keep))

    transform!(groupby(comb, ["model", "cortest"]), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg())) => "qvalue")

    subset!(comb, "qvalue"=> ByRow(<(0.2)))
    sort!(comb, ["model", "qvalue"])
end

#-

res_full = let nact = getneuroactive(replace.(names(gfs_wide, r"UniRef"), "UniRef90_"=>""); consolidate = false)
    mapreduce(vcat, ["peak_amp_N1","peak_latency_N1","peak_amp_P1","peak_latency_P1","peak_amp_N2","peak_latency_N2",
                     "peak_latency_P1_corrected", "peak_latency_N2_corrected", "peak_amp_P1_corrected", "peak_amp_N2_corrected"]) do eeg_feat 
        
        @info "EEG feat: $eeg_feat"
        res = runlms(gfs_wide, "data/outputs/$(eeg_feat)_lms_full.csv", eeg_feat, names(gfs_wide, r"UniRef"))

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
        CSV.write("data/outputs/$(eeg_feat)_gsea_full.csv", fseas)
        fseas
    end
end

res_noage_full = let nact = getneuroactive(replace.(names(gfs_wide, r"UniRef"), "UniRef90_"=>""); consolidate = false)
    mapreduce(vcat, ["peak_amp_N1","peak_latency_N1","peak_amp_P1","peak_latency_P1","peak_amp_N2","peak_latency_N2",
                     "peak_latency_P1_corrected", "peak_latency_N2_corrected", "peak_amp_P1_corrected", "peak_amp_N2_corrected"]) do eeg_feat 

        @info "EEG feat: $eeg_feat"
        res = runlms(gfs_wide, "data/outputs/$(eeg_feat)_noage_lms_full.csv", eeg_feat, names(gfs_wide, r"UniRef"); formula = term(:func) ~ term(eeg_feat) + term(:trials))

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
        CSV.write("data/outputs/$(eeg_feat)_noage_gsea_full.csv", fseas)
        fseas
    end
end

res_predage_full = let nact = getneuroactive(replace.(names(gfs_wide, r"UniRef"), "UniRef90_"=>""); consolidate = false)
    mapreduce(vcat, ["peak_amp_N1","peak_latency_N1","peak_amp_P1","peak_latency_P1","peak_amp_N2","peak_latency_N2",
                     "peak_latency_P1_corrected", "peak_latency_N2_corrected", "peak_amp_P1_corrected", "peak_amp_N2_corrected"]) do eeg_feat 

        @info "EEG feat: $eeg_feat"
        res = runlms(gfs_wide, "data/outputs/$(eeg_feat)_pred_lms_full.csv", eeg_feat, names(gfs_wide, r"UniRef"); 
                    age_col = "predicted_age",
                    formula = term(:func) ~ term(eeg_feat) + term(:predicted_age) + term(:trials))
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
        CSV.write("data/outputs/$(eeg_feat)_pred_gsea_full.csv", fseas)
        fseas
    end
end

#-

res_full.model .= "base"
res_noage_full.model .= "no age"
res_predage_full.model .= "pred. age"

nact_full = getneuroactive(replace.(names(gfs_wide, r"UniRef"), "UniRef90_"=>""); consolidate = false)

res_comb_full = let keep = filter(gs-> length(nact_full[gs]) > 5, unique(res_full.geneset))
    comb = vcat(res_full, res_noage_full, res_predage_full)
    subset!(comb, "geneset"=> ByRow(g-> g ∈ keep))

    transform!(groupby(comb, ["model", "cortest"]), "pvalue"=> (p-> adjust(collect(p), BenjaminiHochberg())) => "qvalue")

    subset!(comb, "qvalue"=> ByRow(<(0.2)))
    sort!(comb, ["model", "qvalue"])
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

amp_N2 = sort(CSV.read("data/outputs/peak_amp_N2_lms.csv", DataFrame), "feature")
lat_N1 = sort(CSV.read("data/outputs/peak_latency_N1_lms.csv", DataFrame), "feature")
lat_N1_noage = sort(CSV.read("data/outputs/peak_latency_N1_noage_lms.csv", DataFrame), "feature")
amp_N2_noage = sort(CSV.read("data/outputs/peak_amp_N2_noage_lms.csv", DataFrame), "feature")
amp_P1_noage = sort(CSV.read("data/outputs/peak_amp_P1_noage_lms.csv", DataFrame), "feature")

# amp_N2_pred = sort(CSV.read("data/outputs/peak_amp_N2_pred_lms.csv", DataFrame), "feature")
amp_P1_pred = sort(CSV.read("data/outputs/peak_amp_P1_pred_lms.csv", DataFrame), "feature")
amp_P1cor_pred = sort(CSV.read("data/outputs/peak_amp_P1_corrected_pred_lms.csv", DataFrame), "feature")
lat_N2cor_pred = sort(CSV.read("data/outputs/peak_latency_N2_corrected_pred_lms.csv", DataFrame), "feature")

#- 

@assert names(gfs_wide, r"UniRef") == amp_N2.feature
@assert names(gfs_wide, r"UniRef") == lat_N1.feature
@assert names(gfs_wide, r"UniRef") == lat_N1_noage.feature
@assert names(gfs_wide, r"UniRef") == amp_N2_noage.feature
@assert names(gfs_wide, r"UniRef") == amp_P1_noage.feature

# @assert names(gfs_wide, r"UniRef") == amp_N2_pred.feature
@assert names(gfs_wide, r"UniRef") == amp_P1_pred.feature
@assert names(gfs_wide, r"UniRef") == amp_P1cor_pred.feature
@assert names(gfs_wide, r"UniRef") == lat_N2cor_pred.feature

#-

fig = Figure();
grid1=GridLayout(fig[1,1])
grid2=GridLayout(fig[2,1])
grid3=GridLayout(fig[1,2])
grid4=GridLayout(fig[2,2])

plot_fsea!(grid1, amp_N2[nact["Tryptophan synthesis"], "z"], amp_N2[Not(nact["Tryptophan synthesis"]), "z"]; label="Tryptophan synthesis - amp N2")
plot_fsea!(grid2, lat_N1[nact["DOPAC synthesis"], "z"], lat_N1[Not(nact["DOPAC synthesis"]), "z"]; label="DOPAC synthesis - lat N1")
plot_fsea!(grid3, lat_N1[nact["Isovaleric acid synthesis"], "z"], lat_N1[Not(nact["Isovaleric acid synthesis"]), "z"]; label="Isovaleric acid synthesis - lat N1")
plot_fsea!(grid4, lat_N1[nact["Menaquinone synthesis"], "z"], lat_N1[Not(nact["Menaquinone synthesis"]), "z"]; label="Menaquinone synthesis - lat N1")

Makie.Label(fig[0,1:2], L"bug \sim eeg + age + trials"; fontsize=30)

save("data/figures/sig_fsea.svg", fig)
fig


#-

#-

fig = Figure(; resolution=(600, 900))
grid1=GridLayout(fig[1,1])
grid2=GridLayout(fig[2,1])


plot_fsea!(grid1, amp_N2_noage[nact["Tryptophan synthesis"], "z"], amp_N2_noage[Not(nact["Tryptophan synthesis"]), "z"]; label="Tryptophan synthesis - amp N2")
plot_fsea!(grid2, amp_P1_noage[nact["Glutamate synthesis"], "z"], amp_P1_noage[Not(nact["Glutamate synthesis"]), "z"]; label="Glutamate synthesis - lat N1")

Makie.Label(fig[0,1], L"bug \sim eeg + trials"; fontsize=30, tellwidth=false)

save("data/figures/sig_fsea_no_age.svg", fig)
fig



#- 

genesets = unique(res_comb.geneset)

gs_gfs = mapreduce(vcat, genesets) do geneset
    urefs = Set(names(gfs_wide, r"UniRef")[nact[geneset]])

    gs_gfs = mapreduce(vcat, subset(eegvep, "file"=> ByRow(!ismissing)).path) do f
        seqprep = replace(basename(f), r"_S\d+_genefamilies\.tsv"=> "")
        df = CSV.read(f, DataFrame; skipto = 2, header=["feature", "abundance"])
        df.seqprep .= seqprep
        df.subject .= prep2sub[seqprep]
        df.geneset .= geneset
        subset!(df, "feature"=> ByRow(f-> first(split(f, '|')) ∈ urefs))
    end

    transform!(gs_gfs, "feature" => ByRow(f-> begin
        spl = String.(split(f, '|'))
        uniref = first(spl)
        length(spl) == 1 && return (; uniref, genus = missing, species = missing)
        spl[2] == "unclassified" && return (; uniref, genus = "unclassified", species = "unclassified")
        (genus, species) = String.(split(spl[2], '.'))
        return (; uniref, genus, species)
    end) => ["uniref", "genus", "species"])

    subset!(gs_gfs, "species"=> ByRow(!ismissing))
end


gf_speccontrib = DataFrames.combine(groupby(gs_gfs, ["geneset", "genus"]), "abundance"=> sum => "abundance_total")
transform!(groupby(gf_contrib, ["geneset"]), "abundance_total"=> (a-> a ./ sum(a)) => "abundance_perc")

#-

top12 = subset(gf_contrib, "abundance_perc" => ByRow(>(0.05))).genus |> unique

top12df = DataFrames.combine(groupby(gf_contrib, "geneset"), AsTable(["genus", "abundance_perc"]) => (nt -> begin
    keep = findall(g-> g ∈ top12, nt.genus)
    other = sum(nt.abundance_perc[Not(keep)])
    
    return (; genus = [nt.genus[keep]; "other"], abundance_perc = [nt.abundance_perc[keep]; other])
end) => ["genus", "abundance_perc"])


#-
using ColorSchemes
using ColorSchemes.ColorTypes

let genera = filter(!=("other"), unique(top12df.genus))
    genus_groups = Dict(g=>  i for (i, g) in enumerate(genera))
    genus_groups["other"] = 13
    transform!(top12df, "genus"=> ByRow(g-> genus_groups[g])=> "genus_group")
    transform!(top12df, "genus_group"=> ByRow(gg-> gg == 13 ? RGB(0.85,0.85,0.85) : ColorSchemes.Paired_12[gg])=> "color")
end

let genesets = unique(top12df.geneset)
    gs_groups = Dict(g=> i for (i, g) in enumerate(genesets))
    transform!(top12df, "geneset"=> ByRow(g-> gs_groups[g])=> "geneset_group")
end

#- 

fig, ax, bp = barplot(top12df.geneset_group, top12df.abundance_perc;
    stack = top12df.genus_group,
    color = top12df.color,
    colormap = :Paired_12,
    axis = (; xticks = (1:maximum(top12df.geneset_group), unique(top12df.geneset)), 
              xticklabelrotation = π / 4,
              title = "Genus contribution to genesets")
)

Legend(fig[1,2],
    [MarkerElement(; color = col, marker = :rect) for col in unique(top12df.color)],
    unique(top12df.genus))

fig
#- 

#-

gf_speccontrib = DataFrames.combine(groupby(gs_gfs, ["geneset", "species"]), "abundance"=> sum => "abundance_total")
transform!(groupby(gf_speccontrib, ["geneset"]), "abundance_total"=> (a-> a ./ sum(a)) => "abundance_perc")

topspec = subset(gf_speccontrib, "abundance_perc" => ByRow(>(0.05))).species |> unique

topspecdf = DataFrames.combine(groupby(gf_speccontrib, "geneset"), AsTable(["species", "abundance_perc"]) => (nt -> begin
    keep = findall(g-> g ∈ topspec, nt.species)
    other = sum(nt.abundance_perc[Not(keep)])
    
    return (; species = [nt.species[keep]; "other"], abundance_perc = [nt.abundance_perc[keep]; other])
end) => ["species", "abundance_perc"])


#-
using ColorSchemes
using ColorSchemes.ColorTypes

let genera = filter(!=("other"), unique(topspecdf.species))
    species_groups = Dict(g=>  i for (i, g) in enumerate(genera))
    species_groups["other"] = 13
    transform!(topspecdf, "species"=> ByRow(g-> species_groups[g])=> "species_group")
    transform!(topspecdf, "species_group"=> ByRow(gg-> gg == 13 ? RGB(0.85,0.85,0.85) : ColorSchemes.Paired_12[gg])=> "color")
end

let genesets = unique(topspecdf.geneset)
    gs_groups = Dict(g=> i for (i, g) in enumerate(genesets))
    transform!(topspecdf, "geneset"=> ByRow(g-> gs_groups[g])=> "geneset_group")
end

#- 

fig, ax, bp = barplot(topspecdf.geneset_group, topspecdf.abundance_perc;
    stack = topspecdf.species_group,
    color = topspecdf.color,
    colormap = :Paired_12,
    axis = (; xticks = (1:maximum(topspecdf.geneset_group), unique(topspecdf.geneset)), 
              xticklabelrotation = π / 4,
              title = "species contribution to genesets")
)

Legend(fig[1,2],
    [MarkerElement(; color = col, marker = :rect) for col in unique(topspecdf.color)],
    unique(topspecdf.species))

fig
#- 
