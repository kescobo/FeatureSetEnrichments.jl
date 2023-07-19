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

txdf = let
    txfiles = filter(readdir("/grace/sequencing/processed/mgx/metaphlan/", join=true)) do path
        file = basename(path)
        m = match(r"(SEQ\d+)_S\d+_profile.tsv", file)
        isnothing(m) && return false
        return m[1] in seqs.seqprep
    end

    DataFrame(map(txfiles) do f
        m = match(r"(SEQ\d+)_(S\d+)_profile.tsv", f)
        (seqprep, S_well) = m.captures
        dir = dirname(f)
        file = basename(f)
        (; seqprep, S_well, dir, file, path=f)
    end)
end

leftjoin!(seqs, txdf; on="seqprep")

#-

leftjoin!(eegbase, subset(seqs, "visit"=> ByRow(==("3mo"))); on="subject_id")
subset!(eegbase, "seqprep"=> ByRow(!ismissing))

leftjoin!(eegvep, subset(seqs, "visit"=> ByRow(==("3mo"))); on="subject_id")
subset!(eegvep, "seqprep"=> ByRow(!ismissing))

#-

prep2sub = Dict(p => s for (p,s) in zip(seqs.seqprep, seqs.subject_id))

txs = mapreduce(vcat, subset(eegvep, "file"=> ByRow(!ismissing)).path) do f
    seqprep = replace(basename(f), r"_S\d+_profile\.tsv"=> "")
    df = CSV.read(f, DataFrame; skipto = 5, header=["feature", "taxid", "abundance", "additional"])
    df.seqprep .= seqprep
    df.subject .= prep2sub[seqprep]
    select!(df, Not(["taxid", "additional"]))
end

txs_wide = unstack(txs, "feature", "abundance",)

foreach(n-> txs_wide[!,n] = coalesce.(txs_wide[!,n], 0.) ./ 100, names(txs_wide)[3:end])
leftjoin!(txs_wide, 
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

select!(txs_wide, "subject", "seqprep","sex","age_weeks", "trials",
                  "peak_amp_N1","peak_latency_N1",
                  "peak_amp_P1","peak_latency_P1",
                  "peak_amp_N2","peak_latency_N2", 
                  Cols(Not("UNKNOWN"))
)

#-

using Statistics
using ThreadsX
using MultipleTesting
using GLM
using HypothesisTests

#-

lmresults = let
    eeg_features = [
        "peak_amp_N1","peak_latency_N1",
        "peak_amp_P1","peak_latency_P1",
        "peak_amp_N2","peak_latency_N2"
    ]

    bugs = filter(names(txs_wide, r"s__")) do bug
        prevalence(txs_wide[!, bug]) > 0.1
    end


    res = mapreduce(vcat, eeg_features) do feature
        @info feature
        DataFrame(ThreadsX.map(bugs) do bug
            # @debug "Feature: $feature"

            df = select(txs_wide, feature, "age_weeks", "trials")
                       
            df.bug = asin.(sqrt.(txs_wide[!, bug]))
            # @debug "DataFrame: $df"

            mod = lm(term(:bug) ~ term(feature) + term(:age_weeks) + term(:trials), df)
            ct = DataFrame(coeftable(mod))
            ct.bug .= last(split(bug, "|"))
            rename!(ct, "Name"=> "feature", "Pr(>|t|)"=>"pvalue", "Lower 95%"=> "lower_95", "Upper 95%"=> "upper_95", "Coef."=> "coef", "Std. Error"=>"std_err")
            select!(ct, Cols(:feature, :bug, :))
            return NamedTuple(only(filter(row-> row.feature == feature, eachrow(ct))))
        end)
    end
end

subset!(lmresults, "pvalue"=> ByRow(!isnan))
DataFrames.transform!(groupby(lmresults, "feature"), :pvalue => (col-> MultipleTesting.adjust(collect(col), BenjaminiHochberg())) => :qvalue)
DataFrames.transform!(lmresults, :pvalue => (col-> MultipleTesting.adjust(collect(col), BenjaminiHochberg())) => "qvalue₀")
sort!(lmresults, :qvalue)

lmresults

#- 

using CairoMakie


#-

fig = Figure(;resolution=(600, 1100))

ax1 = Axis(fig[1,1]; xlabel="E. faecalis (asin sqrt)", ylabel = "P1 latency")
scatter!(ax1, asin.(sqrt.(txs_wide[!, r"s__Enterococcus_faecalis"][!, 1])), txs_wide."peak_latency_P1")

ax2 = Axis(fig[2,1]; xlabel="R. lactiformans (asin sqrt)", ylabel = "P1 amp")
scatter!(ax2, asin.(sqrt.(txs_wide[!, r"s__Ruthenibacterium_lactatiformans"][!, 1])), txs_wide."peak_amp_P1")

ax3 = Axis(fig[3,1]; xlabel="W.lenta (asin sqrt)", ylabel = "P1 amp")
scatter!(ax3, asin.(sqrt.(txs_wide[!, r"s__Eggerthella_lenta"][!, 1])), txs_wide."peak_amp_P1")

fig
