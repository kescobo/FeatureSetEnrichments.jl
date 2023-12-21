module FeatureSetEnrichments

export fsea,
       pvalue,
       enrichment_score,
       neuroactive_map

import HypothesisTests: MannWhitneyUTest
import HypothesisTests: pvalue
import InvertedIndices: Not
import CodecZlib: GzipDecompressorStream
import Dictionaries: Dictionary, set!
import ThreadsX

include("fsea.jl")
include("neurofeatures.jl")

end # module FeatureSetEnrichments