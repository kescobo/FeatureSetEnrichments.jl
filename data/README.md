# Data

- `eeg_3mo_baseline.csv`
- `eeg_3mo_vep.csv`

Both exported from excel file from Cara Bosco / Laurel Gabard-Durnam (email, 5/22/23).

## Notes

Me: 

> Are the features I'm primarily interested in as response variables `age_3m_eeg`, `visual_Average_aper_offset`, and  `visual_Average_aper_exponent`?
> The columns to the right of those largely look like QC - should I be filtering on any of those, or including any of them as controls?

Laurel:

> yup yup, 
>
> for the first set of models: Please filter out anyone with visual_average_r_squared < 0.9000 ?
> That indicates insufficient model fit. 
> you'll need to control for age or look at microbiomexage interactions (up to you),
> and we also often covary for number of retained baseline EEG trials/segments in our models often.
> The primary brain "outcomes" are the ones you listed:  `visual_Average_aper_offset` (the intercept of the EEG power function)
> and the `visual_Average_aper_exponent` (the slope of the EEG power function....my personal most interested feature). 
>
> For the second set of analyses on the VEP (visual evoked potential) response data,
> everyone we sent you should have good data so no need to filter for model fit or anything here.
> You'll control/covary for age (or look at microbiome interactions), and covary for number of retained VEP trials as well.
> Then you can look at  the VEP latencies (N1, P1, N2) as brain outcomes for your models!
