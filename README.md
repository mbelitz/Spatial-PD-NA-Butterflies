# Spatial-PD-NA-Butterflies
Code for Earl et al. "Spatial phylogenetics of butterflies in relation to environmental drivers and angiosperm diversity across North America"


## In the data repository, there are four subdirectories:
1) AllSeedPlants - Biodiverse outputs for all seep plant data used in Mishler et al. 2020
2) Angiosperms - Biodiverse outputs for the angiosperm plant data adpated from Mishler et al. 2020
3) Biodiverse_Inputs - geographic data used for Biodiverse analyses to generate phylodiversity metrics for the butterflies of North America
4) Biodiverse_Outputs - Biodiverse outputs for the butterflies of North America
5) Model_Covariates - Environmental layers used to predict butterfly phylodiversity metrics 

## In the figure_outputs directory, there are the raw .png's of the Figures used in the main text:
1) Figure 2
2) Figure 3
3) Figure 4
4) Figure 5
5) Figure 6

There are also .png's of Figures shown in the supplemental information:
1) bfly_rich_PD_RPD_PE
2) seed_vs_angiosperms

## In the scripts directory, there are two subdirectories:
1) Figures - code to reproduce main text figures (Fig 2 - Fig 6) and two supplemental figures.
    --- note that correlations referred to in the main text are calculated in the Figure 5 script ---
2) Models - code to reproduce the models presented in the paper
    - Code to run the non spatial GLM model selection and summaries are in the generate_nonSpatial_GLMs script
    - Code to run the spatial GLMs are divided for each phylodiversity metric 
        * Note fitting the spatial GLMs are computational expensive -- gausian models took ~ 1 hour to fit and binomial logit link took ~6 hours to fit
