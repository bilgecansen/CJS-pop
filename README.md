# CJS-pop

This document explains the R scripts used in Şen and Akçakaya (2019), which describes CJS-pop, applies it to simulated and MAPS data, and uses the estimated parameters for population projections. 

The code in general are presented as "notebooks" which are R Markdown files. You can knit the whole file or you can run the code chunks seperately. Each notebook includes detailed notes about the code for data wrangling and analysis.

We can't share MAPS data. Brow Creeper analysis code will not work unless you have access to full MAPS data. However, we share an example analysis with simulated data under **sim_example_CJSpop_run**. Clone the repository and run **sim_models_notebook.rmd** generate the .jags files. Then run **sim_example_CJSpop_run** to obtain example results.

If you have access to data the main analysis has three parts: MAPS, Simulations, and Projections and Plotting.

## 1. MAPS

Run the **brcr_data_wrangling_notebook** to prepare Brown Creeper MAPS data for CJS-pop in JAGS. It creates a new directory named **data** if it is already not created, and stores capture history data.

Run **brcr_models_notebook** to generate the .jags files used in applying CJS-pop to BRCR capture history data. Easiest way is to knit the whole notebook.

Run **brcr_models.R** to apply 4 different mark-recapture models (3 different IRMs and one CJS) to Brown Creeper capture history data. Results from this analysis are stored in the **results** sub-directory.

## 2. Simulations

Run **data_generation_notebook** to simulate capture history data as explained in Şen and Akçakaya (2019). It creates a new directory named **data** if it is already not created, that stores all of the generated data.

Run **sim_models_notebook** to generate the .jags files used in applying CJS-pop to simulated data. Easiest way is to knit the whole notebook.

We analysied the simulated data sets using the high performance comuputing cluster of Stony Brook University, which is named SeaWulf. We provide a sample code, **sim_example_IRM_run.R**, that can be run locally and analyses a single simulation set. This sample code can easily be extended to analyse all simulated sets in a single script in another cluster. The example code saves results locally under **results** sub-directory.

## 3. Projections and Plotting

The plots and population projections in Şen and Akçakaya et al. (2019) can be reproduced with **Plots.R**.


