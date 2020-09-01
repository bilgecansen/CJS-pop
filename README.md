# CJS-pop

This document explains the R scripts used in Şen and Akçakaya (2020; available as preprint also under review in Ecology), which describes CJS-pop, applies it to simulated and MAPS data, and uses the estimated parameters for population projections. 

The code in general are presented as "notebooks" which are R Markdown files. You can knit the whole file or you can run the code chunks seperately. Each notebook includes detailed notes about the code for data wrangling and analysis.

We don't have permission to share MAPS data directly. Brown Creeper analysis code will not work unless you have access to full MAPS data. However, we share an example analysis with simulated data under **sim_example_CJSpop_run.R**. Clone the repository and run **sim_models_notebook.rmd** to generate the .jags files. Then run **sim_example_CJSpop_run.R** to obtain example results. We'll update this reposttiory if we get permisson to share Brown Creeper data. Meanwhile the results of the Brown Creeper models can be accessed at **BRCR_results.rds**. 

If you have access to data, the main analysis has three parts: MAPS, Simulations, and Projections and Plotting.

## 1. MAPS

Run the **brcr_data_wrangling_notebook.rmd** to prepare Brown Creeper MAPS data for CJS-pop in JAGS.

Run **brcr_models_notebook.rmd** to generate the .jags files used in applying CJS-pop to BRCR capture history data. Easiest way is to knit the whole notebook.

Run **brcr_models.R** to apply 4 different mark-recapture models (3 different CJS-pops and one CJS) to Brown Creeper capture history data. These results are also available in the file **BRCR_results.rds**.

## 2. Simulations

Run **data_generation_notebook.rmd** to simulate capture history data as explained in Şen and Akçakaya (in review). 

Run **sim_models_notebook.rmd** to generate the .jags files used in applying CJS-pop to simulated data. Easiest way is to knit the whole notebook.

We analysied the simulated data sets using the high performance comuputing cluster of Stony Brook University, which is named SeaWulf. We provide a sample code, **sim_example_CJSpop_run.R**, that can be run locally and analyses a single simulation set. This sample code can easily be extended to analyse all simulated sets in a single script in another cluster. 

## 3. Projections and Plotting

The plots and population projections in Şen and Akçakaya (in review) can be reproduced with **Plots.R**.

## Preprint
Şen, B., and H. R. Akçakaya. 2020. Estimating fecundity and density dependence from mark-recapture data for making population projections. bioRxiv:2020.08.26.268656. https://www.biorxiv.org/node/1501472.external-links.html



