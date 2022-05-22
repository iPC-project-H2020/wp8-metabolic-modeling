# D8.3 Metabolic models
This repository provide the codes of the pipeline developped in WP8.3 for creating context-specific metabolic models from scRNA-seq data.
To run the codes you will need to have intalled in your system R, Matlab and the Cobra toolbox. 

The pipeline consists of three main steps. A first pre-processing step where dimensionality reduction, clustering and cluster annotation is performed on the scRNA-seq and a unique gene expression profile per cluster is retrieved. A second step where gene expression is integrated into a human genome-scale model by means of the GIMME algorithm. The third step consists in tuning the parameters from the previous steps with the rationale that clusters with the same annotation across multiple datasets should show a similar distribution of metabolic fluxes.

## Code organization
The repository is organised as follows:
1. Rcodes and MATLABcodes folders contain the R and MATLAB codes for the construction of context-specific metabolic models from scRNA-seq data.
2. scRNA-seq_data folder should contain the scRNA-seq datasets. Here we provide the PDX scRNA-seq data from Aynaud et al. 2019. 
3. GEM folder contains the initial GEM models. Here, we provide the GEMs generated in Masid et al. 2020. 
4. IC_genes folder contains the gene set signatures used for clustering annotation. 

To run the first two steps of the pipeline, one would need to setup the 'job_profile.xlsx' file in the parent folder.
This file must include a:
1. ‘file.name’ column with the scRNA-seq file names.
2. ‘is.cell_line’ column with logicals indicating if the scRNA-seq data is from induced cell lines or not;
3. ‘model_name’ column with the initial GEM file name; 
4. ‘ic_score’ column with the method to calculate gene signature scores (e.g. median);
5. ‘aggregation_by’ column with the aggregation method to construct a unique gene expression profile per cluster; 
6. ‘clustering_resolution column with the clustering resolution parameter value.  
7. 'input_type_clustering' (only needed for inducible cell-lines) which is set to 'temporal' or to 'ic_signature' to select to construct context-specific models by time or a specific gene signature called 'ic.xlsx' in the IC_genes folder. 

The algorithms will generate 5 folders organized in sub-folders by cell-line or pdx, the initial GEM model name, clustering resolution, gene signature scoring method and clustering aggregation method. The folders that are generated are:
1. ‘cluster_profiles’ which contains the gene expression profiles. 
2. ‘cluster_info’ which contains the gene signature scores. 
3. ‘ValToReact’ which contains the mapped cluster gene expression values to reactions. 
4. ‘FBA_solution’ which contains the FBA solutions. 
5. ‘csGEM’ which contains the context-specific GEMs.

## Running the first two steps of the pipeline
To run the first two steps of the pipeline you will need to:
1. clone the folder and unzip it;
2. set the working directory in the Rcodes/work_flow.R to the unzipped folder;
3. set the job_profiles.xlsx file;
4. upload the scRNA-seq data in the scRNAseq_data folder;
5. upload to the IC_score folder the gene signatures you would like to use;  
6. locate the directory of the Rscript.exe file on your system
7. Open the command prompt
- go to the unzipped folder
- run: '/dir_to_Rscript/Rscript.exe' Rcodes/work_flow.R

## Application to scRNA-seq data from Ewing Sarcoma
we considered the dataset from (Aynaud et al. 2020) which comprised temporal scRNA-seq data from A63 cell line upon EWSR1-FLI1 silencing, scRNA-seq from patient-derived xenografts (PDX) under DOX+/DOX- treatment and 8 patient-derived xenografts. In addition, we consider an independent corpus of 7 unpublished scRNA-seq data from patient-derived xenograft. 
Temporal scRNA-seq data from A63 cell lines and patient-derived xenografts under DOX+/DOX- treatment were exploited to study the effect of EWSR1-FLI1 on the metabolism by clustering cells based on gene signature associated to EwS activity.

The 'job_profiles.xlsx' here provided includes the set of jobs we have runned with the pipeline. In folder scRNAseq_data, IC_genes and GEM, you can find the scRNA-seq data for the 7 PDXs from Aynaud et al. 2019, the genes signatures associated to EWSR1-FLI1 activity (ic10), cell cycle (G1/S and G2/M phases, ic1 and ic2), oxidative phosphorylation (ic4) and glucose catabolism (ic14) they identified in their study and the metabolic models used in this work.  

We provide in folders cluster_info, cluster_profiles, csGEM, FBA_solution, ValToReact results for which the silhouette scores where the maximum. 
File XX contains the results obtained using as initial model Recon3_redHuman, a thermodynamically curated version of Recon3. 
