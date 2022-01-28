# Power analysis for spatial omics

This repo contains code related to _Power analysis for spatial omics_. 

A description of the contents:

+ `codex_spleen/` : 
    - `Spleen_IST_generation.ipynb` : IST generation for mouse spleen, **Figure 2k**
    - `Spleen_NegBinom.ipynb` : Sampling experiments for spleen, **Figure 2l**
    - `ProspectivePower_2mn.ipynb`: Prospective power analysis, **Figure 2m-n**
    - `FigS7_2021-12-03.ipynb` : FOV size experiments and visualization, **Supplementary Figure 7**
    - `Spleen_Cohort_Comparison.ipynb` : Interaction enrichment statistic, **Supplementary Figure 8**
    - `spleen_data/` : Contains support files. 
+ `simulation/` :
    - `FigureS3.ipynb` : Sampling experiments for synthetic data. **Supplementary Figure 3**
    - `FigureS4Heatmaps-2021-12-02.ipynb`: Clustering experiments for ISTs, **Supplementary Figure 4**
+ `osmfish_cortex/` 
    - `osmfish_generation.ipynb` : IST generation for mouse cortex, **Figure 2g**
    - `NB_Cell_Discov_clean.ipynb`: Cell type discovery sampling experiments, **Figure 2h, Supplementary Figure 6d**
    - `data/` : Contains support files
+ `hdst_breastcancer/`
    - `BreastCancer_IST_Generation.ipynb`: IST generation for breast cancer, **Figure 2c**
    - `BreastCancer_NB.ipynb` : Cell type discovery sampling experiments, **Figure 2d**
    - `data/` : Contains support files
+ `spatialpower/` : Package for IST generation and supporting analysis
+ `scripts/` : Contains support scripts for other analyses
    - `generate_tiles.py` : Generates tiles for shuffling analysis corresponding to Figure 2m-n
    - `random_self_pref_cluster.py` :   Generates ISTs for clustograms in Supplementary Figure 4.  
+ `simulated_tissue.ipynb` : An example of fully synthetic tissue generation/optimization. 
+ `sample_results\` : Contains the results and support files for the tissue generation example in `simualted_tissue.ipynb`

## Requirements
We provide a clone of the conda environment used to generate these results in the `env.yml` file. To install the environment, use `conda env create -n <environment name> --file env.yml`.

## Using the `spatialpower` package
To install our software, add the `spatialpower/` directory to your Python path, or run a Jupyter notebook from a directory containing `spatialpower/`. 

## Example usage
We provide example usage of our approach for tissue generation. See the `osmfish_cortex/osmfish_generation.ipynb` file for a complete example with generated tissue from raw osmFISH data. 

Additionally, we provide `simulated_tissue.ipynb` as an example of how to use our tissue generation method on a theoretical tissue (e.g. for testing methods to recover a specific spatial feature). 
