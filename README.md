# Power analysis for spatial omics

This repo contains code related to [_Power analysis for spatial omics_](https://www.biorxiv.org/content/10.1101/2022.01.26.477748v2)

If you use this in your work, please cite: 
**Power analysis for spatial omics.** Ethan Alexander García Baker, Denis Schapiro, Bianca Dumitrascu, Sanja Vickovic, Aviv Regev
*bioRxiv* 2022.01.26.477748; doi: https://doi.org/10.1101/2022.01.26.477748

## Contents

+ `codex_spleen/` : 
    - `Spleen_IST_generation.ipynb` : IST generation for mouse spleen, **Figure 2k**
    - `Spleen_NegBinom.ipynb` : Sampling experiments for spleen, **Figure 2l**
    - `ProspectivePower_2mn.ipynb`: Prospective power analysis, **Figure 2m-n**
    - `FigS7_2021-12-03.ipynb` : FOV size experiments and visualization, **Supplementary Figure 7**
    - `Spleen_Cohort_Comparison.ipynb` : Interaction enrichment statistic, **Supplementary Figure 8**
    - `spleen_binning.ipynb` : Resolution analysis for spleen dataset, **Figure 3**
    - `spleen_multisample.ipynb` : Analysis on the impact of inclusion of multiple samples on power, **Figure 3**
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
    - `HDST_binning.ipynb` : Resolution analysis for HDST breast cancer data, **Figure 3**
    - `BreastCancer_multisample.ipynb` : Analysis on the impact of inclusion of multiple samples on power, **Figure 3**
    - `data/` : Contains support files
+ `spatialpower/` : Package for IST generation and supporting analysis
+ `scripts/` : Contains support scripts for other analyses
    - `generate_tiles.py` : Generates tiles for shuffling analysis corresponding to Figure 2m-n
    - `random_self_pref_cluster.py` :   Generates ISTs for clustograms in Supplementary Figure 4.  

## Requirements
We provide a clone of the conda environment used to generate these results in the env.yml file. To install the environment, use `conda env create -n <environment name> --file env.yml`.

## Generating a tissue _in silico_
We provide a command line Python tool to generate tissue _in silico_. 

Generating an IST requires knowledge of two parameters: a vector describing the abundance of the _k_ cell types, p, and the _k x k_ matrix describing probability that two cell types are directly adjacent, H. These objects are described in the paper. We suggest that _p_ and _H_ are estimated from pilot data; the IST generation notebooks above illustrate how one might do this. 

The generalized steps for the construction of the IST are:
 1. Generate a tissue scaffold
 2. Estimate _p_ and _H_
 3. Label the tissue scaffold. 

#### Generate a tissue scaffold
To construct a tissue scaffold, execute `random_circle_packing.py`:
`python spatialpower/tissue_generation/random_circle_packing.py -x 1000 -y 1000 -o sample_results`

Key arguments that can be adjusted to tune the circle packing are `-x` and `-y`, which control the width and height, respectively, of the rectangular tissue area, and `--rmin` and `--rmax` which control the minimum and maximum radius, respectively, of the circles that are packed within the bounding rectangle to generate the random planar graph. 

A full enumeration of the arguments, including controls for visualization and export, is available using `-h` flag. 

#### Estimate _p_ and _H_
We suggest obtaining values for _p_ and _H_ from a pilot experiment or estimating them by prior knowledge. 

We provide a function for the efficient computation of _H_ given the adjacency matrix, _A_ and one-hot encoded assignment matrix _B_ for a graph representation of a tissue (real or simulated). To calculate _H_:

```python
import spatialpower.neighborhoods.permutationtest as perm_test
H = perm_test.calculate_neighborhood_distribution(A, B)
```

The cell type abundance _p_ can be easily computed using the one-hot encoded assignment matrix _B_:
`p = np.sum(B, axis=0)/np.sum(np.sum(B, axis=0))`

#### Label the tissue scaffold:
To perform a labeling of the tissue scaffold using the optimization approach: 
```python
import spatialpower.tissue_generation.assign_labels as assign_labels
cell_assignments = assign_labels.optimize(A, p, H, learning_rate=1e-5, iterations = 10)
```

To perform a labeling of the tissue scaffold using the heuristic approach: 

```python
import networkx as nx
import spatialpower.tissue_generation.assign_labels as assign_labels

G = nx.from_numpy_array(A)
cell_assignments = assign_labels.heuristic_assignment(G, p, H, mode='graph', dim=1000, position_dict=position_dict, grid_size=50, revision_iters=100, n_swaps=25)
```

### Example notebooks
We provide `simulated_tissue.ipynb` as an example of how to use our tissue generation method on a theoretical tissue (e.g. for testing methods to recover a specific spatial feature) following this approach above.

Additionally, we provide example usage of our approach for tissue generation. See the `osmfish_cortex/osmfish_generation.ipynb` file for a complete example with generated tissue from raw osmFISH data. The run time should be below 30 minutes for the full notebook on a modern laptop.

## Performing power analysis 
Our work provides a general framework for the considerations that should be taken into account in spatial experimental design. In our manuscript, we consider experiments to detect several spatial features, including the discovery of a cell type of interest and the detection of cell-cell interactions. 

We examine experiments to discover these spatial features as illustrative examples of our general framework; we encourage individual users to adapt these approaches for their particular question of interest. 

#### Cell type discovery

In general, the overall procedure for cell type discovery is: 
1. Obtain pilot data
2. Estimate model parameters
3. Calculate probability of detecting cell type of interest given some level of sampling (e.g. number of cells or FOVs sampled)

We provide notebooks implementing this framework for the three differently-structured data sets discussed in our manuscript: 
- `FigureS3.ipynb` : Sampling experiments for synthetic data. **Supplementary Figure 3**
- `Spleen_NegBinom.ipynb` : Sampling experiments for spleen, **Figure 2l**
- `NB_Cell_Discov_clean.ipynb`: Cell type discovery sampling experiments, **Figure 2h, Supplementary Figure 6d**
- `BreastCancer_NB.ipynb` : Cell type discovery sampling experiments, **Figure 2d**

#### Cell-cell interactions
As discussed in the manuscript, we suggest the detection of cell-cell interactions via a permutation test, which we provide code for in `spatialpower/neighborhoods/permutationtest.py`. 

Additionally, we show an illustrative example of this analysis in:
- `FigS7_2021-12-03.ipynb` : FOV size experiments and visualization, **Supplementary Figure 7**
- `ProspectivePower_2mn.ipynb`: Prospective power analysis, **Figure 2m-n**

#### Comparing differently structured tissues
We introduce the interaction enrichment statistic, which is implemented in the functions `calculate_enrichment_statistic` and `z_test` in the module `spatialpower/neighborhoods/permutationtest.py`. 

We provide an illustrative example of the IES test as implemented in our manuscript :
 - `Spleen_Cohort_Comparison.ipynb` : Interaction enrichment statistic, **Supplementary Figure 8**









