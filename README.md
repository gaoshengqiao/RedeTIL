# RedeTIL Features
## Project Introduction
A Python package for analyzing single-cell data using RedeTIL features. This package provides tools to calculate abundance, spatial, and dynamic features of single-cell data, to predict clinical response to 
Immune checkpoint blockade (ICB) therapy for individual patients. 


## Installation

You can install the `rede_til` package using pip:

```bash
pip install git+https://github.com/gaoshengqiao/RedeTIL

Usage
Initialization
First, you need to import the RedeTIL_Features class and initialize an instance with single-cell data. The single-cell data should be provided in AnnData format, where data.obs must contain the 'label' column for cell annotations.

from rede_til import RedeTIL_Features
import scanpy as sc

# Load your single-cell data
adata = sc.read_h5ad('path_to_your_data.h5ad')

# Initialize RedeTIL_Features
redetil = RedeTIL_Features(adata, target='PDCD1', perturbation='block',
                           T_cells='T-cell', Cancer_cells='malignant',
                           outdir='./results')
Abundance Features
To calculate the abundance features of specific gene+ cell subsets:

redetil.Abundance_features()
This will generate a file Abundance_features.txt in the specified output directory, containing the frequency of each gene+ cell subset in different cell types.
Spatial Features
To reconstruct the spatial topology of cells and calculate spatial features:

redetil.Spatial_features(plot=True)
This will generate several files in the Spatial_features directory, including:
cellinfo_tbl.txt: Cell information table with coordinates and labels.
Cell-Cell_connections.txt: Cell-cell interaction pairs.
Cell-Cell_connections_qvalue.txt: Adjusted q-values for cell-cell interactions.
Detailed_connections.txt: Detailed connections between cells.
Rank Cell-Cell_Distance.txt: Rank distance matrix.
Euclidean Cell-Cell_Distance.txt: Euclidean distance matrix.
Ligand-receptor_contribution.txt: Contribution of ligand-receptor pairs.
If plot=True, it will also generate 3D and 2D plots of the cell spatial distribution and interaction maps.
Dynamic Features
To calculate the dynamic features, which simulate the proximity of T cells to cancer cells after receiving drug perturbation:

redetil.Dynamic_features(plot=True)
This will generate several files in the Dynamic_features directory, including:
Perturbated Spatial_features/: Directory containing perturbated spatial features.
Dynamic_features.txt: Infiltration change of T cells to cancer cells.
If plot=True, it will also generate plots of the perturbated cell spatial distribution and interaction maps.
Example
Here is a complete example to demonstrate the usage of RedeTIL_Features:

if __name__ == '__main__':
    # Load single-cell data
    adata = sc.read_h5ad('./demo/demo.h5ad')

    # One target perturbation
    redetil = RedeTIL_Features(adata, target='PDCD1', perturbation='block',
                               T_cells='T-cell', Cancer_cells='malignant',
                               outdir='./results')
    redetil.Abundance_features()
    redetil.Spatial_features(plot=True)
    redetil.Dynamic_features(plot=True)

    # Targets combo
    redetil = RedeTIL_Features(adata, target='PDCD1', combo_target='VEGFA', perturbation='block',
                               T_cells='T-cell', Cancer_cells='malignant',
                               outdir='./results')
    redetil.Dynamic_features(plot=True)
Contact
For any questions or issues, please contact [your_email@example.com].
License
This project is licensed under the MIT License.

