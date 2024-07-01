# AF2RAVE_Glide-kinase

## Purpose:
Please note that this is intended to be a code and data repository for a specific set of systems (i.e. protein kinases)
For a view of the full pipeline excluding multiple walker learning and static bias, please use the repository for [AF2RAVE](https://github.com/tiwarylab/alphafold2rave/)

## A simple colab notebook for tAF2 can be found here -> [Link](https://colab.research.google.com/github/tiwarylab/AF2RAVE_Glide-kinase/blob/main/tAlphaFold2.ipynb)
The above mentioned notebook is based on [Colabfold](https://github.com/sokrypton/ColabFold)

## Motivation
[https://elifesciences.org/reviewed-preprints/99702]
Here, we demonstrate an AlphaFold2 based framework combined with all-atom enhanced sampling molecular dynamics and induced fit docking, named AF2RAVE-Glide, 
to conduct computational model based small molecule binding of metastable protein kinase conformations, initiated from protein sequences. 
We demonstrate the AF2RAVE-Glide workflow on three different protein kinases and their type I and II inhibitors, 
with special emphasis on binding of known type II kinase inhibitors which target the metastable classical DFG-out state. 
These states are not easy to sample from AlphaFold2. Here we demonstrate how with AF2RAVE these metastable conformations can be sampled for different kinases 
with high enough accuracy to enable subsequent docking of known type II kinase inhibitors with more than 50\% success rates across docking calculations. 
We believe the protocol should be deployable for other kinases and more proteins generally.

## Components

Data for this manuscript can be found [in this Drive](https://drive.google.com/drive/folders/1hSsnhNsF2uXZQ7SN6DxEgcLxDCHmVECj?usp=sharing)

In folder scripts
* `US_utils.py` - functions required to run this instance of our methodology
* `basicmd.py` - to run basic unbiased MD
* `equilprotocol.py` - to run equilibration protocol as in the paper
* `distanceUS.py` - to run a umbrella sampling simulation using OP computed with SPIB with distances as input CVs
* `kinaseCVs.py` - to extract the CVs we use
* `CalcLigRMSD.py` - to find the max common substructure of two ligands and compute the ligand RMSD

## Citation

Please cite the following reference if using this protocol with or without the provided code:

* "AlphaFold2-RAVE: From sequence to Boltzmann ensemble"
Bodhi P. Vani, Akashnathan Aranganathan, Dedi Wang, Pratyush Tiwary
J. Chem. Theory Comput. 2023; doi: https://doi.org/10.1021/acs.jctc.3c00290

* "Exploring kinase asp-phe-gly (dfg) loop conformational stability with alphafold2-rave." Vani BP, Aranganathan A, Tiwary P.  Journal of chemical information and modeling. 2023 Nov 20;64(7):2789-97. https://doi.org/10.1021/acs.jcim.3c01436

* "State predictive information bottleneck", Dedi Wang and Pratyush Tiwary, J. Chem. Phys. 154, 134111 (2021) https://doi.org/10.1063/5.0038198

* "Empowering AlphaFold2 for protein conformation selective drug discovery with AlphaFold2-RAVE." Gu X, Aranganathan A, Tiwary P.  ArXiv. 2024 Apr 10. https://elifesciences.org/reviewed-preprints/99702
