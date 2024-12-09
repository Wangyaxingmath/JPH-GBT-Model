
# JPH-GBT model
The codes are for implementing the paper "Join persistent homology (JPH)-based machine learning model for metalloprotein-ligand binding affinity prediction."

## Join persistent homology (JPH)
Join Persistent Homology (JPH) employs a set of filtration functions to generate a multistage filtration for the join of the original simplicial complex and a specially-designed test simplicial complex. For instance, we consider four filtration functions on a pentagon K, where the specially-esigned test simplicial complex is regarded as {w1, w2, w1w2}. An illustration of this is shown in the figure below.

![F004](https://github.com/user-attachments/assets/dcce1405-f8bf-455d-b391-6eb5f1098722)

## Metalloprotein-ligand complex dataset
We need to extract the atomic coordinates from the metalloprotein-ligand dataset provided by PBDBind-v2020. This dataset is organized into training, validation, and test sets based on various types of metal ions, as described by Jiang et al. (2023). For more details, please refer to the Supplementary Files of the article titled "MetalProGNet: A Structure-Based Deep Graph Model for Metalloprotein–Ligand Interaction Predictions" by Jiang et al.

## Molecular element-specific representation and JPH featurization
For each metalloprotein-protein complex, the element-specific atom combinations are considered. The types of elements considered for metalloprotein at the binding core regions are categorized into 15 groups based on their types, including five non-metal ions C, N, O, S, H, and ten metal ions K, Zn, Ca, Na, Mg, Mn, Fe, Cu, Ni, Co. The ligand atoms are decomposed into ten groups based on their types of elements, including C, N, O, S, H, P,  F, Cl, Br, and I. The interactions between metalloproteins and ligands are categorized into protein-ligand interactions, metal ions-protein interactions, metal ions-ligand interactions, respectively, with an illusration shown in the following Figure. In this work, we consider the binning approach.

![F003](https://github.com/user-attachments/assets/24a6afd6-1783-41a5-8552-3e07b2d1330f)


## How to use ?
You can modify the file path and then run JPH-GBT-Train.py and JPH-GBT-Test.py to generate the train and test features. After that, run JPH-GBT-Model.py to obtain the results. 
