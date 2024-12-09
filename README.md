## Join persistent homology (JPH)
Join persistent homology (JPH) involves creating a multistage filtration for the join of the original simplicial complex and a specially-designed test simplicial complex. For example, we consider four filtration functions on pentagon K and the specially-designed test simplicial complex is considered as {w1 w2, w1w2}, with an illusration shown in the following Figure.

![F004](https://github.com/user-attachments/assets/dcce1405-f8bf-455d-b391-6eb5f1098722)

## JPH-GBT model
The codes are for the implementation of paper "Join persistent homology (JPH)-based machine learning model for metalloprotein-ligand binding affinity prediction"

## Metalloprotein-ligand complex dataset
We need to get the atom coordinates from the metalloprotein-ligand dataset from PBDBind-v2020. The dataset is organized into training, validation, and test sets based on the different types of metal ions, as outlined by Jiang et al.(2023). For further details, please refer to the Supplementary Files of the article: "MetalProGNet: a structure-based deep graph model for
metalloprotein–ligand interaction predictions, Jiang et al." 

## Molecular element-specific representation and JPH featurization
For each metalloprotein-protein complex, the element-specific atom combinations are considered. The types of elements considered for metalloprotein at the binding core regions are categorized into 15 groups based on their types, including five non-metal ions C, N, O, S, H, and ten metal ions K, Zn, Ca, Na, Mg, Mn, Fe, Cu, Ni, Co. The ligand atoms are decomposed into ten groups based on their types of elements, including C, N, O, S, H, P,  F, Cl, Br, and I. The interactions between metalloproteins and ligands are categorized into protein-ligand interactions, metal ions-protein interactions, metal ions-ligand interactions, respectively, with an illusration shown in the following Figure. In this work, we consider the binning approach.

![F003](https://github.com/user-attachments/assets/24a6afd6-1783-41a5-8552-3e07b2d1330f)


## How to use 
You can replace the file path by your own and run the JPH-GBT-Train.py and JPH-GBT-Test.py to get the train and test features, and then run the JPH-GBT-Model.py to get the result.
