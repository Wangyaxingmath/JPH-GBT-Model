## Join persistent homology (JPH)
Join persistent homology (JPH) involves creating a multistage filtration for the join of the original simplicial complex and a specially designed test simplicial complex. JPH can provide a more comprehensive characterization of the intrinsic topological information. Details can be found in Methods of the paper.

## JPH-GBT model
The codes are for the implementation of paper "Join persistent homology (JPH)-based machine learning model for metalloprotein-ligand binding affinity prediction"

## Metalloprotein-ligand complex dataset
We need to get the atom coordinates from the metalloprotein-ligand dataset from PBDBind-v2020. The dataset is organized into training, validation, and test sets based on the different types of metal ions, as outlined by Jiang et al.(2023). For further details, please refer to the Supplementary Files of the article: "MetalProGNet: a structure-based deep graph model for
metalloprotein–ligand interaction predictions, Jiang et al." 

## Algebraic representation and persistent Tor-algebra featurization

# How to use 
You can replace the file path by your own and run the JPH-GBT-Train.py and JPH-GBT-Test.py to get the train and test features, and then run the JPH-GBT-Model.py to get the result.
