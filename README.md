# JPH-GBT model
The codes are for the implementation of paper "Join persistent homology (JPH)-based machine learning model for metalloprotein-ligand binding affinity prediction"

Join persistent homology (JPH) involves creating a multistage filtration for the join of the original simplicial complex and a specially designed test simplicial complex (see Methods). In our JPH-GBT model, we use the bipartite graph (BG) to model the interactions for metalloprotein-ligand complexes, and the test simplicial complex $L$ is considered as $\{w\}$. The join persistent homology is constructed by combining the two filtration processes together through the join of bipartite graph and $\{w\}$. We consider two specially-designed filtration functions,
\begin{gather*}
\mathcal{F}^0:~BG^0_{f_0}\hookrightarrow BG^0_{f_1}\hookrightarrow \cdots \hookrightarrow BG^0_{f_{n-1}}\hookrightarrow BG^0_{f_n}=BG;\\
\mathcal{F}^1:~BG^1_{g_0}\hookrightarrow BG^1_{g_1}\hookrightarrow \cdots \hookrightarrow BG^1_{g_{m-1}}\hookrightarrow BG^1_{g_m}=BG.
\end{gather*}
where $(f_0, f_1 , \cdots,  f_n)$ and $(g_0 , g_1, \cdots, g_m)$ are filtration values.  

# How to use 
You can replace the file path by your own and run the JPH-GBT-Train.py and JPH-GBT-Test.py to get the train and test features, and then run the JPH-GBT-Model.py to get the result.
