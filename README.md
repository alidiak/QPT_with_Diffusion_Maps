# QPT_with_Diffusion_Maps
Quantum Phase Transition Discovery with Diffusion Maps. Supporting code for:
https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.225701. Unfortunately, most of the data files for the large system sizes used to generate the figures were too large to upload, but these are the codes that generated them. 

The Exact_Diagonalization folder contains ED simulations of quantum-many-body models and will use the generated samples to perform diffusion maps (among other unsupervised dimension reduction and clustering techniques such as PCA and autoencoders with k-means). The MPS_Postprocessing uses a set of optimized MPS, draws samples, and performs diffusion maps. The Utils folder contains the supporting functions used for ED, MPS processing, dimension reduction, clustering, and plotting.

The Edgemats folder contains a Mathematica notebook that generates the adjacency matrices for certain lattices. 
