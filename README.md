# BindingFreeEnergy

Calculate the binding free energy based on the method initiated by [Woo and Roux, 2005](https://doi.org/10.1073/pnas.0409005102) and generalised by  [James C. Gumbart et al., 2013](https://pubs.acs.org/doi/10.1021/ct3008099).

The equilibium constant is calculating using the following equation:
![](equations.png "output" )

Picture took from the tutotial:[Protein:ligand standard binding free energies: A tutorial for alchemical and geometrical transformations](https://www.ks.uiuc.edu/Training/Tutorials/namd/PLB/tutorial-protein-ligand.pdf)


Once all your simulations done, you can use this script to calculate the binding free energy. Please change the parameters according to your simulation parameters (force constant, temperature, ...).


All the results (each contribution and the final binding free energy) will be written in the file : "binding_free_energy.nrj"
