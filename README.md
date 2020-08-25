# BindingFreeEnergy

Calculate the binding free energy based on the method developed by [James C. Gumbart et al., 2013](https://doi.org/10.1021/ct400273t), [James C. Gumbart et al., 2014](https://pubs.acs.org/doi/10.1021/ct3008099) & [Hong Zhang et al., 2018](https://www.mdpi.com/1420-3049/23/2/228).

The equilibium constant is calculating using the following equation:
![](equations.png "output" )

Picture took from the tutotial:[Protein:ligand standard binding free energies: A tutorial for alchemical and geometrical transformations](https://www.ks.uiuc.edu/Training/Tutorials/namd/PLB/tutorial-protein-ligand.pdf)


Once all your simulations done, you can use this script to calculate the binding free energy. Please change the parameters according to your simulation parameters (force constant, temperature, ...).


##output example (file: binding_free_energy.nrj):

#Contribution	numerator(N)	denominator(D)	Ratio(N/D)	Energy (kcal/mol)
RMSDb 	 x 	 x 	 x 	 x
Theta 	 x	 x 	 x 	 x
Phi 	 x 	 x 	 x 	 x
Psi 	 x 	 x 	 x 	 x
theta (PA) 	 x 	 x 	 x 	 x
phi (PA) 	 x 	 x 	 x x
Theta/Phi/Psi_bulk 	 - 	 - 	 - x
I* 	 - 	 - 	 - x
S* 	 - 	 - 	 - x
RMSDf 	 x 	 x 	 x x

K_eq =  x

G_bind =  x kcal/mol
