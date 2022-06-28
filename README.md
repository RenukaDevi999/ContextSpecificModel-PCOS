# ContextSpecificModel-PCOS
The ContextSpecificModel-PCOS algorithm takes in gene expression data for three cases: Control, PCOS without insulin resistance and PCOS with insulin resistance and creates a disease metabolic network model. It can also determine differentially expressed genes,reactions, biomass flux, metabolic fluxes and reaction subsystems. MATLAB, cobra toolbox and gurobi solver softwares are required to execute the code.

RUNNING ContextSpecificModel-PCOS

To run ContextSpecificModel-PCOS, save the reaction gene expression data in excel format and the relevant reaction list names that are obtained after performing Flux Span Ratio into the "data" folder.Then run main.m from the "Main" folder. Refer to the Folder information given below for more information. 

FOLDER INFORMATION 

Data: 
The data folder consists of 2 excel files. First one consists of gene expression data obtained from NCBI GEO datasets where there are four columns: Gene number that is in obtained by mapping genes to RECON 3D genes, entrez ID values for genes, the corresponding expression values for control, expression values for PCOS without insulin resistance and expression values for PCOS with insulin resistance. The second dataset is the list of reactions for PCOS with insulin resistance (in the form of identifiers) obtained after Flux span ratio analysis.
 
Main:

Main.m is the file to be used to run Contextspecificmodel-PCOS. Intially the context specific model is built using iMAT algortihm(default). So, 3 models are obtained representing 3 cases: Control, PCOS without insulin resistance and PCOS with insulin resistance.Then using the model, sanity checks are performed. To the model, flux balance analysis is performed following which flux variability analysis is performed.Finally Flux enrichment analysis is performed for the given list of reactions in data. By default, the reaction list is for PCOS with insulin resistance.



MATLAB FILE INFORMATION

Given below are discriptions for various ".m files" used in the repository. 


File: main.m
Cobra toolbox and gurobi solver are required to run this file. RECON3DModel_301.mat is used to map the genes to reactions in the process of building the model. By default, iMAT models are built.


Outputs:

	•	result.mat with vectors:
	⁃	mets: The standard reaction Gibbs free energy estimates obtained from data reconciliation
	⁃	genes: The standard formation Gibbs free energy estimates obtained from data reconciliation
	⁃	rxns: the CIDS corresponding to the compounds for which the modified group contribution method could not provide estimates ( due to insufficient measured data )
	⁃	f: the RIDS corresponding to the reactions for which  the modified group contribution method could not provide estimates ( due to insufficient measured data )
	- hellooo
	⁃	resultCell : all RIDS present in the reaction set

File: parseKeggModel.m (&  reaction2sparse.m)

Parses the reaction file with each line representing a reaction in KEGG format to obtain S, CIDS and RIDS

Inputs:

	•	ReactionStrings - cell array of each reaction in KEGG format
	•	arrow - arrow in the reaction (default is '=')

Outputs:

	•	S - Stoichiometric matrix
	•	cids - CIDS of all compounds in the reaction set
	•	rids - RIDS of all compounds in the reaction set

File: legendretransformF.m 

Applies the inverse Legendre transform to the transformed formation Gibbs energies. Maxstar.m and transform.m are used to perform the transform. 

Inputs:
 
	•	m - number of formation energies
	•	kegg_pKa - contains the relevant pKa, nH, charge etc information for each compound
	•	cids_m - CIDS of the formation energies that are to be inverse transformed
	•	pH - vector containing the pH values for each apparent Gibbs energy of formation measurement.
	•	T - vector containing the temperature values for each apparent Gibbs energy of formation measurement.
	•	I - vector containing the ionic strength values for each apparent Gibbs energy of formation measurement.
	•	dG0f_prime - the transformed Gibbs free energy of formation value that is to be converted to standard form ( or "chemical energy".

Outputs:

	•	dG0 - The measured standard Gibbs free energies of formation obtained after performing the inverse Legendre transform
	•	reverse_ddG0s - The difference between the standard  transformed ( or apparent) Gibbs free energies of formation and theGibbs free energies of formation.


File: legendretransformR.m 

Applies the inverse Legendre transform to the transformed reaction Gibbs energies. Maxstar.m and transform.m are used to perform the transform. 

Inputs:
	 
	•	S - stoichiometric matrix
	•	kegg_pKa - contains the relevant pKa, nH, charge etc information for each compound
	•	cids_m - CIDS of the formation energies that are to be inverse transformed
	•	pH - vector containing the pH values for each apparent Gibbs energy of reaction measurement.
	•	T - vector containing the temperature values for each apparent Gibbs energy of reaction measurement.
	•	I - vector containing the ionic strength values for each apparent Gibbs energy of reaction measurement.
	•	dG0r_prime - the transformed Gibbs free energy of reaction value that is to be converted to standard form 

Outputs:

	•	dG0 - The measured standard Gibbs free energies of reaction obtained after performing the inverse Legendre transform
	•	reverse_ddG0s - The difference between the standard  transformed ( or apparent) Gibbs free energies of reaction and theGibbs free energies of reaction.


File: recon_l.m

This file is used to perform the linear data reconciliation for Gibbs free energies. It calculates reconciled estimates for all measured Gibbs energies and computes prediction estimates for the unmeasured observable Gibbs energies. If prompted, it imputes the group contribution estimates for the unobservable Gibbs energies. It requires Gibbs free energy measurements, thermodynamic constraints (Stoichiometric matrix S), variance matrix (Sigma) and indices of the measured reaction and formation Gibbs energies corresponding to the columns and rows of S.  

Inputs: 
	•	 S - The (m x n) stoichiometric matrix for the n reactions between m compounds obtained from the reactions.txt file that represents the thermodynamic constraints for the available data 
	•	 y - vector of consisting of all measured reaction and formation Gibbs energies.
	•	index_rm -  indices of the "n" reaction Gibbs energies corresponding to the columns of S that are measured
	•	index_fm -  indices of the "m" formation Gibbs energies corresponding to the rows of S that are measured

Outputs:
	•	Recon_var= A vector of dimensions (m+n) x 1 representing the reconciled estimates. When the estimate is not available, the corresponding value is "NaN"
	•	index_unobservable= indices from the range 1:(m+n) that represent the unobservable Gibbs energies. 

File: Groupcontributions.m

This file is used to perform the modified group contribution method. By using all the available measured data, the thermodynamic constraints between them, estimates are obtained for all the group contributions using the constrained optimisation function fmincon. 

Inputs:

	•	A -  the reduced constraint matrix representing the thermodynamic constraints between the measured variables (obtained from the QR decomposition of S)
	•	G_r - measured reaction Gibbs energies
	•	G_f - measured formation Gibbs energies
	•	G - the group decomposition matrix
	•	index_fest - the indices of the unmeasured formation energies

File: costFunction.m
The optimisation function for fmincon in the group contribution method

Inputs:

	•	x - vector consisting of reconciled estimates of reaction, formation Gibbs energies and group contributions.
	•	b - vector consisting of measured reaction Gibbs energies
	•	Gf_m - vector consisting of measured formation Gibbs energies
	•	G_m - group decomposition matrix for the measured formation Gibbs energies

Outputs:

	•	J - The average sum squared error between measurements and estimates






