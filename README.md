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

	•	iMAT_model_normal_30n70_minmax.mat with struct field:
	⁃	mets: The metabolites of the control model for threshold 30 and 70(default threshold).
	⁃	genes: The genes involved in the control model for threshold 30 and 70.
	⁃	rxns: The reactions of the control model for threshold 30 and 70.
	•	iMAT_model_PCOSnonIR_30n70_minmax.mat with struct field:
	⁃	mets: The metabolites of the PCOS without insulin resistance disease model for threshold 30 and 70(default threshold).
	⁃	genes: The genes involved in the PCOS without insulin resistance disease model for threshold 30 and 70.
	⁃	rxns: The reactions of the PCOS without insulin resistance disease model for threshold 30 and 70.
	•	iMAT_model_PCOSIR_30n70_minmax.mat with struct field:
	⁃	mets: The metabolites of the PCOS wit insulin resistance disease model for threshold 30 and 70(default threshold).
	⁃	genes: The genes involved in the PCOS with insulin resistance disease model for threshold 30 and 70.
	⁃	rxns: The reactions of the PCOS with insulin resistance disease model for threshold 30 and 70.

	•	fba_30n70_normal_minmax.mat with values:
	⁃	f: The biomass flux of the control for the 'iMAT_model_normal_30n70_minmax' model .
	•	fba_30n70_PCOSnonIR_minmax.mat with values:
	⁃	f: The biomass flux of PCOS without insulin resistance disease for the input 'iMAT_model_PCOSnonIR_30n70_minmax' model.
	•	fba_30n70_PCOSIR_minmax.mat with values:
	⁃	f: The biomass flux of  PCOS with insulin resistance disease for the  input 'iMAT_model_PCOSIR_30n70_minmax' model.

	•	fva_30n70_normal_minmax.mat with vector:
	⁃	minimum and maximum flux of reactions of control model.
	•	fva_30n70_PCOSnonIR_minmax.mat with vector:
	⁃	minimum and maximum flux of reactions of PCOS without insulin resistance model.
	•	fva_30n70_PCOSIR_minmax.mat with vector:
	⁃	minimum and maximum flux of reactions of PCOS with insulin resistance model.

	•	TableChecks: is a table with 16 fields of test results for a model.

	⁃	resultCell : consists of adjusted P values, reaction subsystems, total set size(number of total reactions), enriched set(total number of differentially expressed reactions) for the input reaction list excel file.


File: PCOSmodel1_control.m

Creates context specific model for control, PCOS with/ without insulin resistance using iMAT(default) and GIMME. 

Inputs:

	•	T- table with gene expression data for control, PCOS with/without insulin resistance.
	•   For iMAT models,
	⁃   lb- lower threshold.
	⁃   ub- upper threshold.
	•   For GIMME models,
	⁃   lb- lower threshold,
	•   a=1 for iMAT models and a=0 for GIMME models.


Outputs:

	•	iMAT_model_normal_30n70_minmax.mat with struct field:
	⁃	mets: The metabolites of the control model for threshold 30 and 70(default threshold).
	⁃	genes: The genes involved in the control model for threshold 30 and 70.
	⁃	rxns: The reactions of the control model for threshold 30 and 70.
	•	iMAT_model_PCOSnonIR_30n70_minmax.mat with struct field:
	⁃	mets: The metabolites of the PCOS without insulin resistance disease model for threshold 30 and 70(default threshold).
	⁃	genes: The genes involved in the PCOS without insulin resistance disease model for threshold 30 and 70.
	⁃	rxns: The reactions of the PCOS without insulin resistance disease model for threshold 30 and 70.
	•	iMAT_model_PCOSIR_30n70_minmax.mat with struct field:
	⁃	mets: The metabolites of the PCOS wit insulin resistance disease model for threshold 30 and 70(default threshold).
	⁃	genes: The genes involved in the PCOS with insulin resistance disease model for threshold 30 and 70.
	⁃	rxns: The reactions of the PCOS with insulin resistance disease model for threshold 30 and 70.


File: sanitychck.m 

Performs model checks like leak test, single deletion functions etc.

Inputs:

	•	model: iMAT_model_PCOSIR_30n70_minmax model(default)


Outputs:

	•	TableChecks: is a table with 16 fields of test results for a model.

File: FBA.m
Applies flux balance analysis for the context specific models.

Inputs:

	•	iMAT_model_normal_30n70_minmax model
	•	iMAT_model_PCOSnonIR_30n70_minmax model
	•	iMAT_model_PCOSIR_30n70_minmax model


Outputs:

	•	fba_30n70_normal_minmax.mat with values:
	⁃	f: The biomass flux of the control for the 'iMAT_model_normal_30n70_minmax' model .
	•	fba_30n70_PCOSnonIR_minmax.mat with values:
	⁃	f: The biomass flux of PCOS without insulin resistance disease for the input 'iMAT_model_PCOSnonIR_30n70_minmax' model.
	•	fba_30n70_PCOSIR_minmax.mat with values:
	⁃	f: The biomass flux of  PCOS with insulin resistance disease for the  input 'iMAT_model_PCOSIR_30n70_minmax' model.

File: FVA.m
Applies flux variability analysis for the context specific models.

Inputs:

	•	iMAT_model_normal_30n70_minmax model
	•	iMAT_model_PCOSnonIR_30n70_minmax model
	•	iMAT_model_PCOSIR_30n70_minmax model


Outputs:

	•	fva_30n70_normal_minmax.mat with vector:
	⁃	minimum and maximum flux of reactions of control model.
	•	fva_30n70_PCOSnonIR_minmax.mat with vector:
	⁃	minimum and maximum flux of reactions of PCOS without insulin resistance model.
	•	fva_30n70_PCOSIR_minmax.mat with vector:
	⁃	minimum and maximum flux of reactions of PCOS with insulin resistance model.


File: FEACode.m

This file is used to perform flux enrichment analysis for the results from Flux span ratio analysis. 

Inputs: 
	•	 T- reaction names table
	•    model: iMAT_model_PCOSIR_30n70_minmax model(default)

Outputs:
    •    resultCell : consists of adjusted P values, reaction subsystems, total set size(number of total reactions), enriched set(total number of differentially expressed reactions) for the input reaction list excel file.



File: Recon3DModel_301.m

This file has information on 4,140 metabolites, 13,543 Reactions of 3,288 genes. This is used to build models in 'PCOSmodel1_control.m'






