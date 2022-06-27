%make sure CobraToolbox and Gurobi solver is installed
initCobraToolbox(false);
gene_expression_matrix_file= '../data/no_entrezduplicates (1).xlsx';
% loading expression matrix
T=readtable(gene_expression_matrix_file);

% Enter different threshold values to test
% building a iMAT model
lb = 30;% percentile value
ub = 70;

% a=1 for iMAT model and a=0 for GIMME model
a=1;
[iMAT_model_normal_30n70_minmax, iMAT_model_PCOSnonIR_30n70_minmax, iMAT_model_PCOSIR_30n70_minmax] = PCOSmodel1_control(a,T, lb, ub);

%give a model from PCOSmodel1_control for sanity checks
model= iMAT_model_PCOSIR_30n70_minmax;%default
[TableChecks]=sanitychck(model)

%Flux balance analysis for models from PCOSmodel1_control function
[fba_30n70_normal_minmax,fba_30n70_PCOSnonIR_minmax,fba_30n70_PCOSIR_minmax]=FBA(iMAT_model_normal_30n70_minmax,iMAT_model_PCOSnonIR_30n70_minmax,iMAT_model_PCOSIR_30n70_minmax);

%Flux Variability Analysis for models from PCOSmodel1_control function
[fva_30n70_normal_minmax,fva_30n70_PCOSnonIR_minmax,fva_30n70_PCOSIR_minmax]=FVA(iMAT_model_normal_30n70_minmax,iMAT_model_PCOSnonIR_30n70_minmax,iMAT_model_PCOSIR_30n70_minmax);

%Flux Enrichment Analysis for models from PCOSmodel1_control function
%T here is the list of reaction names from Flux span ratio analysis
%By default the resultCell has output from iMAT derived model iMAT_model_PCOSIR_30n70_minmax
[resultCell]= FEACode(T,model);

%Displaying results from models
for v= [iMAT_model_normal_30n70_minmax, iMAT_model_PCOSnonIR_30n70_minmax, iMAT_model_PCOSIR_30n70_minmax]
    result_contextspecificmodel.gene_no=size(v.genes);
    result_contextspecificmodel.mets_no=size(v.mets);
    result_contextspecificmodel.rxns_no=size(v.rxns);
    result_contextspecificmodel
end

%Displaying biomass flux values from Flux Balance Analysis
for x= [fba_30n70_normal_minmax,fba_30n70_PCOSnonIR_minmax,fba_30n70_PCOSIR_minmax]
    result_fba= x.f;
    result_fba
end