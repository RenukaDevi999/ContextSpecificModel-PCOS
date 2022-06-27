
% initCobraToolbox(false)
% loading expression matrix
% T=readtable("no_entrezduplicates (1).xlsx");
% float -> char array -> cell -> removing space
function [iMAT_model_normal_30n70_minmax, iMAT_model_PCOSnonIR_30n70_minmax, iMAT_model_PCOSIR_30n70_minmax] = PCOSmodel1_control(a,T, lb, ub)
if a==1
    T.gene_number = strrep(cellstr(num2str(T.gene_number)),' ','');
    % loading recon model
    load('Recon3DModel_301.mat')
    Recon3DModel='../data/Recon3DModel_301.mat';
    parsedGPR = GPRparser(Recon3DModel);
    % use sum for OR / max for OR
    minSum = 0;
    
    
    rxn_exp= selectGeneFromGPR(Recon3DModel, T.gene_number, T.avg_N, parsedGPR, minSum);
    rxn_exp_PCOSnonIR = selectGeneFromGPR(Recon3DModel, T.gene_number, T.avg_nonIR, parsedGPR, minSum);
    rxn_exp_PCOSIR = selectGeneFromGPR(Recon3DModel, T.gene_number, T.avg_IR, parsedGPR, minSum);
    
    
    rxn_exp(isnan(rxn_exp)) = -1;
    rxn_exp_PCOSnonIR(isnan(rxn_exp_PCOSnonIR)) = -1;
    rxn_exp_PCOSIR(isnan(rxn_exp_PCOSIR)) = -1;
    
    
    % building a iMAT model
    % lb = 30;% percentile value
    % ub = 70;
    
    % get thresholds for the given percentiles ignoring the zeros
    threshold_normal_lb = prctile(rxn_exp(rxn_exp>0), lb,'ALL');
    threshold_normal_ub = prctile(rxn_exp(rxn_exp>0), ub, 'ALL');
    threshold_PCOSnonIR_lb = prctile(rxn_exp_PCOSnonIR(rxn_exp_PCOSnonIR>0), lb,'ALL');
    threshold_PCOSnonIR_ub = prctile(rxn_exp_PCOSnonIR(rxn_exp_PCOSnonIR>0),ub, 'ALL');
    threshold_PCOSIR_lb = prctile(rxn_exp_PCOSIR(rxn_exp_PCOSIR>0), lb,'ALL');
    threshold_PCOSIR_ub = prctile(rxn_exp_PCOSIR(rxn_exp_PCOSIR>0), ub, 'ALL');
    
    iMAT_model_normal_30n70_minmax = removeUnusedGenes(iMAT(Recon3DModel, rxn_exp, threshold_normal_lb, threshold_normal_ub,1e-5,Recon3DModel.rxns(find(Recon3DModel.c))));
    iMAT_model_PCOSnonIR_30n70_minmax =removeUnusedGenes(iMAT(Recon3DModel, rxn_exp_PCOSnonIR, threshold_PCOSnonIR_lb, threshold_PCOSnonIR_ub,1e-5,Recon3DModel.rxns(find(Recon3DModel.c))));
    iMAT_model_PCOSIR_30n70_minmax =removeUnusedGenes(iMAT(Recon3DModel, rxn_exp_PCOSIR , threshold_PCOSIR_lb, threshold_PCOSIR_ub,1e-5,Recon3DModel.rxns(find(Recon3DModel.c))));
    %save('iMAT_20n80_minmax','iMAT_model_normal_20n80_minmax','iMAT_model_PCOSnonIR_20n80_minmax','iMAT_model_PCOSIR_20n80_minmax')
elseif a==0
    T.gene_number = strrep(cellstr(num2str(T.gene_number)),' ','');
    % loading recon model
    Recon3DModel='../data/Recon3DModel_301.mat';
    parsedGPR = GPRparser(Recon3DModel);
    % use sum for OR / max for OR
    minSum = 0;
    
    
    rxn_exp= selectGeneFromGPR(Recon3DModel, T.gene_number, T.avg_N, parsedGPR, minSum);
    rxn_exp_PCOSnonIR = selectGeneFromGPR(Recon3DModel, T.gene_number, T.avg_nonIR, parsedGPR, minSum);
    rxn_exp_PCOSIR = selectGeneFromGPR(Recon3DModel, T.gene_number, T.avg_IR, parsedGPR, minSum);
    
    
    rxn_exp(isnan(rxn_exp)) = -1;
    rxn_exp_PCOSnonIR(isnan(rxn_exp_PCOSnonIR)) = -1;
    rxn_exp_PCOSIR(isnan(rxn_exp_PCOSIR)) = -1;
    % building a GIMME model
    % lb = 30;% percentile value
    
    % get thresholds for the given percentiles ignoring the zeros
    threshold_normal_lb = prctile(rxn_exp(rxn_exp>0), lb,'ALL');
    threshold_PCOSnonIR_lb = prctile(rxn_exp_PCOSnonIR(rxn_exp_PCOSnonIR>0), lb,'ALL');
    threshold_PCOSIR_lb = prctile(rxn_exp_PCOSIR(rxn_exp_PCOSIR>0), lb,'ALL');
    
    GIMME_model_normal_30_minmax = removeUnusedGenes(GIMME(Recon3DModel, rxn_exp, threshold_normal_lb,Recon3DModel.rxns(find(Recon3DModel.c))));
    GIMME_model_PCOSnonIR_30_minmax =removeUnusedGenes(GIMME(Recon3DModel, rxn_exp_PCOSnonIR, threshold_PCOSnonIR_lb,Recon3DModel.rxns(find(Recon3DModel.c))));
    GIMME_model_PCOSIR_30_minmax =removeUnusedGenes(GIMME(Recon3DModel, rxn_exp_PCOSIR , threshold_PCOSIR_lb,Recon3DModel.rxns(find(Recon3DModel.c))));
end

end
