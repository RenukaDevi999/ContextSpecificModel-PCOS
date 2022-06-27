function [resultCell]= FEACode(T,model)
file= '../data/Rxns_Names_FEA.xlsx';
T= readtable(file);
A=cellstr(T.Rxns);
model=iMAT_model_PCOSIR_30n70_minmax;
X=findRxnIDs(model,A);
resultCell = FEA(model,X,'subSystems');
end