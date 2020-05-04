function [real_out, perm_out]=lme_perm(table,predictor,formula,totperm)
% table=GO_table;
% 
% formula='RT~1+pred+Task+(1|SubID)';
% permvars={'Task','SubID'};
% totperm=100;

% run real model
eval(sprintf('table.pred=table.%s;',predictor));
model= fitlme(table,formula);
real_out=[double(model.Coefficients(3,2)) double(model.Coefficients(3,4)) double(model.Coefficients(3,6))];

uniqueIDs=unique(table.SubID);
uniqueTasks=unique(table.Task);
perm_out=nan(totperm,4);
fprintf('%4.0f/%4.0f\n',0,totperm)
for np=1:totperm
    pred_perm=table.pred;
    for nT=1:length(uniqueTasks)
        for nS=1:length(uniqueIDs)
            idx=find(table.SubID==uniqueIDs(nS) & table.Task==uniqueTasks(nT));
            temp=pred_perm(idx);
            pred_perm(idx)=temp(randperm(length(temp)));
        end
    end
    
    table2=table;
    table2.pred=pred_perm;
    model= fitlme(table2,formula);
    perm_out(np,:)=[double(model.Coefficients(3,2)) double(model.Coefficients(3,4)) double(model.Coefficients(3,6)) np];
    fprintf('\b\b\b\b\b\b\b\b\b\b%4.0f/%4.0f\n',np,totperm)
end
fprintf('\n');