%%
clear all;
close all;
run ../localdef_wanderIM.m

addpath(genpath(path_export))
addpath(genpath(path_RainCloudPlot))
addpath(genpath(lscpTools_path))

%% Initialize variables.
recompute=0;
np=0;
totperm=1000;
Elec_Labels=layout.label;
Task_Labels={'Face','Digit'};
Mat_Param=[];
vBias_est=cell(1,2);
if recompute

for nE=1:63
    fprintf('... ... Electrode %s\n',Elec_Labels{nE})
    filename=(['/Users/tand0009/Data/WanderIM/hddm/Models_20s/model_W_' Elec_Labels{nE} filesep 'model_stats_W_' Elec_Labels{nE} '.csv']);
    if exist(filename)~=0
        
    else
        fprintf('... ... Electrode %s MISSING\n',Elec_Labels{nE})
        continue;
    end
    delimiter = ',';
    startRow = 2;
    formatSpec = '%s%f%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    
    model = table;
    model.VarName1 = cellstr(dataArray{:, 1});
    model.mean = dataArray{:, 2};
    model.std = dataArray{:, 3};
    model.q = dataArray{:, 4};
    model.q1 = dataArray{:, 5};
    model.q2 = dataArray{:, 6};
    model.q3 = dataArray{:, 7};
    model.q4 = dataArray{:, 8};
    model.mcerr = dataArray{:, 9};
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    
    %%%
    var_2factor={'v'};    
    mat_LME=[];
    for n=1:length(var_2factor)
        rows=find(~cellfun(@isempty,regexp(model.VarName1,sprintf('^%s_subj*',var_2factor{n}))));
        rownames=strvcat(model.VarName1);
        mat_LME{n}(:,1)=model.mean(rows);
        mat_LME{n}(:,2)=str2num(rownames(rows,8));
        mat_LME{n}(:,3)=str2num(rownames(rows,10));
        mat_LME{n}(:,4)=str2num(rownames(rows,12));
        mat_LME{n}(:,5)=str2num(rownames(rows,end-2:end));
        Mat_Param=[Mat_Param ; [repmat([1 nE],size(mat_LME{n},1),1) mat_LME{n}]];
    end
    
    
    uniqueSubs=unique(mat_LME{1}(:,5));
    mat_vbias=[];
    for nS=1:length(uniqueSubs)
        for nTask=1:2
            for nSW=0:1
                tempA=mat_LME{1}(mat_LME{1}(:,5)==uniqueSubs(nS) & mat_LME{1}(:,3)==nSW & mat_LME{1}(:,2)==nTask  & mat_LME{1}(:,4)==1,1);
                tempB=mat_LME{1}(mat_LME{1}(:,5)==uniqueSubs(nS) & mat_LME{1}(:,3)==nSW & mat_LME{1}(:,2)==nTask  & mat_LME{1}(:,4)==0,1);
                if ~isempty(tempA) && ~isempty(tempB)
                    mat_vbias=[mat_vbias ; [tempA+tempB nTask nSW uniqueSubs(nS)]];
                end
            end
        end
    end
    tbl_bv=array2table(mat_vbias,'VariableNames',{'bv','Task','sw','SubID'});
    tbl_bv.Task=categorical(tbl_bv.Task);
    tbl_bv.SubID=categorical(tbl_bv.SubID);
    
    
    tic;
    [real_out, perm_out]=lme_perm(tbl_bv,'sw','bv~1+pred+Task+(1|SubID)',totperm);
    vBias_est{1}=[vBias_est{1} ; [nE real_out]];
    vBias_est{2}=[vBias_est{2} ; [nE*ones(totperm,1) perm_out]];
    toc;
    fprintf('\n')
    
    
end
save('/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_HDDM_LMEandPerm_vBias_20s.mat','vBias_est')
else
load('/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_HDDM_LMEandPerm_vBias_20s.mat');
end
%% Compute monte cartlo p
clus_alpha=0.025;
montecarlo_alpha=0.05/6;

addpath((path_fieldtrip)); ft_defaults;
cfg_neighb=[];
cfg_neighb.method = 'tri';
cfg_neighb.layout=path_PsychFTlayout;
neighbours = ft_prepare_neighbours(cfg_neighb);

[vBias_clus]=get_clusterperm_lme(vBias_est,clus_alpha,montecarlo_alpha,totperm,neighbours);

%%
cmap=cbrewer('div','RdBu',256); cmap=flipud(cmap);
Titles={'vBias'};
for nP=1:length(Titles)
    figure;
    if nP==1
        temp_topo=vBias_est{1}(:,3);
        temp_clus=vBias_clus; 
    end
    % temp_topo(match_str(layout.label,{'TP9','TP10'}))=0;
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    colorbar; caxis([-1 1]*6); %caxis([-1 1]*max(abs(temp_topo)))
        colormap(cmap);

    load(path_PsychFTlayout);
    if ~isempty(temp_clus)
        for nclus=1:length(temp_clus)
        ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
        fprintf('... ... found %s cluster (%g) of %g electrodes (tval cluster=%g, Pmc=%g)\n',temp_clus{nclus}{1},nclus,length(temp_clus{nclus}{2}),temp_clus{nclus}{3},temp_clus{nclus}{4})
        end
    end
    title([Titles{nP}])
%     export_fig([path_fig filesep 'LocalSleep_HDDM_LMEbyTrial_' Titles{nP} '.fig'])
%     export_fig([path_fig filesep 'LocalSleep_HDDM_LMEbyTrial_' Titles{nP} '.eps'],'-r 300')
end
