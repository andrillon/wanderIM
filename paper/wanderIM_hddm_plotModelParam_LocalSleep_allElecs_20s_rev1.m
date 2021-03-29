%%
clear all;
close all;
run ../localdef_wanderIM.m

addpath(genpath(path_export))
addpath(genpath(path_RainCloudPlot))
addpath(genpath(lscpTools_path))

recompute=1;
%% Initialize variables.
np=0;
totperm=1000;
Elec_Labels=layout.label;
Task_Labels={'Face','Digit'};
Mat_Param=[];
a_est=cell(1,2);
t_est=cell(1,2);
z_est=cell(1,2);
vGO_est=cell(1,2);
vNOGO_est=cell(1,2);
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
    var_1factor={'a','t','z'};
    
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
    tbl_v=array2table(mat_LME{1},'VariableNames',{'v','Task','sw','stimtype','SubID'});
    tbl_v.stimtype=categorical(tbl_v.stimtype);
    tbl_v.Task=categorical(tbl_v.Task);
    tbl_v.SubID=categorical(tbl_v.SubID);
    
    tic;
    [real_out, perm_out]=lme_perm(tbl_v(tbl_v.stimtype=='1',:),'sw','v~1+pred+Task+(1|SubID)',totperm);
    vGO_est{1}=[vGO_est{1} ; [nE real_out]];
    vGO_est{2}=[vGO_est{2} ; [nE*ones(totperm,1) perm_out]];
    toc;
    fprintf('\n')
    [real_out, perm_out]=lme_perm(tbl_v(tbl_v.stimtype=='0',:),'sw','v~1+pred+Task+(1|SubID)',totperm);
    vNOGO_est{1}=[vNOGO_est{1} ; [nE real_out]];
    vNOGO_est{2}=[vNOGO_est{2} ; [nE*ones(totperm,1) perm_out]];
    
    %%%
    mat_LME2=[];
    for n=1:length(var_1factor)
        rows=find(~cellfun(@isempty,regexp(model.VarName1,sprintf('^%s_subj*',var_1factor{n}))));
        rownames=strvcat(model.VarName1);
        mat_LME2{n}(:,1)=model.mean(rows);
        mat_LME2{n}(:,2)=str2num(rownames(rows,8));
        mat_LME2{n}(:,3)=str2num(rownames(rows,10));
        mat_LME2{n}(:,4)=str2num(rownames(rows,13:15));
        Mat_Param=[Mat_Param ; [repmat([2+n nE],size(mat_LME2{n},1),1) mat_LME2{n}(:,1:3) nan(size(mat_LME2{n},1),1) mat_LME2{n}(:,4)]];
    end
    tbl_a=array2table(mat_LME2{1},'VariableNames',{'a','Task','sw','SubID'});
    tbl_a.SubID=categorical(tbl_a.SubID);
    tbl_a.Task=categorical(tbl_a.Task);
    
    tbl_t=array2table(mat_LME2{2},'VariableNames',{'t','Task','sw','SubID'});
    tbl_t.SubID=categorical(tbl_t.SubID);
    tbl_t.Task=categorical(tbl_t.Task);
    
    tbl_z=array2table(mat_LME2{3},'VariableNames',{'z','Task','sw','SubID'});
    tbl_z.SubID=categorical(tbl_z.SubID);
    tbl_z.Task=categorical(tbl_z.Task);
    
    [real_out, perm_out]=lme_perm(tbl_a,'sw','a~1+pred+Task+(1|SubID)',totperm);
    a_est{1}=[a_est{1} ; [nE real_out]];
    a_est{2}=[a_est{2} ; [nE*ones(totperm,1) perm_out]];
    
    [real_out, perm_out]=lme_perm(tbl_t,'sw','t~1+pred+Task+(1|SubID)',totperm);
    t_est{1}=[t_est{1} ; [nE real_out]];
    t_est{2}=[t_est{2} ; [nE*ones(totperm,1) perm_out]];
    
    [real_out, perm_out]=lme_perm(tbl_z,'sw','z~1+pred+Task+(1|SubID)',totperm);
    z_est{1}=[z_est{1} ; [nE real_out]];
    z_est{2}=[z_est{2} ; [nE*ones(totperm,1) perm_out]];
end
save('/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_HDDM_LMEandPerm_20s.mat','vGO_est','vNOGO_est','a_est','t_est','z_est')
else
    load('/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_HDDM_LMEandPerm_20s.mat')
end
%% Compute monte cartlo p
clus_alpha=0.025;
montecarlo_alpha=0.05/6;

addpath((path_fieldtrip)); ft_defaults;
cfg_neighb=[];
cfg_neighb.method = 'tri';
cfg_neighb.layout=path_PsychFTlayout;
neighbours = ft_prepare_neighbours(cfg_neighb);

[vGO_clus]=get_clusterperm_lme(vGO_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[vNOGO_clus]=get_clusterperm_lme(vNOGO_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[a_clus]=get_clusterperm_lme(a_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[t_clus]=get_clusterperm_lme(t_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[z_clus]=get_clusterperm_lme(z_est,clus_alpha,montecarlo_alpha,totperm,neighbours);


%%
cmap=cbrewer('div','RdBu',256); cmap=flipud(cmap);
Titles={'vGO','vNOGO','a','t','z'};
for nP=1:length(Titles)
    figure;
    if nP==1
        temp_topo=vGO_est{1}(:,3);
        temp_clus=vGO_clus; 
    elseif nP==2
        temp_topo=vNOGO_est{1}(:,3);
        temp_clus=vNOGO_clus;
    elseif nP==3
        temp_topo=a_est{1}(:,3);
        temp_clus=a_clus;
    elseif nP==4
        temp_topo=t_est{1}(:,3);
        temp_clus=t_clus;
    elseif nP==5
        temp_topo=z_est{1}(:,3);
        temp_clus=z_clus;
    end
    % temp_topo(match_str(layout.label,{'TP9','TP10'}))=0;
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    colorbar; 
    if nP<3
        caxis([-1 1]*6); %caxis([-1 1]*max(abs(temp_topo)))
    else
        caxis([-1 1]*12); %caxis([-1 1]*max(abs(temp_topo)))
    end
    colormap(cmap);
    
    load(path_PsychFTlayout);
    if ~isempty(temp_clus)
        for nclus=1:length(temp_clus)
        ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
        fprintf('... ... found %s cluster (%g) of %g electrodes (tval cluster=%g, Pmc=%g)\n',temp_clus{nclus}{1},nclus,length(temp_clus{nclus}{2}),temp_clus{nclus}{3},temp_clus{nclus}{4})
        end
    end
    title([Titles{nP}])
    export_fig([path_fig filesep 'LocalSleep_HDDM_LMEbyTrial_' Titles{nP} '.fig'])
    export_fig([path_fig filesep 'LocalSleep_HDDM_LMEbyTrial_' Titles{nP} '.eps'],'-r 300')
end

%%
figure; hold on;
traj=[];
for k=1:3
    m=0.65;
    while m(end)<1 && m(end)>0
        m(end+1)=m(end)+(0.02*rand-0.01);
    end
    traj{k}=m;
    if traj{k}(end)>1
        plot(traj{k},'Color','b','LineWidth',2)
    else
        plot(traj{k},'Color','r','LineWidth',2)
    end
end
format_fig;
line([-1 -1]*1000,[0 1],'Color','k','LineWidth',1,'LineStyle','--')
line(xlim,[0 0],'Color','r','LineWidth',2)
line(xlim,[1 1],'Color','b','LineWidth',2)
line([0 0],[0 1],'Color','k','LineWidth',1,'LineStyle','--')
xlim([-1200 10000])
set(gca,'box','off','visible','off')
scatter(0,0.65,'Marker','o','MarkerFaceColor',[1 1 1]*0.7,'MarkerEdgeColor','k','SizeData',300)
% export_fig([path_fig filesep 'HDDM_ExDrift4.eps'],'-r 300')
% export_fig([path_fig filesep 'HDDM_ExDrift4.fig'])