%%
clear all
close all

run ../localdef_wanderIM.m
addpath(genpath(lscpTools_path))
addpath(genpath(path_RainCloudPlot))
addpath(genpath(path_export))

%% load data
filename = '/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_localsleep_amp_pup_thrE90P2P_Dec21_v5.txt';
% filename = '/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_localsleep_amp_pup_thrE90neg_Dec21_v5.txt';

delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

fclose(fileID);

res_table = table;
res_table.SubID = dataArray{:, 1};
res_table.BlockN = dataArray{:, 2};
res_table.TrialN = dataArray{:, 3};
res_table.ProbeN = dataArray{:, 4};
res_table.DistProbe = dataArray{:, 5};
res_table.Task = dataArray{:, 6};
res_table.StimCat = dataArray{:, 7};
res_table.Perf = dataArray{:, 8};
res_table.RCode = dataArray{:, 9};
res_table.RT = dataArray{:, 10};
res_table.W_Fp1 = dataArray{:, 11};
res_table.W_Fz = dataArray{:, 12};
res_table.W_F3 = dataArray{:, 13};
res_table.W_F7 = dataArray{:, 14};
res_table.W_FT9 = dataArray{:, 15};
res_table.W_FC5 = dataArray{:, 16};
res_table.W_FC1 = dataArray{:, 17};
res_table.W_C3 = dataArray{:, 18};
res_table.W_T7 = dataArray{:, 19};
res_table.W_TP9 = dataArray{:, 20};
res_table.W_CP5 = dataArray{:, 21};
res_table.W_CP1 = dataArray{:, 22};
res_table.W_Pz = dataArray{:, 23};
res_table.W_P3 = dataArray{:, 24};
res_table.W_P7 = dataArray{:, 25};
res_table.W_O1 = dataArray{:, 26};
res_table.W_Oz = dataArray{:, 27};
res_table.W_O2 = dataArray{:, 28};
res_table.W_P4 = dataArray{:, 29};
res_table.W_P8 = dataArray{:, 30};
res_table.W_TP10 = dataArray{:, 31};
res_table.W_CP6 = dataArray{:, 32};
res_table.W_CP2 = dataArray{:, 33};
res_table.W_Cz = dataArray{:, 34};
res_table.W_C4 = dataArray{:, 35};
res_table.W_T8 = dataArray{:, 36};
res_table.W_FT10 = dataArray{:, 37};
res_table.W_FC6 = dataArray{:, 38};
res_table.W_FC2 = dataArray{:, 39};
res_table.W_F4 = dataArray{:, 40};
res_table.W_F8 = dataArray{:, 41};
res_table.W_Fp2 = dataArray{:, 42};
res_table.W_AF7 = dataArray{:, 43};
res_table.W_AF3 = dataArray{:, 44};
res_table.W_F1 = dataArray{:, 45};
res_table.W_F5 = dataArray{:, 46};
res_table.W_FT7 = dataArray{:, 47};
res_table.W_FC3 = dataArray{:, 48};
res_table.W_C1 = dataArray{:, 49};
res_table.W_C5 = dataArray{:, 50};
res_table.W_TP7 = dataArray{:, 51};
res_table.W_CP3 = dataArray{:, 52};
res_table.W_P1 = dataArray{:, 53};
res_table.W_P5 = dataArray{:, 54};
res_table.W_PO7 = dataArray{:, 55};
res_table.W_PO3 = dataArray{:, 56};
res_table.W_POz = dataArray{:, 57};
res_table.W_PO4 = dataArray{:, 58};
res_table.W_PO8 = dataArray{:, 59};
res_table.W_P6 = dataArray{:, 60};
res_table.W_P2 = dataArray{:, 61};
res_table.W_CPz = dataArray{:, 62};
res_table.W_CP4 = dataArray{:, 63};
res_table.W_TP8 = dataArray{:, 64};
res_table.W_C6 = dataArray{:, 65};
res_table.W_C2 = dataArray{:, 66};
res_table.W_FC4 = dataArray{:, 67};
res_table.W_FT8 = dataArray{:, 68};
res_table.W_F6 = dataArray{:, 69};
res_table.W_AF8 = dataArray{:, 70};
res_table.W_AF4 = dataArray{:, 71};
res_table.W_F2 = dataArray{:, 72};
res_table.W_FCz = dataArray{:, 73};
res_table.A_Fp1 = dataArray{:, 74};
res_table.A_Fz = dataArray{:, 75};
res_table.A_F3 = dataArray{:, 76};
res_table.A_F7 = dataArray{:, 77};
res_table.A_FT9 = dataArray{:, 78};
res_table.A_FC5 = dataArray{:, 79};
res_table.A_FC1 = dataArray{:, 80};
res_table.A_C3 = dataArray{:, 81};
res_table.A_T7 = dataArray{:, 82};
res_table.A_TP9 = dataArray{:, 83};
res_table.A_CP5 = dataArray{:, 84};
res_table.A_CP1 = dataArray{:, 85};
res_table.A_Pz = dataArray{:, 86};
res_table.A_P3 = dataArray{:, 87};
res_table.A_P7 = dataArray{:, 88};
res_table.A_O1 = dataArray{:, 89};
res_table.A_Oz = dataArray{:, 90};
res_table.A_O2 = dataArray{:, 91};
res_table.A_P4 = dataArray{:, 92};
res_table.A_P8 = dataArray{:, 93};
res_table.A_TP10 = dataArray{:, 94};
res_table.A_CP6 = dataArray{:, 95};
res_table.A_CP2 = dataArray{:, 96};
res_table.A_Cz = dataArray{:, 97};
res_table.A_C4 = dataArray{:, 98};
res_table.A_T8 = dataArray{:, 99};
res_table.A_FT10 = dataArray{:, 100};
res_table.A_FC6 = dataArray{:, 101};
res_table.A_FC2 = dataArray{:, 102};
res_table.A_F4 = dataArray{:, 103};
res_table.A_F8 = dataArray{:, 104};
res_table.A_Fp2 = dataArray{:, 105};
res_table.A_AF7 = dataArray{:, 106};
res_table.A_AF3 = dataArray{:, 107};
res_table.A_F1 = dataArray{:, 108};
res_table.A_F5 = dataArray{:, 109};
res_table.A_FT7 = dataArray{:, 110};
res_table.A_FC3 = dataArray{:, 111};
res_table.A_C1 = dataArray{:, 112};
res_table.A_C5 = dataArray{:, 113};
res_table.A_TP7 = dataArray{:, 114};
res_table.A_CP3 = dataArray{:, 115};
res_table.A_P1 = dataArray{:, 116};
res_table.A_P5 = dataArray{:, 117};
res_table.A_PO7 = dataArray{:, 118};
res_table.A_PO3 = dataArray{:, 119};
res_table.A_POz = dataArray{:, 120};
res_table.A_PO4 = dataArray{:, 121};
res_table.A_PO8 = dataArray{:, 122};
res_table.A_P6 = dataArray{:, 123};
res_table.A_P2 = dataArray{:, 124};
res_table.A_CPz = dataArray{:, 125};
res_table.A_CP4 = dataArray{:, 126};
res_table.A_TP8 = dataArray{:, 127};
res_table.A_C6 = dataArray{:, 128};
res_table.A_C2 = dataArray{:, 129};
res_table.A_FC4 = dataArray{:, 130};
res_table.A_FT8 = dataArray{:, 131};
res_table.A_F6 = dataArray{:, 132};
res_table.A_AF8 = dataArray{:, 133};
res_table.A_AF4 = dataArray{:, 134};
res_table.A_F2 = dataArray{:, 135};
res_table.A_FCz = dataArray{:, 136};
res_table.State = dataArray{:, 137};
res_table.Vig = dataArray{:, 138};
res_table.Pup = dataArray{:, 139};
res_table.pcPup = dataArray{:, 140};
res_table.pcRT = dataArray{:, 141};
res_table.stimulus = dataArray{:, 142};
res_table.response = dataArray{:, 143};

res_table.W_all=nanmean(table2array(res_table(:,[11:73])),2);

%%
% clean 334
warning('removing 334')
res_table(res_table.SubID==334,:)=[];

% clean RTs
res_table(res_table.RT<0.3,:)=[];
% res_table(res_table.SubID==29,:)=[];
res_mat=table2array(res_table);

res_table.State=categorical(res_table.State);
res_table.SubID=categorical(res_table.SubID);
res_table.StimCat=categorical(res_table.StimCat);
res_table.Task=categorical(res_table.Task);

% res_table.W_TP9=nan(size(res_table,1),1);
% res_table.W_TP10=nan(size(res_table,1),1);
res_table(res_table.DistProbe>-20,:)=[];

%% LME models
mdl_0= fitlme(res_table(res_table.stimulus==1,:),sprintf('RT~1+(1|SubID)'));
mdl_1= fitlme(res_table(res_table.stimulus==1,:),sprintf('RT~1+Task+(1|SubID)'));
mdl_2= fitlme(res_table(res_table.stimulus==1,:),sprintf('RT~1+W_all+Task+pcPup+State+(1|SubID)'));
mdl_3= fitlme(res_table(res_table.stimulus==1,:),sprintf('RT~1+W_all*Task+(1|SubID)'));


%% No SPLIT TASK
GO_table=res_table(res_table.stimulus==1,:);
GO_RT=double(GO_table.RT);
GO_Perf=double(GO_table.Perf==0);
GO_table.Perf=double(GO_table.Perf==0);
NOGO_table=res_table(res_table.stimulus==0,:);
NOGO_RT=double(NOGO_table.RT);
NOGO_Perf=double(NOGO_table.Perf==0);
NOGO_table.Perf=double(NOGO_table.Perf==0);

rtGO_est=cell(1,2);
perfGO_est=cell(1,2);
perfNOGO_est=cell(1,2);

totperm=1000;

fprintf('E:%2.0f/63\n',0)
myS=unique(res_table.SubID);
redo=0;
if redo==1
    for nE=1:63
        fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
        
        %     rtGO_mdl2= fitlme(GO_table,sprintf('RT~1+%s+Task+(1|SubID)',GO_table.Properties.VariableNames{10+nE}));%,GO_table.Properties.VariableNames{10+nE}));
        %     perfGO_mdl2= fitlme(GO_table,sprintf('Perf~1+%s+Task+(1|SubID)',GO_table.Properties.VariableNames{10+nE}));%,GO_table.Properties.VariableNames{10+nE}));
        %     perfNOGO_mdl2= fitglme(NOGO_table,sprintf('Perf~1+Task+%s+(1|SubID)',NOGO_table.Properties.VariableNames{10+nE}));%,NOGO_table.Properties.VariableNames{10+nE}));
        %
        %     rtGO_est{1}=[rtGO_est{1} ; [nE double(rtGO_mdl2.Coefficients(3,2)) double(rtGO_mdl2.Coefficients(3,4)) double(rtGO_mdl2.Coefficients(3,6))]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
        %     perfGO_est{1}=[perfGO_est{1} ; [nE double(perfGO_mdl2.Coefficients(3,2)) double(perfGO_mdl2.Coefficients(3,4)) double(perfGO_mdl2.Coefficients(3,6))]];% double(perfGO_mdl.Coefficients(4,2)) double(perfGO_mdl.Coefficients(4,4)) double(perfGO_mdl.Coefficients(4,6))]];
        %     perfNOGO_est{1}=[perfNOGO_est{1} ; [nE double(perfNOGO_mdl2.Coefficients(3,2)) double(perfNOGO_mdl2.Coefficients(3,4)) double(perfNOGO_mdl2.Coefficients(3,6))]];% double(perfNOGO_mdl.Coefficients(4,2)) double(perfNOGO_mdl.Coefficients(4,4)) double(perfNOGO_mdl.Coefficients(4,6))]];
        
        [real_out, perm_out]=lme_perm(GO_table,GO_table.Properties.VariableNames{10+nE},'RT~1+pred+Task+(1|SubID)',totperm,0);
        rtGO_est{1}=[rtGO_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
        rtGO_est{2}=[rtGO_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
        
        [real_out, perm_out]=lme_perm(GO_table,GO_table.Properties.VariableNames{10+nE},'Perf~1+pred+Task+(1|SubID)',totperm,0);
        perfGO_est{1}=[perfGO_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
        perfGO_est{2}=[perfGO_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
        
        
        [real_out, perm_out]=lme_perm(NOGO_table,NOGO_table.Properties.VariableNames{10+nE},'Perf~1+pred+Task+(1|SubID)',totperm,0);
        perfNOGO_est{1}=[perfNOGO_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
        perfNOGO_est{2}=[perfNOGO_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    end
    
    save('/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_RT_Perf_LMEandPerm_20s.mat','perfNOGO_est','rtGO_est','perfGO_est')
else
    load('/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_RT_Perf_LMEandPerm_20s.mat');
end

%% Pupil by Effect
figure; set(gcf,'Position',[1241         151         387         663]);
format_fig;
Pup_TimeEffect_bySubj=[];
for nBlock=1:6
        Pup_TimeEffect_bySubj(nBlock,:)=grpstats(res_table.pcPup(res_table.BlockN==nBlock),...
            res_table.SubID(res_table.BlockN==nBlock));
end
line(1:6,nanmean(Pup_TimeEffect_bySubj,2),'Color',[0 0.8 0],'LineWidth',3);
for nBlock=1:6
simpleDotPlot(nBlock,Pup_TimeEffect_bySubj(nBlock,:),366,[0 0.8 0],0.65,'k','o',[],3,1);
end
xlabel('Block')
xlim([0.5 6.5])
ylim([1 5])
ylabel('Pupil Size')
set(gca,'FontSize',26)
% export_fig([path_fig filesep 'LocalSleep_compERPs_Pup_BlockEffect.fig'])
% export_fig([path_fig filesep 'LocalSleep_compERPs_Pup_BlockEffect.png'],'-r 300')

%% Vig by Effect
figure; set(gcf,'Position',[1241         151         387         663]);
format_fig;
Vig_TimeEffect_bySubj=[];
for nBlock=1:6
        Vig_TimeEffect_bySubj(nBlock,:)=grpstats(res_table.Vig(res_table.BlockN==nBlock),...
            res_table.SubID(res_table.BlockN==nBlock));
end
line(1:6,mean(Vig_TimeEffect_bySubj,2),'Color',[0 0 0.8],'LineWidth',3);
for nBlock=1:6
simpleDotPlot(nBlock,Vig_TimeEffect_bySubj(nBlock,:),366,[0 0 0.8],0.65,'k','o',[],3,1);
end
xlabel('Block')
xlim([0.5 6.5])
ylim([1 4])
ylabel('Vigilance Ratings')
set(gca,'FontSize',26)
% export_fig([path_fig filesep 'LocalSleep_compERPs_Vig_BlockEffect.fig'])
% export_fig([path_fig filesep 'LocalSleep_compERPs_Vig_BlockEffect.png'],'-r 300')
