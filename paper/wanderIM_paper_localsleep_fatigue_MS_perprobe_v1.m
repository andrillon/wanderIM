%%
clear all
close all

run ../localdef_wanderIM.m
addpath(genpath(lscpTools_path))
addpath(genpath(path_RainCloudPlot))

%% load data
% filename = '/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_localsleep_pup_perprobe_thrE_P2P_v5.txt';
filename = '/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_localsleep_pup_perprobe_thrE_P2P_v5.txt';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
res_table = table;
res_table.SubID = dataArray{:, 1};
res_table.ContProbeN = dataArray{:, 2};
res_table.BlockN = dataArray{:, 3};
res_table.ProbeN = dataArray{:, 4};
res_table.Task = dataArray{:, 5};
res_table.W_Fp1 = dataArray{:, 6};
res_table.W_Fz = dataArray{:, 7};
res_table.W_F3 = dataArray{:, 8};
res_table.W_F7 = dataArray{:, 9};
res_table.W_FT9 = dataArray{:, 10};
res_table.W_FC5 = dataArray{:, 11};
res_table.W_FC1 = dataArray{:, 12};
res_table.W_C3 = dataArray{:, 13};
res_table.W_T7 = dataArray{:, 14};
res_table.W_TP9 = dataArray{:, 15};
res_table.W_CP5 = dataArray{:, 16};
res_table.W_CP1 = dataArray{:, 17};
res_table.W_Pz = dataArray{:, 18};
res_table.W_P3 = dataArray{:, 19};
res_table.W_P7 = dataArray{:, 20};
res_table.W_O1 = dataArray{:, 21};
res_table.W_Oz = dataArray{:, 22};
res_table.W_O2 = dataArray{:, 23};
res_table.W_P4 = dataArray{:, 24};
res_table.W_P8 = dataArray{:, 25};
res_table.W_TP10 = dataArray{:, 26};
res_table.W_CP6 = dataArray{:, 27};
res_table.W_CP2 = dataArray{:, 28};
res_table.W_Cz = dataArray{:, 29};
res_table.W_C4 = dataArray{:, 30};
res_table.W_T8 = dataArray{:, 31};
res_table.W_FT10 = dataArray{:, 32};
res_table.W_FC6 = dataArray{:, 33};
res_table.W_FC2 = dataArray{:, 34};
res_table.W_F4 = dataArray{:, 35};
res_table.W_F8 = dataArray{:, 36};
res_table.W_Fp2 = dataArray{:, 37};
res_table.W_AF7 = dataArray{:, 38};
res_table.W_AF3 = dataArray{:, 39};
res_table.W_F1 = dataArray{:, 40};
res_table.W_F5 = dataArray{:, 41};
res_table.W_FT7 = dataArray{:, 42};
res_table.W_FC3 = dataArray{:, 43};
res_table.W_C1 = dataArray{:, 44};
res_table.W_C5 = dataArray{:, 45};
res_table.W_TP7 = dataArray{:, 46};
res_table.W_CP3 = dataArray{:, 47};
res_table.W_P1 = dataArray{:, 48};
res_table.W_P5 = dataArray{:, 49};
res_table.W_PO7 = dataArray{:, 50};
res_table.W_PO3 = dataArray{:, 51};
res_table.W_POz = dataArray{:, 52};
res_table.W_PO4 = dataArray{:, 53};
res_table.W_PO8 = dataArray{:, 54};
res_table.W_P6 = dataArray{:, 55};
res_table.W_P2 = dataArray{:, 56};
res_table.W_CPz = dataArray{:, 57};
res_table.W_CP4 = dataArray{:, 58};
res_table.W_TP8 = dataArray{:, 59};
res_table.W_C6 = dataArray{:, 60};
res_table.W_C2 = dataArray{:, 61};
res_table.W_FC4 = dataArray{:, 62};
res_table.W_FT8 = dataArray{:, 63};
res_table.W_F6 = dataArray{:, 64};
res_table.W_AF8 = dataArray{:, 65};
res_table.W_AF4 = dataArray{:, 66};
res_table.W_F2 = dataArray{:, 67};
res_table.W_FCz = dataArray{:, 68};
res_table.A_Fp1 = dataArray{:, 69};
res_table.A_Fz = dataArray{:, 70};
res_table.A_F3 = dataArray{:, 71};
res_table.A_F7 = dataArray{:, 72};
res_table.A_FT9 = dataArray{:, 73};
res_table.A_FC5 = dataArray{:, 74};
res_table.A_FC1 = dataArray{:, 75};
res_table.A_C3 = dataArray{:, 76};
res_table.A_T7 = dataArray{:, 77};
res_table.A_TP9 = dataArray{:, 78};
res_table.A_CP5 = dataArray{:, 79};
res_table.A_CP1 = dataArray{:, 80};
res_table.A_Pz = dataArray{:, 81};
res_table.A_P3 = dataArray{:, 82};
res_table.A_P7 = dataArray{:, 83};
res_table.A_O1 = dataArray{:, 84};
res_table.A_Oz = dataArray{:, 85};
res_table.A_O2 = dataArray{:, 86};
res_table.A_P4 = dataArray{:, 87};
res_table.A_P8 = dataArray{:, 88};
res_table.A_TP10 = dataArray{:, 89};
res_table.A_CP6 = dataArray{:, 90};
res_table.A_CP2 = dataArray{:, 91};
res_table.A_Cz = dataArray{:, 92};
res_table.A_C4 = dataArray{:, 93};
res_table.A_T8 = dataArray{:, 94};
res_table.A_FT10 = dataArray{:, 95};
res_table.A_FC6 = dataArray{:, 96};
res_table.A_FC2 = dataArray{:, 97};
res_table.A_F4 = dataArray{:, 98};
res_table.A_F8 = dataArray{:, 99};
res_table.A_Fp2 = dataArray{:, 100};
res_table.A_AF7 = dataArray{:, 101};
res_table.A_AF3 = dataArray{:, 102};
res_table.A_F1 = dataArray{:, 103};
res_table.A_F5 = dataArray{:, 104};
res_table.A_FT7 = dataArray{:, 105};
res_table.A_FC3 = dataArray{:, 106};
res_table.A_C1 = dataArray{:, 107};
res_table.A_C5 = dataArray{:, 108};
res_table.A_TP7 = dataArray{:, 109};
res_table.A_CP3 = dataArray{:, 110};
res_table.A_P1 = dataArray{:, 111};
res_table.A_P5 = dataArray{:, 112};
res_table.A_PO7 = dataArray{:, 113};
res_table.A_PO3 = dataArray{:, 114};
res_table.A_POz = dataArray{:, 115};
res_table.A_PO4 = dataArray{:, 116};
res_table.A_PO8 = dataArray{:, 117};
res_table.A_P6 = dataArray{:, 118};
res_table.A_P2 = dataArray{:, 119};
res_table.A_CPz = dataArray{:, 120};
res_table.A_CP4 = dataArray{:, 121};
res_table.A_TP8 = dataArray{:, 122};
res_table.A_C6 = dataArray{:, 123};
res_table.A_C2 = dataArray{:, 124};
res_table.A_FC4 = dataArray{:, 125};
res_table.A_FT8 = dataArray{:, 126};
res_table.A_F6 = dataArray{:, 127};
res_table.A_AF8 = dataArray{:, 128};
res_table.A_AF4 = dataArray{:, 129};
res_table.A_F2 = dataArray{:, 130};
res_table.A_FCz = dataArray{:, 131};
res_table.MA_Fp1 = dataArray{:, 132};
res_table.MA_Fz = dataArray{:, 133};
res_table.MA_F3 = dataArray{:, 134};
res_table.MA_F7 = dataArray{:, 135};
res_table.MA_FT9 = dataArray{:, 136};
res_table.MA_FC5 = dataArray{:, 137};
res_table.MA_FC1 = dataArray{:, 138};
res_table.MA_C3 = dataArray{:, 139};
res_table.MA_T7 = dataArray{:, 140};
res_table.MA_TP9 = dataArray{:, 141};
res_table.MA_CP5 = dataArray{:, 142};
res_table.MA_CP1 = dataArray{:, 143};
res_table.MA_Pz = dataArray{:, 144};
res_table.MA_P3 = dataArray{:, 145};
res_table.MA_P7 = dataArray{:, 146};
res_table.MA_O1 = dataArray{:, 147};
res_table.MA_Oz = dataArray{:, 148};
res_table.MA_O2 = dataArray{:, 149};
res_table.MA_P4 = dataArray{:, 150};
res_table.MA_P8 = dataArray{:, 151};
res_table.MA_TP10 = dataArray{:, 152};
res_table.MA_CP6 = dataArray{:, 153};
res_table.MA_CP2 = dataArray{:, 154};
res_table.MA_Cz = dataArray{:, 155};
res_table.MA_C4 = dataArray{:, 156};
res_table.MA_T8 = dataArray{:, 157};
res_table.MA_FT10 = dataArray{:, 158};
res_table.MA_FC6 = dataArray{:, 159};
res_table.MA_FC2 = dataArray{:, 160};
res_table.MA_F4 = dataArray{:, 161};
res_table.MA_F8 = dataArray{:, 162};
res_table.MA_Fp2 = dataArray{:, 163};
res_table.MA_AF7 = dataArray{:, 164};
res_table.MA_AF3 = dataArray{:, 165};
res_table.MA_F1 = dataArray{:, 166};
res_table.MA_F5 = dataArray{:, 167};
res_table.MA_FT7 = dataArray{:, 168};
res_table.MA_FC3 = dataArray{:, 169};
res_table.MA_C1 = dataArray{:, 170};
res_table.MA_C5 = dataArray{:, 171};
res_table.MA_TP7 = dataArray{:, 172};
res_table.MA_CP3 = dataArray{:, 173};
res_table.MA_P1 = dataArray{:, 174};
res_table.MA_P5 = dataArray{:, 175};
res_table.MA_PO7 = dataArray{:, 176};
res_table.MA_PO3 = dataArray{:, 177};
res_table.MA_POz = dataArray{:, 178};
res_table.MA_PO4 = dataArray{:, 179};
res_table.MA_PO8 = dataArray{:, 180};
res_table.MA_P6 = dataArray{:, 181};
res_table.MA_P2 = dataArray{:, 182};
res_table.MA_CPz = dataArray{:, 183};
res_table.MA_CP4 = dataArray{:, 184};
res_table.MA_TP8 = dataArray{:, 185};
res_table.MA_C6 = dataArray{:, 186};
res_table.MA_C2 = dataArray{:, 187};
res_table.MA_FC4 = dataArray{:, 188};
res_table.MA_FT8 = dataArray{:, 189};
res_table.MA_F6 = dataArray{:, 190};
res_table.MA_AF8 = dataArray{:, 191};
res_table.MA_AF4 = dataArray{:, 192};
res_table.MA_F2 = dataArray{:, 193};
res_table.MA_FCz = dataArray{:, 194};
res_table.State = dataArray{:, 195};
res_table.Vig = dataArray{:, 196};
res_table.Pup = dataArray{:, 197};
res_table.pcPup = dataArray{:, 198};

clearvars filename delimiter startRow formatSpec fileID dataArray ans;

res_table.W_all=nanmean(table2array(res_table(:,[6:68])),2);
res_table.A_all=nanmean(table2array(res_table(:,[69:131])),2);
res_table.MA_all=nanmean(table2array(res_table(:,[132:194])),2);

%%
% clean 334
% warning('removing 334')
% res_table(res_table.SubID==334,:)=[];
% res_table(res_table.SubID==334,:)=[];

% res_table(res_table.SubID==29,:)=[];
res_mat=table2array(res_table);

res_table.State=categorical(res_table.State);
res_table.SubID=categorical(res_table.SubID);
res_table.Task=categorical(res_table.Task);

%% LME models
mdl_0= fitlme(res_table,sprintf('W_all~1+(1|SubID)'));
mdl_1= fitlme(res_table,sprintf('W_all~1+Task+(1|SubID)'));
% mdl_2= fitlme(res_table,sprintf('W_all~1+Task+ContProbeN+(1|SubID)'));
% mdl_3= fitlme(res_table,sprintf('W_all~1+Task+ContProbeN+Vig+(1|SubID)'));
% mdl_4= fitlme(res_table,sprintf('W_all~1+Task+ContProbeN+Vig+State+(1|SubID)'));
mdl_2a= fitlme(res_table,sprintf('W_all~1+Task+State+(1|SubID)'));

mdl_2b= fitlme(res_table,sprintf('W_all~1+Task+Vig+(1|SubID)'));
mdl_2c= fitlme(res_table,sprintf('W_all~1+Task+Pup+(1|SubID)'));

% mdl_3b= fitlme(res_table,sprintf('W_all~1+Task*State+(1|SubID)'));
% mdl_3c= fitlme(res_table,sprintf('W_all~1+Task*pcPup+(1|SubID)'));

%%
mySubsID=unique(res_table.SubID);
datplot=[];
for nTask=1:2
    for nVig=1:4
        tempVal=res_table.W_all(res_table.Task==num2str(nTask) & res_table.Vig==nVig);
        tempSub=res_table.SubID(res_table.Task==num2str(nTask) & res_table.Vig==nVig);
        tempCell=[];
        for nS=1:length(mySubsID)
            tempCell(nS)=nanmean(tempVal(tempSub==mySubsID(nS)));
        end
        datplot{nVig,nTask}=tempCell(~isnan(tempCell));
    end
end

figure;
rm_raincloud(datplot, [1 0 0; 0 0 1]);
xlim([0 4])
%%

% %%
% for nTask=1:2
%     Titles={'RT','GO','NOGO'};
%     for nP=1:length(Titles)
%         figure;
%         addpath((path_fieldtrip)); ft_defaults;
%         if nP==1
%             temp_topo=rtGO_est{nTask}(:,3);
%             temp_pV=rtGO_est{nTask}(:,4); %temp_topo(temp_pV>=fdr(temp_pV,0.05))=0;
%         elseif nP==2
%             temp_topo=perfGO_est{nTask}(:,3);
%             temp_pV=perfGO_est{nTask}(:,4); %temp_topo(temp_pV>=fdr(temp_pV,0.05))=0;
%         elseif nP==3
%             temp_topo=perfNOGO_est{nTask}(:,3);
%             temp_pV=perfNOGO_est{nTask}(:,4); %temp_topo(temp_pV>=fdr(temp_pV,0.05))=0;
%         end
%         % temp_topo(match_str(layout.label,{'TP9','TP10'}))=0;
%         simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
%         colorbar; caxis([-1 1]*max(abs(temp_topo)))
%         title(Titles{nP})
%
%         load(path_PsychFTlayout);
%         if ~isempty(find(temp_pV<fdr(temp_pV,0.05)))
%             ft_plot_lay_me(layout, 'chanindx',find(temp_pV<fdr(temp_pV,0.05)),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','no')
%         end
%     end
% end
% %%
% figure;
% % subplot(1,3,1)
% addpath((path_fieldtrip)); ft_defaults;
% temp_topo=nanmean(table2array(res_table(:,[11:73])),1)';
% simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'labels',[],0,1);
% colorbar; caxis([0 1]*max(abs(temp_topo)))
% title('')
%
% figure;
% % subplot(1,3,1)
% addpath((path_fieldtrip)); ft_defaults;
% temp_topo=nanmean(table2array(res_table(:,[74:136])),1)';
% simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'labels',[],0,1);
% colorbar; caxis([-1 1]*max(abs(temp_topo)))
% title('')

%%

res_table.MWvsON=nan(size(res_table,1),1);
res_table.MWvsON(res_table.State=='2')=1;
res_table.MWvsON(res_table.State=='1')=0;

res_table.MBvsON=nan(size(res_table,1),1);
res_table.MBvsON(res_table.State=='3')=1;
res_table.MBvsON(res_table.State=='1')=0;

res_table.MWvsMB=nan(size(res_table,1),1);
res_table.MWvsMB(res_table.State=='2')=1;
res_table.MWvsMB(res_table.State=='3')=0;

MWvsON_est=[];
MBvsON_est=[];
MWvsMB_est=[];
State_est=[];
Vig_est=[];
pcPup_est=[];
fprintf('E:%2.0f/63\n',0)
myS=unique(res_table.SubID);
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    for nVar=1:3
        try
            MWvsON_mdl0= fitlme(res_table,sprintf('MWvsON~1+Task+(1|SubID)'));
            MWvsON_mdl1= fitlme(res_table,sprintf('MWvsON~1+Task+%s+(1|SubID)',res_table.Properties.VariableNames{5+63*(nVar-1)+nE}));%,GO_table.Properties.VariableNames{10+nE}));
            tempcomp=compare(MWvsON_mdl0,MWvsON_mdl1);
            MWvsON_est=[MWvsON_est ; [nVar nE double(MWvsON_mdl1.Coefficients(3,2)) double(MWvsON_mdl1.Coefficients(3,4)) double(MWvsON_mdl1.Coefficients(3,6)) tempcomp.AIC(2) tempcomp.pValue(2)]];
        catch
            MWvsON_est=[MWvsON_est ; nan(1,7)];
        end
        
        try
            MBvsON_mdl0= fitlme(res_table,sprintf('MBvsON~1+Task+(1|SubID)'));
            MBvsON_mdl1= fitlme(res_table,sprintf('MBvsON~1+Task+%s+(1|SubID)',res_table.Properties.VariableNames{5+63*(nVar-1)+nE}));%,GO_table.Properties.VariableNames{10+nE}));
            tempcomp=compare(MBvsON_mdl0,MBvsON_mdl1);
            MBvsON_est=[MBvsON_est ; [nVar nE double(MBvsON_mdl1.Coefficients(3,2)) double(MBvsON_mdl1.Coefficients(3,4)) double(MBvsON_mdl1.Coefficients(3,6)) tempcomp.AIC(2) tempcomp.pValue(2)]];
        catch
            MBvsON_est=[MBvsON_est ; nan(1,7)];
        end
        
        try
            MWvsMB_mdl0= fitlme(res_table,sprintf('MWvsMB~1+Task+(1|SubID)'));
            MWvsMB_mdl1= fitlme(res_table,sprintf('MWvsMB~1+Task+%s+(1|SubID)',res_table.Properties.VariableNames{5+63*(nVar-1)+nE}));%,GO_table.Properties.VariableNames{10+nE}));
            tempcomp=compare(MWvsMB_mdl0,MWvsMB_mdl1);
            MWvsMB_est=[MWvsMB_est ; [nVar nE double(MWvsMB_mdl1.Coefficients(3,2)) double(MWvsMB_mdl1.Coefficients(3,4)) double(MWvsMB_mdl1.Coefficients(3,6)) tempcomp.AIC(2) tempcomp.pValue(2)]];
        catch
            MWvsMB_est=[MWvsMB_est ; nan(1,7)];
        end
        
        try
            State_mdl0= fitlme(res_table,sprintf('%s~1+Task+(1|SubID)',res_table.Properties.VariableNames{5+63*(nVar-1)+nE}));
            State_mdl1= fitlme(res_table,sprintf('%s~1+Task+State+(1|SubID)',res_table.Properties.VariableNames{5+63*(nVar-1)+nE}));%,GO_table.Properties.VariableNames{10+nE}));
            tempcomp=compare(State_mdl0,State_mdl1);
            State_est=[State_est ; [nVar nE double(State_mdl1.Coefficients(3,2)) double(State_mdl1.Coefficients(3,4)) double(State_mdl1.Coefficients(3,6)) double(State_mdl1.Coefficients(4,2)) double(State_mdl1.Coefficients(4,4)) double(State_mdl1.Coefficients(4,6)) tempcomp.AIC(2) tempcomp.pValue(2)]];
        catch
            State_est=[State_est ; nan(1,10)];
        end
        
        try
            Vig_mdl0= fitlme(res_table,sprintf('%s~1+Task+(1|SubID)',res_table.Properties.VariableNames{5+63*(nVar-1)+nE}));
            Vig_mdl1= fitlme(res_table,sprintf('%s~1+Task+Vig+(1|SubID)',res_table.Properties.VariableNames{5+63*(nVar-1)+nE}));%,GO_table.Properties.VariableNames{10+nE}));
            tempcomp=compare(Vig_mdl0,Vig_mdl1);
            Vig_est=[Vig_est ; [nVar nE double(Vig_mdl1.Coefficients(3,2)) double(Vig_mdl1.Coefficients(3,4)) double(Vig_mdl1.Coefficients(3,6)) tempcomp.AIC(2) tempcomp.pValue(2)]];
        catch
            Vig_est=[Vig_est ; nan(1,7)];
        end
        
        try
            pcPup_mdl0= fitlme(res_table,sprintf('%s~1+Task+(1|SubID)',res_table.Properties.VariableNames{5+63*(nVar-1)+nE}));
            pcPup_mdl1= fitlme(res_table,sprintf('%s~1+Task+pcPup+(1|SubID)',res_table.Properties.VariableNames{5+63*(nVar-1)+nE}));%,GO_table.Properties.VariableNames{10+nE}));
            tempcomp=compare(pcPup_mdl0,pcPup_mdl1);
            pcPup_est=[pcPup_est ; [nVar nE double(pcPup_mdl1.Coefficients(3,2)) double(pcPup_mdl1.Coefficients(3,4)) double(pcPup_mdl1.Coefficients(3,6)) tempcomp.AIC(2) tempcomp.pValue(2)]];
        catch
            pcPup_est=[pcPup_est ; nan(1,7)];
        end
        
    end
end


%%
for nVar=1:3
    Titles={'MWvsON','MBvsON','MWvsMB','State2','State3'};
    for nP=4:length(Titles)
        figure;
        addpath((path_fieldtrip)); ft_defaults;
        if nP==1
            temp_topo=MWvsON_est(MWvsON_est(:,1)==nVar,4);
            temp_pV=MWvsON_est(MWvsON_est(:,1)==nVar,5); %temp_topo(temp_pV>=fdr(temp_pV,0.05))=0;
        elseif nP==2
            temp_topo=MBvsON_est(MWvsON_est(:,1)==nVar,4);
            temp_pV=MBvsON_est(MWvsON_est(:,1)==nVar,5); %temp_topo(temp_pV>=fdr(temp_pV,0.05))=0;
        elseif nP==3
            temp_topo=MWvsMB_est(MWvsON_est(:,1)==nVar,4);
            temp_pV=MWvsMB_est(MWvsON_est(:,1)==nVar,5); %temp_topo(temp_pV>=fdr(temp_pV,0.05))=0;
        elseif nP==4
            temp_topo=State_est(State_est(:,1)==nVar,4);
            temp_pV=State_est(State_est(:,1)==nVar,5); %temp_topo(temp_pV>=fdr(temp_pV,0.05))=0;
        elseif nP==5
            temp_topo=State_est(State_est(:,1)==nVar,7);
            temp_pV=State_est(State_est(:,1)==nVar,8); %temp_topo(temp_pV>=fdr(temp_pV,0.05))=0;
        end
        % temp_topo(match_str(layout.label,{'TP9','TP10'}))=0;
        simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
        colorbar; caxis([-1 1]*max(abs(temp_topo)))
        title(Titles{nP})
        
        load(path_PsychFTlayout);
        if ~isempty(find(temp_pV<fdr(temp_pV,0.05)))
            ft_plot_lay_me(layout, 'chanindx',find(temp_pV<fdr(temp_pV,0.05)),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','no')
        end
    end
end

%%
for nVar=1:3
    Titles={'Vig','pcPup'};
    for nP=1:length(Titles)
        figure;
        addpath((path_fieldtrip)); ft_defaults;
        if nP==1
            temp_topo=Vig_est(Vig_est(:,1)==nVar,4);
            temp_pV=Vig_est(Vig_est(:,1)==nVar,5); %temp_topo(temp_pV>=fdr(temp_pV,0.05))=0;
        elseif nP==2
            temp_topo=pcPup_est(pcPup_est(:,1)==nVar,4);
            temp_pV=pcPup_est(pcPup_est(:,1)==nVar,5); %temp_topo(temp_pV>=fdr(temp_pV,0.05))=0;
        end
        % temp_topo(match_str(layout.label,{'TP9','TP10'}))=0;
        simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
        colorbar; caxis([-1 1]*max(abs(temp_topo)))
        title(Titles{nP})
        
        load(path_PsychFTlayout);
        if ~isempty(find(temp_pV<fdr(temp_pV,0.05)))
            ft_plot_lay_me(layout, 'chanindx',find(temp_pV<fdr(temp_pV,0.05)),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','no')
        end
    end
end

%%
nVar=1;
figure;
for nState=1:3
    subplot(1,3,nState)
    addpath((path_fieldtrip)); ft_defaults;
    temp_topo=nanmean(table2array(res_table((res_table.State)==num2str(nState),5+63*(nVar-1)+(1:63))),1)';
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    colorbar; %caxis([0 75]);
    title('')
end

%%
nVar=1;
figure;
for nState=1:4
    subplot(1,4,nState)
    addpath((path_fieldtrip)); ft_defaults;
    temp_topo=nanmean(table2array(res_table((res_table.Vig)==(nState),5+63*(nVar-1)+(1:63))),1)';
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    colorbar; caxis([0 1]*0.75);
    title('')
end

%%
figure;
%     addpath((path_fieldtrip)); ft_defaults;
temp_topo=nanmean(table2array(res_table(:,5+63*(nVar-1)+(1:63))),1)';
simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
colorbar; %caxis([0 75]);
