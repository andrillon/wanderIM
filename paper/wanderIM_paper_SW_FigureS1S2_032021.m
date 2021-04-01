%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
% close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_ica'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_infEEG_S3*.mat']);


%%
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
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

res_table(res_table.RT<0.3,:)=[];
res_table(res_table.DistProbe<-20,:)=[];

res_mat=table2array(res_table);

uniqueSubID=unique(res_table.SubID);
res_table2=[];
for nS=1:length(uniqueSubID)
    for nP=1:60
        res_table2=[res_table2 ; mean(res_mat(res_table.SubID==uniqueSubID(nS) & res_table.ProbeN==nP,:),1)];
    end
end
res_table2=array2table(res_table2,'VariableNames',res_table.Properties.VariableNames);



%% loop across trials for baseline blocks
prticle_Thr=90; % 80 or 90 or 95
LimFrqW=[1 4]; % [1 4] or [4 10]
AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
% fixThr=75;
fixThr=[];
art_ampl=150;
max_posampl=75;
max_Freq=7;
frontalElecs=[1 32 33 60];
WavesAndVig=[];
nSc=0;

elecs_ERP=1:63;

ppties_Waves_Vig=[];
ERP_Waves=cell(1,4);

for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    SubID=filename;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    if ~ismember(SubID,GoodSubID)
        continue;
    end
    
    D=spm_eeg_load([bsl_files(n).folder filesep filename]);
    these_times=D.indsample(D.time);
    temp_data=D(1:63,:,:); % D contains the data with channels * time * trials
    temp_data=temp_data-repmat(mean(temp_data([10 21],:,:),1),[size(temp_data,1) 1 1]);
    
    % load behavioural results
    if exist([eeg_path filesep 'wanderIM_twa5_noica_bigwnd_' SubID '.mat'])~=0
        fprintf('... load local sleep detection for subject %s\n',SubID)
        
        load([eeg_path filesep 'wanderIM_twa5_noica_bigwnd_' SubID]);
        %         load([eeg_path2 filesep 'wanderIM_twa4_' SubID]);
    else
        fprintf('... load local sleep detection subject %s DOES NOT EXIST\n',SubID)
        continue;
    end
    
    nSc=nSc+1;
    GoodSudID{nSc}=SubID;
    if AmpCriterionIdx==9
        all_Waves(:,AmpCriterionIdx)=-all_Waves(:,AmpCriterionIdx);
    end
    %     hdr=ft_read_header([eeg_path filesep filename]);
    
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    
    all_Waves(all_Waves(:,7)>0,:)=[];
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./D.fsample);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,AmpCriterionIdx)>art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>max_posampl | all_Waves(:,14)>art_ampl| abs(all_Waves(:,15))>art_ampl)*100)
    all_Waves(all_freq>max_Freq | all_Waves(:,AmpCriterionIdx)>art_ampl | all_Waves(:,11)>max_posampl| all_Waves(:,14)>art_ampl| abs(all_Waves(:,15))>art_ampl,:)=[];
    
    all_Wpos=[];
    slow_Waves=[];
    vicinity_matrix2(nSc,:)=zeros(1,63);
    for nE=1:63
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./D.fsample;
        temp_freq=1./temp_len;
        temp_abs=1./temp_len;
        temp_p2p=thisE_Waves(:,AmpCriterionIdx);
        temp_p2p2=thisE_Waves(:,4);
        
        
        if ~isempty(fixThr)
            thr_Wave(nSc,nE)=fixThr;
        else
            %             thr_Wave2(n,nE)=prctile(all_Waves((temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),AmpCriterionIdx),prticle_Thr);
            thr_Wave(nSc,nE)=prctile(thisE_Waves(:,AmpCriterionIdx),prticle_Thr);
        end
        dens_Waves(nSc,nE)=sum(temp_p2p>thr_Wave(nSc,nE));
        %         dens_Waves(n,nE)=sum(temp_p2p>thr_Wave(n,nE) & temp_p2p2>75)/(hdr.nSamples/D.fsample);
        %           len_Waves(n,nE)=sum(all_Waves(all_Waves(:,3)==nE,4)>thr_Wave2(n,nE));
        upSlope_Waves(nSc,nE)=mean(thisE_Waves(temp_p2p>thr_Wave(nSc,nE),13));
        downSlope_Waves(nSc,nE)=mean(thisE_Waves(temp_p2p>thr_Wave(nSc,nE),12));
        
        for nPr=1:60
            thisE_Waves_Probes=thisE_Waves(temp_p2p>thr_Wave(nSc,nE) & thisE_Waves(:,2)==nPr & thisE_Waves(:,5)>-20*D.fsample & thisE_Waves(:,5)<0,:);
            
            amp=nanmean(thisE_Waves_Probes(:,4));
            dens=size(thisE_Waves_Probes(:,4),1)/20*60;
            uS=nanmean(thisE_Waves_Probes(:,13));
            dS=nanmean(thisE_Waves_Probes(:,12));
            
            temp_Pup=nanmean(res_table.pcPup(res_table.SubID==str2num(SubID) & res_table.ProbeN==nPr & res_table.StimCat==0));
            ppties_Waves_Vig=[ppties_Waves_Vig ; [nSc nE nPr amp dens uS dS probe_res(nPr,38) temp_Pup]];
            
            
        end
        if nE==24
            for nVig=1:4
                low_Blocks=find(probe_res(:,end)==nVig);
                sel_Waves=thisE_Waves(temp_p2p>thr_Wave(nSc,nE) & thisE_Waves(:,5)>-20*D.fsample & thisE_Waves(:,5)<0 & ismember(thisE_Waves(:,2),low_Blocks),:);
                onset_Waves=sel_Waves(:,5)/D.fsample;
                temp_ERP_Waves=[];
                for nW=1:length(onset_Waves)
                    if onset_Waves(nW)-1<-20 || onset_Waves(nW)+2>0
                        temp_ERP_Waves(nW,:)=nan(1,length((-1*D.fsample:2*D.fsample)));
                        continue;
                    end
                    temp=temp_data(nE,D.indsample(onset_Waves(nW))+(-1*D.fsample:2*D.fsample),sel_Waves(nW,2));
                    temp_ERP_Waves(nW,:)=temp;
                    temp_ERP_Waves(nW,:)=temp_ERP_Waves(nW,:)-mean(temp_ERP_Waves(nW,:));
                end
%                 if size(temp_ERP_Waves,1)<5
%                     continue;
%                 end
                ERP_Waves{nVig}=[ERP_Waves{nVig} ; (temp_ERP_Waves)];
                
            end
            
        end
    end
end

%%
ppties_Waves_Vig(:,8)=abs(ppties_Waves_Vig(:,8)-5);
%%
dens_high=[];
dens_low=[];

amp_high=[];
amp_low=[];

uS_high=[];
uS_low=[];

dS_high=[];
dS_low=[];

for nVig=1:4
    for nE=1:63
        dens_low(nE,nVig)=nanmean(ppties_Waves_Vig(ppties_Waves_Vig(:,2)==nE & ppties_Waves_Vig(:,8)==nVig,5));
        
        amp_low(nE,nVig)=nanmean(ppties_Waves_Vig(ppties_Waves_Vig(:,2)==nE & ppties_Waves_Vig(:,8)==nVig,4));
        
        uS_low(nE,nVig)=nanmean(ppties_Waves_Vig(ppties_Waves_Vig(:,2)==nE & ppties_Waves_Vig(:,8)==nVig,6));
        
        dS_low(nE,nVig)=nanmean(ppties_Waves_Vig(ppties_Waves_Vig(:,2)==nE & ppties_Waves_Vig(:,8)==nVig,7));
    end
end


%%
cmap=cbrewer('seq','YlOrRd',5);
colormap(cmap);

figure;
set(gcf,'Position',[ 61     8   957   790]);
tight_subplot(4,5);
for nVig=1:4
    subplot(4,5,nVig+1);
    temp_topo=dens_low(:,nVig);
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    % title('Density (low)')
    format_fig;
    caxis([3.5 8.5])
end
%
hb=colorbar;
set(hb,'Position',[0.92 0.82 0.015 0.05])

subplot(4,5,1); hold on;
temp_mean=[]; temp_sem=[];
for nVig=1:4
    A=ppties_Waves_Vig(ppties_Waves_Vig(:,8)==nVig,5);
    B=ppties_Waves_Vig(ppties_Waves_Vig(:,8)==nVig,1);
    C=grpstats(A,B);
    temp_mean(nVig)=nanmean(C);
    temp_sem(nVig)=sem(C);
    errorbar(nVig,temp_mean(nVig),temp_sem(nVig),'Marker','o','MarkerSize',10,'Color',cmap(nVig+1,:),'MarkerFaceColor',cmap(nVig+1,:),'LineWidth',2);
end
ylim([min(temp_mean-temp_sem) max(temp_mean+temp_sem)])
xlim([.5 4.5])
format_fig;
% ylabel('density')
set(gca,'XColor','w');

for nVig=1:4
    subplot(4,5,6+nVig);
    temp_topo=amp_low(:,nVig);
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    % title('Amplitude (low)')
    format_fig;
    caxis([10 65])
end
%
hb=colorbar;
set(hb,'Position',[0.92 0.6 0.015 0.05])

subplot(4,5,6); hold on;
for nVig=1:4
    A=ppties_Waves_Vig(ppties_Waves_Vig(:,8)==nVig,4);
    B=ppties_Waves_Vig(ppties_Waves_Vig(:,8)==nVig,1);
    C=grpstats(A,B);
    temp_mean(nVig)=nanmean(C);
    temp_sem(nVig)=sem(C);
    errorbar(nVig,temp_mean(nVig),temp_sem(nVig),'Marker','o','MarkerSize',10,'Color',cmap(nVig+1,:),'MarkerFaceColor',cmap(nVig+1,:),'LineWidth',2);
end
xlim([.5 4.5])
ylim([min(temp_mean-temp_sem) max(temp_mean+temp_sem)])
format_fig;
% ylabel('amplitude')
set(gca,'XColor','w');

for nVig=1:4
    subplot(4,5,11+nVig);
    temp_topo=dS_low(:,nVig);
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    % title('dSlope (low)')
    format_fig;
    caxis([150 620])
end
hb=colorbar;
set(hb,'Position',[0.92 0.38 0.015 0.05])

subplot(4,5,11); hold on;
for nVig=1:4
    A=ppties_Waves_Vig(ppties_Waves_Vig(:,8)==nVig,7);
    B=ppties_Waves_Vig(ppties_Waves_Vig(:,8)==nVig,1);
    C=grpstats(A,B);
    temp_mean(nVig)=nanmean(C);
    temp_sem(nVig)=sem(C);
    errorbar(nVig,temp_mean(nVig),temp_sem(nVig),'Marker','o','MarkerSize',10,'Color',cmap(nVig+1,:),'MarkerFaceColor',cmap(nVig+1,:),'LineWidth',2);
end
xlim([.5 4.5])
ylim([min(temp_mean-temp_sem) max(temp_mean+temp_sem)])
format_fig;
% ylabel('d-Slope')
set(gca,'XColor','w');

for nVig=1:4
    subplot(4,5,16+nVig);
    temp_topo=uS_low(:,nVig);
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    % title('uSlope (low)')
    format_fig;
    caxis([150 620])
end
hb=colorbar;
set(hb,'Position',[0.92 0.16 0.015 0.05])

subplot(4,5,16); hold on;
for nVig=1:4
    A=ppties_Waves_Vig(ppties_Waves_Vig(:,8)==nVig,6);
    B=ppties_Waves_Vig(ppties_Waves_Vig(:,8)==nVig,1);
    C=grpstats(A,B);
    temp_mean(nVig)=nanmean(C);
    temp_sem(nVig)=sem(C);
    errorbar(nVig,temp_mean(nVig),temp_sem(nVig),'Marker','o','MarkerSize',10,'Color',cmap(nVig+1,:),'MarkerFaceColor',cmap(nVig+1,:),'LineWidth',2);
end
xlim([.5 4.5])
ylim([min(temp_mean-temp_sem) max(temp_mean+temp_sem)])
format_fig;
% ylabel('u-Slope')
set(gca,'XTick',1:4);

% cmap=cbrewer('seq','YlOrRd',64);
% colormap(cmap);
% 
% export_fig([path_fig filesep 'LocalSleep_SWppties_linkFatigue_Topo.fig'])
% export_fig([path_fig filesep 'LocalSleep_SWppties_linkFatigue_Topo.png'],'-r 300')

%%

figure;
xTime=-1:1/D.fsample:2;
for nVig=1:4
    [~,hp(nVig)]=simpleTplot(xTime,ERP_Waves{nVig},0,cmap(nVig+1,:),0,'-',0.5,1,20,1,2);
end
xlim([-0.5 1])
format_fig;
ylabel('Amplitude (\muV)')
xlabel('Time from onset (s)')
% legend(hp,{'low sleepiness','high sleepiness'},'Location','southeast')
title('Average Slow Waves')

% export_fig([path_fig filesep 'LocalSleep_SWppties_linkFatigue_ERP.fig'])
% export_fig([path_fig filesep 'LocalSleep_SWppties_linkFatigue_ERP.png'],'-r 300')

%%
uS=unique(ppties_Waves_Vig(:,1));
av_ppties_Waves_Vig=[];
for nS=1:length(uS)
    for nPr=1:60
        av_ppties_Waves_Vig=[av_ppties_Waves_Vig ; nanmean(ppties_Waves_Vig(ppties_Waves_Vig(:,1)==uS(nS) & ppties_Waves_Vig(:,3)==nPr,:))];
    end
end
% nSc nE nPr amp dens uS dS probe_res(nPr,38)
table=array2table(double(av_ppties_Waves_Vig),'VariableNames',{'SubID','Elec','ProbeN','Amp','Dens','uS','dS','Vig','Pup'});

table(isnan(table.Vig),:)=[];
% table.SubID=categorical(table.SubID);
table.Elec=categorical(table.Elec);

mdl0= fitlme(table,sprintf('Dens~1+ProbeN+(1|SubID)'));
mdl1= fitlme(table,sprintf('Dens~1+ProbeN+Vig+(1|SubID)'));
mdl2= fitlme(table,sprintf('Dens~1+ProbeN+Pup+(1|SubID)'));

mdlb0= fitlme(table,sprintf('Amp~1+ProbeN+(1|SubID)'));
mdlb1= fitlme(table,sprintf('Amp~1+ProbeN+Vig+(1|SubID)'));
mdlb2= fitlme(table,sprintf('Amp~1+ProbeN+Pup+(1|SubID)'));

mdlc0= fitlme(table,sprintf('uS~1+ProbeN+(1|SubID)'));
mdlc1= fitlme(table,sprintf('uS~1+ProbeN+Vig+(1|SubID)'));
mdlc2= fitlme(table,sprintf('uS~1+ProbeN+Pup+(1|SubID)'));

mdld0= fitlme(table,sprintf('dS~1+ProbeN+(1|SubID)'));
mdld1= fitlme(table,sprintf('dS~1+ProbeN+Vig+(1|SubID)'));
mdld2= fitlme(table,sprintf('dS~1+ProbeN+Pup+(1|SubID)'));

% %% Add pupil
% 
% ppties_Waves_Vig(:,9)=nan;
% for k=1:size(ppties_Waves_Vig,1)
% %     ppties_Waves_Vig=[ppties_Waves_Vig ; [nSc nE nPr amp dens uS dS probe_res(nPr,38)]];
%     thisS=GoodSudID{ppties_Waves_Vig(k,1)};
%     thisP=double(ppties_Waves_Vig(k,3));
%     
%     this_line=find(ismember(res_table.SubID,str2num(thisS)) & res_table.ProbeN==thisP);
%     thisPup=abs(round(nanmean(res_table.pcPup(this_line)))-6);
%     
%     ppties_Waves_Vig(k,9)=thisPup;
% end
%%
% ppties_Waves_Vig_ori=ppties_Waves_Vig;
% uS=unique(ppties_Waves_Vig(:,1));
% for nS=1:length(uS)
%     temp_Pup=ppties_Waves_Vig(ppties_Waves_Vig(:,1)==uS(nS),9);
%     bins=prctile(temp_Pup,0:25:100);
%     temp_Pup2=nan(size(temp_Pup));
%     temp_Pup2(temp_Pup<bins(2))=1;
%     temp_Pup2(temp_Pup>=bins(2) & temp_Pup<bins(3))=2;
%     temp_Pup2(temp_Pup>=bins(3) & temp_Pup<bins(4))=3;
%     temp_Pup2(temp_Pup>=bins(4))=4;
%     
%     ppties_Waves_Vig(ppties_Waves_Vig(:,1)==uS(nS),9)=temp_Pup2;
% end
% uS=unique(ppties_Waves_Vig(:,1));
% r_dens=nan(63,length(uS));
% r_amp=nan(63,length(uS));
% r_uS=nan(63,length(uS));
% r_dS=nan(63,length(uS));
% for nE=1:63
%     for nS=1:length(uS)
%         if isempty(ppties_Waves_Vig(ppties_Waves_Vig(:,1)==uS(nS) & ppties_Waves_Vig(:,2)==nE,9))
%             continue;
%         end
%         [r_dens(nE,nS), pV]=corr(ppties_Waves_Vig(ppties_Waves_Vig(:,1)==uS(nS) & ppties_Waves_Vig(:,2)==nE,9),ppties_Waves_Vig(ppties_Waves_Vig(:,1)==uS(nS) & ppties_Waves_Vig(:,2)==nE,5),'type','Spearman','rows','pairwise');
%         [r_amp(nE,nS), pV]=corr(ppties_Waves_Vig(ppties_Waves_Vig(:,1)==uS(nS) & ppties_Waves_Vig(:,2)==nE,9),ppties_Waves_Vig(ppties_Waves_Vig(:,1)==uS(nS) & ppties_Waves_Vig(:,2)==nE,4),'type','Spearman','rows','pairwise');
%         [r_uS(nE,nS), pV]=corr(ppties_Waves_Vig(ppties_Waves_Vig(:,1)==uS(nS) & ppties_Waves_Vig(:,2)==nE,9),ppties_Waves_Vig(ppties_Waves_Vig(:,1)==uS(nS) & ppties_Waves_Vig(:,2)==nE,6),'type','Spearman','rows','pairwise');
%         [r_dS(nE,nS), pV]=corr(ppties_Waves_Vig(ppties_Waves_Vig(:,1)==uS(nS) & ppties_Waves_Vig(:,2)==nE,9),ppties_Waves_Vig(ppties_Waves_Vig(:,1)==uS(nS) & ppties_Waves_Vig(:,2)==nE,7),'type','Spearman','rows','pairwise');
%     end
% end
ppties_Waves_Vig_ori=ppties_Waves_Vig;
%%
dens_pup=[];
amp_pup=[];
uS_pup=[];
% dS_pup=[];
ppties_Waves_Vig(:,9)=round(ppties_Waves_Vig(:,9)); %round(minmax(ppties_Waves_Vig(:,9))*3+1);

for nPup=1:5
    for nE=1:63
        dens_pup(nE,nPup)=nanmean(ppties_Waves_Vig(ppties_Waves_Vig(:,2)==nE & round(ppties_Waves_Vig(:,9))==nPup,5));
        amp_pup(nE,nPup)=nanmean(ppties_Waves_Vig(ppties_Waves_Vig(:,2)==nE & round(ppties_Waves_Vig(:,9))==nPup,4));
        uS_pup(nE,nPup)=nanmean(ppties_Waves_Vig(ppties_Waves_Vig(:,2)==nE & round(ppties_Waves_Vig(:,9))==nPup,6));
        dS_pup(nE,nPup)=nanmean(ppties_Waves_Vig(ppties_Waves_Vig(:,2)==nE & round(ppties_Waves_Vig(:,9))==nPup,7));
    end
end

cmap=cbrewer('seq','YlOrRd',6);
colormap(cmap);

figure;
set(gcf,'Position',[61           8        1158         790]);
tight_subplot(4,6);
for nPup=1:5
    subplot(4,6,nPup+1);
    temp_topo=dens_pup(:,nPup);
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    % title('Density (low)')
    format_fig;
    caxis([3.1 8])
end
%
hb=colorbar;
set(hb,'Position',[0.92 0.82 0.015 0.05])

subplot(4,6,1); hold on;
C_all=[];
for nPup=1:5
    A=ppties_Waves_Vig(round(ppties_Waves_Vig(:,9))==nPup,5);
    B=ppties_Waves_Vig(round(ppties_Waves_Vig(:,9))==nPup,1);
    C=grpstats(A,B);
    C_all=[C_all C];
end
% C_all=zscore(C_all,[],2);
for nPup=1:5
    temp_mean(nPup)=nanmean(C_all(:,nPup));
    temp_sem(nPup)=sem(C_all(:,nPup));
    errorbar(nPup,temp_mean(nPup),temp_sem(nPup),'Marker','o','MarkerSize',10,'Color',cmap(nPup+1,:),'MarkerFaceColor',cmap(nPup+1,:),'LineWidth',2);
end
ylim([min(temp_mean-temp_sem) max(temp_mean+temp_sem)])
xlim([.5 5.5])
format_fig;
% ylabel('density')
set(gca,'XColor','w');

for nPup=1:5
    subplot(4,6,7+nPup);
    temp_topo=amp_pup(:,nPup);
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    % title('Amplitude (low)')
    format_fig;
    caxis([15 58])
end
%
hb=colorbar;
set(hb,'Position',[0.92 0.6 0.015 0.05])

subplot(4,6,7); hold on;
C_all=[];
for nPup=1:5
    A=ppties_Waves_Vig(round(ppties_Waves_Vig(:,9))==nPup,4);
    B=ppties_Waves_Vig(round(ppties_Waves_Vig(:,9))==nPup,1);
    C=grpstats(A,B);
    C_all=[C_all C];
end
% C_all=zscore(C_all,[],2);
for nPup=1:5
    temp_mean(nPup)=nanmean(C_all(:,nPup));
    temp_sem(nPup)=sem(C_all(:,nPup));
    errorbar(nPup,temp_mean(nPup),temp_sem(nPup),'Marker','o','MarkerSize',10,'Color',cmap(nPup+1,:),'MarkerFaceColor',cmap(nPup+1,:),'LineWidth',2);
end
xlim([.5 5.5])
ylim([min(temp_mean-temp_sem) max(temp_mean+temp_sem)])
format_fig;
% ylabel('amplitude')
set(gca,'XColor','w');

for nPup=1:5
    subplot(4,6,13+nPup);
    temp_topo=dS_pup(:,nPup);
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    % title('dSlope (low)')
    format_fig;
    caxis([160 480])
end
hb=colorbar;
set(hb,'Position',[0.92 0.38 0.015 0.05])

subplot(4,6,13); hold on;
C_all=[];
for nPup=1:5
    A=ppties_Waves_Vig(round(ppties_Waves_Vig(:,9))==nPup,7);
    B=ppties_Waves_Vig(round(ppties_Waves_Vig(:,9))==nPup,1);
    C=grpstats(A,B);
    C_all=[C_all C];
end
% C_all=zscore(C_all,[],2);
for nPup=1:5
    temp_mean(nPup)=nanmean(C_all(:,nPup));
    temp_sem(nPup)=sem(C_all(:,nPup));
    errorbar(nPup,temp_mean(nPup),temp_sem(nPup),'Marker','o','MarkerSize',10,'Color',cmap(nPup+1,:),'MarkerFaceColor',cmap(nPup+1,:),'LineWidth',2);
end
xlim([.5 5.5])
ylim([min(temp_mean-temp_sem) max(temp_mean+temp_sem)])
format_fig;
% ylabel('d-Slope')
set(gca,'XColor','w');

for nPup=1:5
    subplot(4,6,19+nPup);
    temp_topo=uS_pup(:,nPup);
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    % title('uSlope (low)')
    format_fig;
    caxis([170 560])
end
hb=colorbar;
set(hb,'Position',[0.92 0.16 0.015 0.05])

subplot(4,6,19); hold on;
C_all=[];
for nPup=1:5
    A=ppties_Waves_Vig(round(ppties_Waves_Vig(:,9))==nPup,6);
    B=ppties_Waves_Vig(round(ppties_Waves_Vig(:,9))==nPup,1);
    C=grpstats(A,B);
    C_all=[C_all C];
end
% C_all=zscore(C_all,[],2);
for nPup=1:5
    temp_mean(nPup)=nanmean(C_all(:,nPup));
    temp_sem(nPup)=sem(C_all(:,nPup));
    errorbar(nPup,temp_mean(nPup),temp_sem(nPup),'Marker','o','MarkerSize',10,'Color',cmap(nPup+1,:),'MarkerFaceColor',cmap(nPup+1,:),'LineWidth',2);
end
xlim([.5 5.5])
ylim([min(temp_mean-temp_sem) max(temp_mean+temp_sem)])
format_fig;
% ylabel('u-Slope')
set(gca,'XTick',1:5);

% export_fig([path_fig filesep 'LocalSleep_SWppties_linkPupil_Topo.fig'])
% export_fig([path_fig filesep 'LocalSleep_SWppties_linkPupil_Topo.png'],'-r 300')