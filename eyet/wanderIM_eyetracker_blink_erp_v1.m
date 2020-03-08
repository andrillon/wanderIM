%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
close all;
clear all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
eeg_path2=[root_path filesep 'preproc_ica'];
behav_path=[root_path filesep 'behav'];
files_EEG=dir([eeg_path filesep 'nfEEG_S3*.mat']);

eyet_path=[root_path filesep 'eyetracker'];
files_EL=dir([eyet_path filesep 'wanderIM_eyelink_s3*.edf']);

%% loop across trials for baseline blocks
prticle_Thr=90; % 80 or 90 or 95
LimFrqW=[1 4]; % [1 4] or [4 10]
AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
% fixThr=75;
fixThr=[];
art_ampl=200;
max_posampl=75;
max_Freq=7;
frontalElecs=[1 32 33 60];
Fs_EEG=500;

window_probes=[-20 0]; %in seconds
step=0.5; maxW=5;
nhist=-maxW:step:maxW;

doperm=0;
totPerm=20;
%% Loop on files
load([pwd filesep '..' filesep 'paper' filesep 'paper_SubID'])

nSc=0;
for n=1:length(files_EL)
    %%% Load EyeLink files
    subID=files_EL(n).name;
    bound=findstr(subID,'wanderIM_eyelink_s');
    subID=subID(length('wanderIM_eyelink_s')+(1:3));
             if ~ismember(subID,GoodSudID)
        continue;
             end
             savename=['wanderIM_eyelink_S' subID];
    fprintf('... load %s from .mat file\n',subID)
    load([eyet_path filesep savename])
    Blinks_Start=EL_events.Blinks.start;
    Blinks_Duration=(EL_events.Blinks.end-EL_events.Blinks.start)/EL_headers.Fs;
    Blinks_Start(Blinks_Duration>2)=[];

    %%% Load EEG
    D=spm_eeg_load([eeg_path2 filesep 'probe_infEEG_S' subID]);
    these_times=D.indsample(-20):D.indsample(0)-1;
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    temp_data=temp_data-repmat(mean(temp_data([10 21],:,:),1),[size(temp_data,1) 1 1]); % D contains the data with channels * time * trials

    %%% Load slow-waves
    if exist([eeg_path filesep 'wanderIM_twa5_noica_' subID '.mat'])~=0 %&& exist([eeg_path2 filesep 'wanderIM_twa4_' SubID '.mat'])~=0
        fprintf('... load local sleep detection for subject %s\n',subID)
        load([eeg_path filesep 'wanderIM_twa5_noica_' subID]);
    else
        fprintf('... load local sleep detection subject %s DOES NOT EXIST\n',subID)
        continue;
    end
    if AmpCriterionIdx==9
        all_Waves(:,AmpCriterionIdx)=-all_Waves(:,AmpCriterionIdx);
    end
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./Fs_EEG);
   fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,2)>art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>max_posampl | all_Waves(:,14)>max_posampl| all_Waves(:,15)>max_posampl)*100)
    all_Waves(all_freq>max_Freq | all_Waves(:,2)>art_ampl | all_Waves(:,11)>max_posampl| all_Waves(:,14)>max_posampl| all_Waves(:,15)>max_posampl,:)=[];
    

    
    slow_Waves=[];
    for nE=1:63
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./500;
        temp_freq=1./temp_len;
        temp_abs=1./temp_len;
        temp_p2p=thisE_Waves(:,AmpCriterionIdx);
        temp_p2p2=thisE_Waves(:,4);
        
        if ~isempty(fixThr)
            thr_Wave(n,nE)=fixThr;
        else
            thr_Wave(n,nE)=prctile(all_Waves(:,AmpCriterionIdx),prticle_Thr);
        end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(n,nE),:)];
    end
    slow_Waves_marked=slow_Waves;
    slow_Waves_marked=[slow_Waves_marked nan(size(slow_Waves_marked,1),1)];
    
    %%% Load slow-waves new ica corrected
    if exist([eeg_path2 filesep 'wanderIM_twa4_' subID '.mat'])~=0 %&& exist([eeg_path2 filesep 'wanderIM_twa4_' SubID '.mat'])~=0
        fprintf('... load local sleep detection for subject %s\n',subID)
        load([eeg_path2 filesep 'wanderIM_twa4_' subID]);
    else
        fprintf('... load local sleep detection subject %s DOES NOT EXIST\n',subID)
        continue;
    end
     if AmpCriterionIdx==9
        all_Waves(:,AmpCriterionIdx)=-all_Waves(:,AmpCriterionIdx);
     end
        all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./Fs_EEG);
 fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,2)>art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>max_posampl)*100)
    all_Waves(all_freq>max_Freq | all_Waves(:,2)>art_ampl | all_Waves(:,11)>max_posampl,:)=[];

    slow_Waves2=[];
    for nE=1:63
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./500;
        temp_freq=1./temp_len;
        temp_abs=1./temp_len;
        temp_p2p=thisE_Waves(:,AmpCriterionIdx);
        temp_p2p2=thisE_Waves(:,4);
        
        if ~isempty(fixThr)
            thr_Wave2(n,nE)=fixThr;
        else
            thr_Wave2(n,nE)=prctile(all_Waves(:,AmpCriterionIdx),prticle_Thr);
        end
        slow_Waves2=[slow_Waves2 ; thisE_Waves(temp_p2p>thr_Wave2(n,nE),:)];
    end
    slow_Waves2_marked=slow_Waves2;
    slow_Waves2_marked=[slow_Waves2_marked nan(size(slow_Waves2_marked,1),1)];
    
    Probe_Onsets=find(~(cellfun(@isempty,regexp(EL_events.Events.type,'^P[1-9]$'))) | ~(cellfun(@isempty,regexp(EL_events.Events.type,'^P10$'))));
    if length(Probe_Onsets)~=60
        warning('... wrong number of probe onsets (%g instead of %g)\n',length(Probe_Onsets),60);
    end
    ET_time=EL_events.Events.time;
    nSc=nSc+1;
    fprintf('%2.0f/%2.0f\n',0,63)
    for nE=[1 32 2 24 13 17]
        fprintf('\b\b\b\b\b\b%2.0f/%2.0f\n',nE,63)
        temp_ERP_W1=[];
        temp_ERP_W2=[];
        temp_ERP_BL=[];
        markedWaves=[];
        markedWaves2=[];
        temp_wves=[];
        for nP=1:60
            ProbeOnset=ET_time(Probe_Onsets(nP));
            thisP_Blinks=(Blinks_Start(Blinks_Start>(ProbeOnset-20*EL_headers.Fs) & Blinks_Start<ProbeOnset)-ProbeOnset)./EL_headers.Fs;
            thisP_Blinks(thisP_Blinks<=-19 | thisP_Blinks>=-1)=[];
            thisP_Waves=(slow_Waves(slow_Waves(:,2)==nP & slow_Waves(:,3)==nE,5)-20*500)/500;
            thisP_Waves(thisP_Waves<=-19 | thisP_Waves>=-1)=[];
            temp_wves=[temp_wves ; slow_Waves(slow_Waves(:,2)==nP & slow_Waves(:,3)==nE,:)];
            thisP_Waves2=(slow_Waves2(slow_Waves2(:,2)==nP & slow_Waves2(:,3)==nE,5)-20*500)/500;
            thisP_Waves2(thisP_Waves2<=-19 | thisP_Waves2>=-1)=[];
            %%% ERP to Waves (no ICA)
            markedWaves=[markedWaves ; zeros(length(thisP_Waves),1)];
            for i=1:length(thisP_Waves)
                temp_eeg=temp_data(nE,round((thisP_Waves(i)+20)*500)+(-1*D.fsample:1*D.fsample),nP);
                temp_eeg=temp_eeg-mean(temp_eeg(1:D.fsample/2));
                temp_ERP_W1=[temp_ERP_W1 ; temp_eeg];
                
                if min(abs(thisP_Waves(i)-thisP_Blinks))<1
                    markedWaves(i)=findclosest(thisP_Waves(i)-thisP_Blinks,0);
                else
                    markedWaves(i)=0;
                end
            end
            %%% ERP to Waves (ICA)
            markedWaves2=[markedWaves2 ; zeros(length(thisP_Waves2),1)];
            for i=1:length(thisP_Waves2)
                temp_eeg=temp_data(nE,round((thisP_Waves2(i)+20)*500)+(-1*D.fsample:1*D.fsample),nP);
                temp_eeg=temp_eeg-mean(temp_eeg(1:D.fsample/2));
                temp_ERP_W2=[temp_ERP_W2 ; temp_eeg];
                
                 if min(abs(thisP_Waves2(i)-thisP_Blinks))<1
                    markedWaves2(i)=findclosest(thisP_Waves2(i)-thisP_Blinks,0);
                else
                    markedWaves2(i)=0;
                end
            end
            %%% ERP to Waves (Blinks)
            for i=1:length(thisP_Blinks)
                temp_eeg=temp_data(nE,round((thisP_Blinks(i)+20)*500)+(-1*D.fsample:1*D.fsample),nP);
                temp_eeg=temp_eeg-mean(temp_eeg(1:D.fsample/2));
                temp_ERP_BL=[temp_ERP_BL ; temp_eeg];
            end
        end
        nArt_W1(nSc,nE)=mean(max(temp_ERP_W1,[],2)>250);
        nArt_W2(nSc,nE)=mean(max(temp_ERP_W2,[],2)>250);
        
        nERP_W1(nSc,nE)=sum(max(abs(temp_ERP_W1),[],2)<500);
        nERP_W2(nSc,nE)=sum(max(abs(temp_ERP_W2),[],2)<500);
        nERP_BL(nSc,nE)=sum(max(abs(temp_ERP_BL),[],2)<500);
        
        nERP_W1_withBL(nSc,nE)=sum(max(abs(temp_ERP_W1),[],2)<500 & abs(markedWaves)<1 & markedWaves~=0);
        nERP_W2_withBL(nSc,nE)=sum(max(abs(temp_ERP_W2),[],2)<500 & abs(markedWaves2)<1 & markedWaves2~=0);
        
        nERP_W1_noBL(nSc,nE)=sum(max(abs(temp_ERP_W1),[],2)<500 & markedWaves==0);
        nERP_W2_noBL(nSc,nE)=sum(max(abs(temp_ERP_W2),[],2)<500 & markedWaves2==0);
        
        if nERP_W1(nSc,nE)>0
            ERP_W1(nSc,nE,:)=mean(temp_ERP_W1(max(abs(temp_ERP_W1),[],2)<500,:),1);
        else
            ERP_W1(nSc,nE,:)=nan(1,1001);
        end
        if nERP_W2(nSc,nE)>0
            ERP_W2(nSc,nE,:)=mean(temp_ERP_W2(max(abs(temp_ERP_W2),[],2)<500,:),1);
        else
            ERP_W2(nSc,nE,:)=nan(1,1001);
        end
        if nERP_BL(nSc,nE)>0
            ERP_BL(nSc,nE,:)=mean(temp_ERP_BL(max(abs(temp_ERP_BL),[],2)<500,:),1);
        else
            ERP_BL(nSc,nE,:)=nan(1,1001);
        end
        
        if nERP_W1_withBL(nSc,nE)>0
            ERP_W1_withBL(nSc,nE,:)=mean(temp_ERP_W1(max(abs(temp_ERP_W1),[],2)<500 & abs(markedWaves)<1 & markedWaves~=0,:),1);
        else
            ERP_W1_withBL(nSc,nE,:)=nan(1,1001);
        end
        if nERP_W2_withBL(nSc,nE)>0
            ERP_W2_withBL(nSc,nE,:)=mean(temp_ERP_W2(max(abs(temp_ERP_W2),[],2)<500 & abs(markedWaves2)<1 & markedWaves2~=0,:),1);
        else
            ERP_W2_withBL(nSc,nE,:)=nan(1,1001);
        end
        
        if nERP_W1_noBL(nSc,nE)>0
            ERP_W1_noBL(nSc,nE,:)=mean(temp_ERP_W1(max(abs(temp_ERP_W1),[],2)<500 & markedWaves==0,:),1);
        else
            ERP_W1_noBL(nSc,nE,:)=nan(1,1001);
        end
        if nERP_W2_noBL(nSc,nE)>0
            ERP_W2_noBL(nSc,nE,:)=mean(temp_ERP_W2(max(abs(temp_ERP_W2),[],2)<500 & markedWaves2==0,:),1);
        else
            ERP_W2_noBL(nSc,nE,:)=nan(1,1001);
        end
    end
end

%%
figure; format_fig;
hold on
plot(squeeze(nanmean(ERP_W1(:,1,:),1)),'r')
plot(squeeze(nanmean(ERP_W2(:,1,:),1)),'r--')


plot(squeeze(nanmean(ERP_W1(:,32,:),1)),'b')
plot(squeeze(nanmean(ERP_W2(:,32,:),1)),'b--')

plot(squeeze(nanmean(ERP_W1(:,24,:),1)),'k')
plot(squeeze(nanmean(ERP_W2(:,24,:),1)),'k--')


plot(squeeze(nanmean(ERP_W1(:,2,:),1)),'m')
plot(squeeze(nanmean(ERP_W2(:,2,:),1)),'m--')

plot(squeeze(nanmean(ERP_W1(:,13,:),1)),'c')
plot(squeeze(nanmean(ERP_W2(:,13,:),1)),'c--')

% plot(squeeze(nanmean(ERP_W1(:,13,:),1)),'c')
% plot(squeeze(nanmean(ERP_W2(:,13,:),1)),'c--')