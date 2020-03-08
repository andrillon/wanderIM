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
max_Freq=7;
frontalElecs=[1 32 33 60];
Fs_EEG=500;

window_probes=[-20 0]; %in seconds
step=0.5; maxW=5;
nhist=-maxW:step:maxW;

doperm=1;
totPerm=20;
%% Loop on files
nSc=0;
for n=1:length(files_EL)
    %%% Load EyeLink files
    subID=files_EL(n).name;
    bound=findstr(subID,'wanderIM_eyelink_s');
    subID=subID(length('wanderIM_eyelink_s')+(1:3));
    savename=['wanderIM_eyelink_S' subID];
    fprintf('... load %s from .mat file\n',subID)
    load([eyet_path filesep savename])
    Sacc_Start=EL_events.Sacc.start;
    Sacc_Duration=(EL_events.Sacc.end-EL_events.Sacc.start)/EL_headers.Fs;
    Sacc_Start(Sacc_Duration>1)=[];
    
    %%% Load slow-waves
    if exist([eeg_path filesep 'wanderIM_twa2_' subID '.mat'])~=0 %&& exist([eeg_path2 filesep 'wanderIM_twa4_' SubID '.mat'])~=0
        fprintf('... load local sleep detection for subject %s\n',subID)
        load([eeg_path filesep 'wanderIM_twa2_' subID]);
    else
        fprintf('... load local sleep detection subject %s DOES NOT EXIST\n',subID)
        continue;
    end
    if AmpCriterionIdx==9
        all_Waves(:,AmpCriterionIdx)=-all_Waves(:,AmpCriterionIdx);
    end
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./Fs_EEG);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>max_Freq)*100)
    all_Waves(all_freq>max_Freq,:)=[];
    
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
    all_Waves(all_freq>max_Freq,:)=[];
    
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
    
    Probe_Onsets=find(~(cellfun(@isempty,regexp(EL_events.Events.type,'^P[1-9]$'))) | ~(cellfun(@isempty,regexp(EL_events.Events.type,'^P10$'))));
    if length(Probe_Onsets)~=60
        warning('... wrong number of probe onsets (%g instead of %g)\n',length(Probe_Onsets),60);
    end
    ET_time=EL_events.Events.time;
    nSc=nSc+1;
    fprintf('%2.0f/%2.0f\n',0,63)
    for nE=1:63
        fprintf('\b\b\b\b\b\b%2.0f/%2.0f\n',nE,63)
        all_distP=[];
        all_distP2=[];
        all_distP_perm=[];
        all_distP2_perm=[];
        for nProbe=1:60
            
            ProbeOnset=ET_time(Probe_Onsets(nProbe));
            thisP_Sacc=(Sacc_Start(Sacc_Start>(ProbeOnset-20*EL_headers.Fs) & Sacc_Start<ProbeOnset)-ProbeOnset)./EL_headers.Fs;
            
            Sacc_PerProbe(nSc,nProbe)=length(thisP_Sacc);
            thisP_Waves=(slow_Waves(slow_Waves(:,2)==nProbe & slow_Waves(:,3)==nE,5)-20*500)/500;
            Waves_PerProbe(nSc,nProbe,nE)=size(thisP_Waves,1);
            if ~isempty(thisP_Sacc) && ~isempty(thisP_Waves)
                for k=1:length(thisP_Sacc)
                    temp=thisP_Waves-thisP_Sacc(k);
                    all_distP=[all_distP ; [temp(abs(temp)<maxW)]];
                end
            end
            
            thisP_Waves=(slow_Waves2(slow_Waves2(:,2)==nProbe & slow_Waves2(:,3)==nE,5)-20*500)/500;
            Waves_PerProbe2(nSc,nProbe,nE)=size(thisP_Waves,1);
            if ~isempty(thisP_Sacc) && ~isempty(thisP_Waves)
                for k=1:length(thisP_Sacc)
                    temp=thisP_Waves-thisP_Sacc(k);
                    all_distP2=[all_distP2 ; [temp(abs(temp)<maxW)]];
                end
            end
        end
        [nout]=histc(all_distP,[nhist]);
        if ~isempty(nout)
            nout(end)=[];
            PSTH_Sacc(nSc,nE,:)=nout;
        else
            PSTH_Sacc(nSc,nE,:)=0;
        end
        
        [nout]=histc(all_distP2,[nhist]);
        if ~isempty(nout)
            nout(end)=[];
            PSTH_Sacc2(nSc,nE,:)=nout;
        else
            PSTH_Sacc2(nSc,nE,:)=0;
        end
    end
    CorrSaccWaves(nSc,:)=corr(squeeze(Waves_PerProbe(nSc,:,:)),Sacc_PerProbe(nSc,:)');
    CorrSaccWaves2(nSc,:)=corr(squeeze(Waves_PerProbe2(nSc,:,:)),Sacc_PerProbe(nSc,:)');
    
    fprintf('\n')
    %%% do perm
    if doperm
        fprintf('... running permutations\n')
        fprintf('%3.0f%%\n',0)
        for nPerm=1:totPerm
            fprintf('\b\b\b\b\b%3.0f%%\n',nPerm/totPerm*100)
            slow_Waves_perm=slow_Waves;
            slow_Waves2_perm=slow_Waves2;
            slow_Waves_perm(:,5)=slow_Waves(randperm(size(slow_Waves,1)),5);
            slow_Waves2_perm(:,5)=slow_Waves2(randperm(size(slow_Waves2,1)),5);
            
            for nE=1:63
                
                all_distP_perm=[];
                all_distP2_perm=[];
                
                for nProbe=1:60
                      ProbeOnset=ET_time(Probe_Onsets(nProbe));
            thisP_Sacc=(Sacc_Start(Sacc_Start>(ProbeOnset-20*EL_headers.Fs) & Sacc_Start<ProbeOnset)-ProbeOnset)./EL_headers.Fs;
            
           thisP_Waves=(slow_Waves_perm(slow_Waves_perm(:,2)==nProbe & slow_Waves_perm(:,3)==nE,5)-20*500)/500;
                    Waves_PerProbe_perm(nSc,nProbe,nE,nPerm)=size(thisP_Waves,1);
                    if ~isempty(thisP_Sacc) && ~isempty(thisP_Waves)
                        for k=1:length(thisP_Sacc)
                            temp=thisP_Waves-thisP_Sacc(k);
                            all_distP_perm=[all_distP_perm ; [temp(abs(temp)<maxW)]];
                        end
                    end
                    
                    thisP_Waves=(slow_Waves2_perm(slow_Waves2_perm(:,2)==nProbe & slow_Waves2_perm(:,3)==nE,5)-20*500)/500;
                    Waves_PerProbe2_perm(nSc,nProbe,nE,nPerm)=size(thisP_Waves,1);
                    if ~isempty(thisP_Sacc) && ~isempty(thisP_Waves)
                        for k=1:length(thisP_Sacc)
                            temp=thisP_Waves-thisP_Sacc(k);
                            all_distP2_perm=[all_distP2_perm ; [temp(abs(temp)<maxW)]];
                        end
                    end
                end
                
                
                [nout]=histc(all_distP_perm,[nhist]);
                if ~isempty(nout)
                    nout(end)=[];
                    PSTH_Sacc_perm(nSc,nE,:,nPerm)=nout;
                else
                    PSTH_Sacc_perm(nSc,nE,:,nPerm)=0;
                end
                
                [nout]=histc(all_distP2_perm,[nhist]);
                if ~isempty(nout)
                    nout(end)=[];
                    PSTH_Sacc2_perm(nSc,nE,:,nPerm)=nout;
                else
                    PSTH_Sacc2_perm(nSc,nE,:,nPerm)=0;
                end
            end
            CorrSaccWaves_perm(nSc,:,nPerm)=corr(squeeze(Waves_PerProbe_perm(nSc,:,:,nPerm)),Sacc_PerProbe(nSc,:)');
            CorrSaccWaves2_perm(nSc,:,nPerm)=corr(squeeze(Waves_PerProbe2_perm(nSc,:,:,nPerm)),Sacc_PerProbe(nSc,:)');
        end
        fprintf('\ndone\n')
        
    end
end


%%
figure; format_fig;
[~,hb(1)]=simpleTplot(1:63,(CorrSaccWaves),0,'r',[3 0.05 0],'-',0.5,1,[],[],2); hold on;
[~,hb(2)]=simpleTplot(1:63,(CorrSaccWaves2),0,'b',[3 0.05 0],'-',0.5,1,[],[],2); hold on;
xlabel('Electrodes')
ylabel('Pearson Correlations with Sacc')
legend(hb,{'raw','ICA'})

%%
cmap=(colormap('parula'));
figure; set(gcf,'Position',[440    42   722   756]);
for nE=1:63
    subplot(1,2,1); format_fig;
    %     tempP=squeeze(mean(PSTH_Sacc_perm(:,nE,:,:),4));
    %     %     tempP=tempP-repmat(mean(tempP(:,abs(nhist(1:end-1))>2),2),1,size(tempP,2));%
    %     tempP=tempP./repmat(sum(tempP,2),1,size(tempP,2))*100;
    %     [~,hb(1)]=simpleTplot(nhist(1:end-1),tempP,0,[1 1 1]*0.5,[0 0.05 0],'-',0,1,[],0,1); hold on;
    %     subplot(1,2,2); format_fig;
    %     tempP=squeeze(mean(PSTH_Sacc2_perm(:,nE,:,:),4));
    %     %     tempP=tempP-repmat(mean(tempP(:,abs(nhist(1:end-1))>2),2),1,size(tempP,2));%
    %     tempP=tempP./repmat(sum(tempP,2),1,size(tempP,2))*100;
    %     [~,hb(2)]=simpleTplot(nhist(1:end-1),tempP,0,[1 1 1]*0.5,[0 0.05 0],'-',0,1,[],0,1); hold on;
    
    subplot(1,2,1); format_fig;
    tempP=squeeze(PSTH_Sacc(:,nE,:));
    %     tempP=tempP-repmat(mean(tempP(:,abs(nhist(1:end-1))>2),2),1,size(tempP,2));%
    tempP=tempP./repmat(nansum(tempP,2),1,size(tempP,2))*100;
    tempP0=squeeze(mean(PSTH_Sacc_perm(:,nE,:,:),4));
    %     tempP=tempP-repmat(mean(tempP(:,abs(nhist(1:end-1))>2),2),1,size(tempP,2));%
    tempP0=tempP0./repmat(nansum(tempP0,2),1,size(tempP0,2))*100;
    tempP=tempP-tempP0;
    [pV,hb(1)]=simpleTplot(nhist(1:end-1),tempP,0,cmap(nE,:),[3 0.05 0],'-',0,1,[],0,1); hold on;
    subplot(1,2,2); format_fig;
    tempP=squeeze(PSTH_Sacc2(:,nE,:));
    %     tempP=tempP-repmat(mean(tempP(:,abs(nhist(1:end-1))>2),2),1,size(tempP,2));%
    tempP=tempP./repmat(sum(tempP,2),1,size(tempP,2))*100;
    tempP0=squeeze(mean(PSTH_Sacc2_perm(:,nE,:,:),4));
    %     tempP=tempP-repmat(mean(tempP(:,abs(nhist(1:end-1))>2),2),1,size(tempP,2));%
    tempP0=tempP0./repmat(sum(tempP0,2),1,size(tempP0,2))*100;
    tempP=tempP-tempP0;
    [~,hb(2)]=simpleTplot(nhist(1:end-1),tempP,0,cmap(nE,:),[3 0.05 0],'-',0,1,[],0,1); hold on;
end
% subplot(1,2,1);
% ylim([-20 50])
% subplot(1,2,2);
% ylim([-20 50])
subplot(1,2,1);
ylim([-4 10])
xlabel('Time from Blink')
ylabel('Number of Slow Wave (real-perm)')


subplot(1,2,2);
ylim([-4 10])
xlabel('Time from Blink')
ylabel('Number of Slow Wave (real-perm)')
