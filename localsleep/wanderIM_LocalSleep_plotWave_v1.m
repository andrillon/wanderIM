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
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);

LimFrqW=[1 8]; %formerly >4 Hz
prticle_Thr=80; %formerly 80

%% loop across trials for baseline blocks
all_Waves_byProbes=[];
all_Waves_byProbes2=[];
all_len=[];
all_Waves_byProbes3=[];
for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    D=spm_eeg_load([eeg_path filesep filename]);
    these_times=D.indsample(-20):D.indsample(0)-1;
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    temp_data=temp_data-repmat(mean(temp_data([10 21],:,:),1),[size(temp_data,1) 1 1]);
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    fprintf('... processing subject %s\n',SubID)
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    
    left_freq(n)=SubjectInfo.FlickerR; % CAREFUL this is inverted
    right_freq(n)=SubjectInfo.FlickerL;
    
    %     param=[];
    %     param.method='fft'; % fast fourier transform
    %     param.mindist=1; % we want to be able to separate peaks separated by at least 1 Hz
    %     these_times=D.indsample(-20):D.indsample(0)-1;
    %     temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    
    % all_Waves
    % nSub nProbe nE P2Pamp negX posX WaveEnd MaxNegPeak MaxPosPeak
    % PaxsPosPeakAmp MaxDownSlope MaxUpSlope
    
    load([eeg_path filesep 'wanderIM_twa2_' SubID])
    fprintf('%2.0f\n',0)
    for nE=1:63
        fprintf('\b\b\b%2.0f\n',nE)
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./500;
        temp_freq=1./temp_len;
        temp_abs=1./temp_len;
        temp_p2p=thisE_Waves(:,4);
        
        %         all_len=[all_len ; [repmat([n nE],length(temp_len),1) temp_len temp_abs]];
        %         thr_Wave1(n,nE)=prctile(all_Waves(all_Waves(:,3)==nE,4),80);
        %         thr_Wave2(n,nE)=prctile(temp_p2p(temp_freq>5),80);
        thr_Wave2(n,nE)=prctile(all_Waves((temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),4),prticle_Thr);
        num_Waves(n,nE)=sum(temp_p2p>thr_Wave2(n,nE));
        %           len_Waves(n,nE)=sum(all_Waves(all_Waves(:,3)==nE,4)>thr_Wave2(n,nE));
        upSlope_Waves(n,nE)=mean(thisE_Waves(temp_p2p>thr_Wave2(n,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),12));
        downSlope_Waves(n,nE)=mean(thisE_Waves(temp_p2p>thr_Wave2(n,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),11));
        
        ppties_Waves(n,nE,:)=mean(thisE_Waves(temp_p2p>thr_Wave2(n,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),[4 7 8 9 10 11 12]),1);
        
        onset_Waves=thisE_Waves(temp_p2p>thr_Wave2(n,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),5);
        probes_Waves=thisE_Waves(temp_p2p>thr_Wave2(n,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),2);
        temp_ERP_Waves=[];
        for nW=1:length(onset_Waves)
            if onset_Waves(nW)-1*D.fsample<1 || onset_Waves(nW)+2*D.fsample>size(temp_data,2)-1
                temp_ERP_Waves(nW,:)=nan(1,length((-1*D.fsample:2*D.fsample)));
                continue;
            end
            temp=temp_data(nE,onset_Waves(nW)+(-1*D.fsample:2*D.fsample),probes_Waves(nW));
            if max(abs(temp))>200
                temp_ERP_Waves(nW,:)=nan(1,length((-1*D.fsample:2*D.fsample)));
                continue;
            end
            temp_ERP_Waves(nW,:)=temp;
            temp_ERP_Waves(nW,:)=temp_ERP_Waves(nW,:)-mean(temp_ERP_Waves(nW,:));
        end
        ERP_Waves(n,nE,:)=nanmean(temp_ERP_Waves);
    end
    %                 num_Waves(nE,npr)=length(cell2mat(twa_results.channels(nE).maxnegpkamp));
    %             amp_Waves(nE,npr)=mean(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))+abs(cell2mat(twa_results.channels(nE).maxpospkamp)));
end

%%
figure; format_fig
simpleTplot(-1:1/D.fsample:2,squeeze((ERP_Waves(:,24,:))),0,'k',[0],'-',0.5,1,2,[],2);
xlim([-0.3 0.6]);
xlabel('Time from wave onset (s)')
ylabel('Amplitude (\muV)')
title('20% higher amplitude')

%% Topo of amplitude and density
figure;
addpath(genpath(path_eeglab));
temp_topo=mean(num_Waves);
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{[],'.','w',24,2});
rmpath(genpath(path_eeglab));
caxis([0 600])
colormap('parula')

% [4 7 8 9 10 11 12]
% P2Pamp WaveEnd MaxNegPeak MaxPosPeak
    % PaxsPosPeakAmp MaxDownSlope MaxUpSlope
figure;
addpath(genpath(path_eeglab));
temp_topo=mean(ppties_Waves(:,:,1));
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{[],'.','w',24,2});
rmpath(genpath(path_eeglab));
% caxis([0 600])
colormap('parula')