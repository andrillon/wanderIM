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

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'lprobe_nfEEG_S3*.mat']);

%% loop across trials for baseline blocks
for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    D=spm_eeg_load([eeg_path filesep filename]);
    fprintf('... processing subject %s\n',D.fname)
    
    
    % load behavioural results
    SubID=D.fname;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    
    left_freq(n)=SubjectInfo.FlickerR; % CAREFUL this is inverted
    right_freq(n)=SubjectInfo.FlickerL;
    
    param=[];
    param.method='fft'; % fast fourier transform
    param.mindist=1; % we want to be able to separate peaks separated by at least 1 Hz
    these_times=D.indsample(-20):D.indsample(0)-1;
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    
    [logSNR, faxis, logpow]=get_logSNR(temp_data,D.fsample,param);
 
    param2=[];
    param2.method='welch'; % fast fourier transform
    param2.mindist=1; % we want to be able to separate peaks separated by at least 1 Hz
    param2.w_overlap=0.75*length(these_times)/2;
    param2.w_window=length(these_times)/2;
    param2.w_df=0.1;
    
    [logSNR2, faxis2, logpow2]=get_logSNR(temp_data,D.fsample,param2);

    onprobe_logSNR(n,:,:)=mean(logSNR,3);
    onprobe_logPow(n,:,:)=mean(logpow,3);
    onprobe_logSNR2(n,:,:)=mean(logSNR2,3);
    onprobe_logPow2(n,:,:)=mean(logpow2,3);
    
    % Split between Face and Digit
    idx_trials=find_trials(D.conditions,'FA');
    onprobe_logSNR_FA(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_FA(n,:,:)=mean(logpow(:,:,idx_trials),3);
    onprobe_logSNR_FA2(n,:,:)=mean(logSNR2(:,:,idx_trials),3);
    onprobe_logPow_FA2(n,:,:)=mean(logpow2(:,:,idx_trials),3);
    
    idx_trials=find_trials(D.conditions,'DT');
    onprobe_logSNR_DG(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_DG(n,:,:)=mean(logpow(:,:,idx_trials),3);
    onprobe_logSNR_DG2(n,:,:)=mean(logSNR2(:,:,idx_trials),3);
    onprobe_logPow_DG2(n,:,:)=mean(logpow2(:,:,idx_trials),3);
    
end

%%
figure;
subplot(1,2,1)
plot(faxis,squeeze(mean(onprobe_logSNR(:,match_str(D.chanlabels,'Oz'),:),1)),'b')
hold on;
plot(faxis2,squeeze(mean(onprobe_logSNR2(:,match_str(D.chanlabels,'Oz'),:),1)),'c')
xlim([1 30])

subplot(1,2,2)
plot(faxis,squeeze(nanmean(onprobe_logPow(:,match_str(D.chanlabels,'Oz'),:),1)),'r')
hold on
plot(faxis2,squeeze(nanmean(onprobe_logPow2(:,match_str(D.chanlabels,'Oz'),:),1)),'m')
xlim([1 30])

%%


subplot(1,2,1)
plot(faxis,squeeze(mean(onprobe_logSNR_FA(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','b')
hold on;
plot(faxis,squeeze(mean(onprobe_logSNR_DG(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','r')
xlim([1 30])
legend({'Face','Digit'})

subplot(1,2,2)
plot(faxis,squeeze(mean(onprobe_logPow_FA(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','b')
hold on;
plot(faxis,squeeze(mean(onprobe_logPow_DG(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','r')
xlim([1 30])
legend({'Face','Digit'})

%%
load('../BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
load(['EasyCap64_layout']);

plotFreqs={[6 7.5],[12 15],[13.5]};
figure;
for ntask=1:2
    for nplot=1:3
        
        subplot(2,3,nplot); format_fig;
        myFreq=plotFreqs{nplot}; idxF=[];
        for i=1:length(myFreq), [closestf, idxF(i)]=findclosest(faxis,myFreq(i)); end
        temp_topo=squeeze(mean(mean(onprobe_logSNR_FA(:,:,idxF),1),3));
        addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
        topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off');
        rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
        caxis([0 3]);
        cmap=colormap('hot'); cmap=flipud(cmap); colormap(cmap);
        colorbar;
        title(sprintf('FACE'));
        
        subplot(2,3,nplot+3); format_fig;
        myFreq=plotFreqs{nplot}; idxF=[];
        for i=1:length(myFreq), [closestf, idxF(i)]=findclosest(faxis,myFreq(i)); end
        temp_topo=squeeze(mean(mean(onprobe_logSNR_DG(:,:,idxF),1),3));
        addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
        topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off');
        rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
        caxis([0 3]);
        cmap=colormap('hot'); cmap=flipud(cmap); colormap(cmap);
        colorbar;
        title(sprintf('DIGIT'));
    end
end

%%
figure;
subplot(2,3,1); format_fig;

temp_topo1=[];
for n=1:length(left_freq)
    [~,idxclosest]=findclosest(faxis,left_freq(n)/2);
    temp_topo1(n,:)=squeeze(mean(onprobe_logSNR(n,:,idxclosest),1));
end
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(temp_topo1,1), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));

caxis([-2 2])
title('Left F')


subplot(2,3,2); format_fig;

temp_topo2=[];
for n=1:length(right_freq)
    [~,idxclosest]=findclosest(faxis,right_freq(n)/2);
    temp_topo2(n,:)=squeeze(mean(onprobe_logSNR(n,:,idxclosest),1));
end
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(temp_topo2,1), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-2 2])
title('Right F')

subplot(2,3,3); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(temp_topo2,1)-mean(temp_topo1,1), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-2 2])
title('Right-Left F')


subplot(2,3,4); format_fig;

temp_topo1=[];
for n=1:length(left_freq)
    [~,idxclosest]=findclosest(faxis,left_freq(n));
    temp_topo1(n,:)=squeeze(mean(onprobe_logSNR(n,:,idxclosest),1));
end
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(temp_topo1,1), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-2 2])
title('Left 2F')


subplot(2,3,5); format_fig;

temp_topo2=[];
for n=1:length(right_freq)
    [~,idxclosest]=findclosest(faxis,right_freq(n));
    temp_topo2(n,:)=squeeze(mean(onprobe_logSNR(n,:,idxclosest),1));
end
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(temp_topo2,1), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-2 2])
title('Right 2F')

subplot(2,3,6); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(temp_topo2,1)-mean(temp_topo1,1), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-2 2])
title('Right-Left 2F')


%%
figure;
format_fig;

temp_topo=zscore(squeeze(mean(mean(onprobe_logPow(:,:,faxis>8 & faxis<11 & (faxis~=9 | faxis~=10.5)),3),1)));

simpleTopoPlot2(mean(temp_topo,1), pos', labels,0,[],0,lay,[]);
% caxis([-2 2])
title('Alpha')

