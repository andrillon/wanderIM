%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Teigane's general notes:
%D.conditions - probe conditions
%D.chantype - shows where each channel is coming from - general)
%D.channels - all 66 channels with information about them
%D.fsample - sampling rate
%D.time - time within the epoch - we have from -32 seconds to 32 seconds. - this is the variable we will need to modify to get the right time window for our analysis

%D.indsample - (stands for index of the sample) is a function that will cut
%the data to the time window I need. To be modified when looking at
%baseline vs around probe (baseline cut into 3 10s blocks)
    
%% Init
clear all;
close all;
run localdef_wanderIM

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
    
    % I probably won't use the left and right flickers individually for now
    % but note: the left and right are showing inverse to logic (needs to
    % be checked)
    left_freq(n)=SubjectInfo.FlickerL;
    right_freq(n)=SubjectInfo.FlickerR;
    
    param=[];
    param.method='fft'; % fast fourier transform
    param.mindist=1; % we want to be able to separate peaks separated by at least 1 Hz
    these_times=D.indsample(-10):D.indsample(0)-1; 
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    
    %What the logsnr function does - We end up with a straight line from the LogSNR function - transforming the data from amplitude*time to power*time, then taking the log of the SNR to get the power chart which looks like an exponential line. In this step you take the log to get a straight line. Next we remove the average signal around the peak to get a straight, horizontal line. We then need to clean the data to remove the noise - we've used 1hz to remove the noise (the frequency chosen needs to be one that will not remove any peaks - i.e. is smaller than the expected peak frequencies - 1hz is good for us.
    %Log SNR function takes 3 inputs: Data, sampling rate, parameters (data structure)
    %getLogSNR function is the FFT function - I need this to get the spectrum (as seen in fig 2.b of Nao's analysis guide)
    
    [logSNR, faxis, logpow]=get_logSNR(temp_data,D.fsample,param);
    
%     %% I've added in this block to try and get the baseline split into 3 blocks of 10s each. Not working just yet... Copied from WanderIM_baselin_tag_v1.m
% 
%     for nW=1:3 % blocks are 30s long so we can have 3 windows of 10s
%             these_times=D.indsample((nW-1)*10):D.indsample(nW*10);
%             temp_data=lapD(1:63,these_times,these_trials); % D contains the data with channels * time * trials
%             temp_data=temp_data; %-repmat(mean(temp_data,1),[size(temp_data,1) 1 1]); % re-reference the data to the average
% 
%             [logSNR, faxis, logpow]=get_logSNR(temp_data,D.fsample,param);
%             
%             baseline_logSNR(n,nC,nW,:,:)=squeeze(mean(logSNR,3));
%             baseline_logPow(n,nC,nW,:,:)=squeeze(mean(logpow,3));
%     end
%         % End of added section
    
    onprobe_logSNR(n,:,:)=mean(logSNR,3);
    onprobe_logPow(n,:,:)=mean(logpow,3);
    
    % Split between Face and Digit
    idx_trials=find_trials(D.conditions,'FA');
    onprobe_logSNR_FA(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_FA(n,:,:)=mean(logpow(:,:,idx_trials),3);
    
    idx_trials=find_trials(D.conditions,'DT');
    onprobe_logSNR_DG(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_DG(n,:,:)=mean(logpow(:,:,idx_trials),3);
end

%% Plot the LogSNR and the Log Power (currently set to Oz channel)
% across both conditions
figure;
subplot(1,2,1)
plot(faxis,squeeze(mean(onprobe_logSNR(:,match_str(D.chanlabels,'Oz'),:),1)),'b')
xlim([1 30])

subplot(1,2,2)
plot(faxis,squeeze(nanmean(onprobe_logPow(:,match_str(D.chanlabels,'Oz'),:),1)),'r')
xlim([1 30])

%% Split the LogSNR and LogPower plots by face and digit task
figure;
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

%% Compute topographic plots for frequencies of interest across all channels
load('BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
myFreq=6;

figure;
subplot(2,2,1); format_fig;
[closestf, idxclosest]=findclosest(faxis,myFreq);
% This is where I specify that I want to reduce the matrix down to a vector
temp_topo=squeeze(mean(onprobe_logSNR(:,:,idxclosest),1));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]); colorbar;

myFreq=7.5;

subplot(2,2,2); format_fig;
[closestf, idxclosest]=findclosest(faxis,myFreq);
temp_topo=squeeze(mean(onprobe_logSNR(:,:,idxclosest),1));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]); colorbar;

myFreq=7.5-6;

subplot(2,2,3); format_fig;
[closestf, idxclosest]=findclosest(faxis,myFreq);
temp_topo=squeeze(mean(onprobe_logSNR(:,:,idxclosest),1));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]); colorbar;

myFreq=7.5+6;

subplot(2,2,4); format_fig;
[closestf, idxclosest]=findclosest(faxis,myFreq);
temp_topo=squeeze(mean(onprobe_logSNR(:,:,idxclosest),1));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]); colorbar;

%% Topographic plots targeting the frequency on the left or right (depending on participant
% Note: I likely won't use this but if I do need to - check results very
% carefully. The left and right appear on the seemingly opposite side
% as-is. Potentially a coding issue?
figure;
subplot(2,2,1); format_fig;

temp_topo=[];
for n=1:length(left_freq)
    [~,idxclosest]=findclosest(faxis,left_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(onprobe_logSNR(n,:,idxclosest),1));
end
simpleTopoPlot2(mean(temp_topo,1), pos', labels,0,[],0,lay,[]);
caxis([-2 2])
title('Left F')


subplot(2,2,2); format_fig;

temp_topo=[];
for n=1:length(right_freq)
    [~,idxclosest]=findclosest(faxis,right_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(onprobe_logSNR(n,:,idxclosest),1));
end
simpleTopoPlot2(mean(temp_topo,1), pos', labels,0,[],0,lay,[]);
caxis([-2 2])
title('Right F')

subplot(2,2,3); format_fig;

temp_topo=[];
for n=1:length(left_freq)
    [~,idxclosest]=findclosest(faxis,left_freq(n));
    temp_topo(n,:)=squeeze(mean(onprobe_logSNR(n,:,idxclosest),1));
end
simpleTopoPlot2(mean(temp_topo,1), pos', labels,0,[],0,lay,[]);
caxis([-2 2])
title('Left 2F')


subplot(2,2,4); format_fig;

temp_topo=[];
for n=1:length(right_freq)
    [~,idxclosest]=findclosest(faxis,right_freq(n));
    temp_topo(n,:)=squeeze(mean(onprobe_logSNR(n,:,idxclosest),1));
end
simpleTopoPlot2(mean(temp_topo,1), pos', labels,0,[],0,lay,[]);
caxis([-2 2])
title('Right 2F')