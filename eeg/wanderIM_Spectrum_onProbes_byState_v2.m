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
    
    % I probably won't use the left and right flickers individually for now
    % but note: the left and right are showing inverse to logic (needs to
    % be checked)
    left_freq(n)=SubjectInfo.FlickerL;
    right_freq(n)=SubjectInfo.FlickerR;
    
    param=[];
    param.method='fft'; % fast fourier transform
    param.mindist=1; % we want to be able to separate peaks separated by at least 1 Hz
    these_times=D.indsample(-20):D.indsample(0)-1;
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
    
    % Split between attention state
    idx_trials=find_trials(D.conditions,'ON');
    onprobe_logSNR_ON(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_ON(n,:,:)=mean(logpow(:,:,idx_trials),3);
    
    idx_trials=find_trials(D.conditions,'MW');
    onprobe_logSNR_MW(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_MW(n,:,:)=mean(logpow(:,:,idx_trials),3);
    
    idx_trials=find_trials(D.conditions,'MB');
    onprobe_logSNR_MB(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_MB(n,:,:)=mean(logpow(:,:,idx_trials),3);
    
    % Split between attention state AND stimulus type - face condition
    idx_trials=find_trials(D.conditions,'FA_ON');
    onprobe_logSNR_ON_FA(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_ON_FA(n,:,:)=mean(logpow(:,:,idx_trials),3);
    
    idx_trials=find_trials(D.conditions,'FA_MW');
    onprobe_logSNR_MW_FA(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_MW_FA(n,:,:)=mean(logpow(:,:,idx_trials),3);
    
    idx_trials=find_trials(D.conditions,'FA_MB');
    onprobe_logSNR_MB_FA(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_MB_FA(n,:,:)=mean(logpow(:,:,idx_trials),3);
    
    % Split between attention state AND stimulus type - digit condition
    idx_trials=find_trials(D.conditions,'DT_ON');
    onprobe_logSNR_ON_DT(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_ON_DT(n,:,:)=mean(logpow(:,:,idx_trials),3);
    
    idx_trials=find_trials(D.conditions,'DT_MW');
    onprobe_logSNR_MW_DT(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_MW_DT(n,:,:)=mean(logpow(:,:,idx_trials),3);
    
    idx_trials=find_trials(D.conditions,'DT_MB');
    onprobe_logSNR_MB_DT(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_MB_DT(n,:,:)=mean(logpow(:,:,idx_trials),3);
end
%% Plot the LogSNR and the Log Power (currently set to Oz channel)
% across both conditions
figure; format_fig;
subplot(1,2,2)
plot(faxis,squeeze(mean(onprobe_logSNR(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','b','LineWidth',2)
xlim([1 30])
xlabel('Frequency (Hz)')
ylabel('log(SNR)')
format_fig;
title('Oz')

subplot(1,2,1); format_fig;
plot(faxis,squeeze(nanmean(onprobe_logPow(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','r','LineWidth',2)
xlim([1 30])
xlabel('Frequency (Hz)')
ylabel('log(Power)')
format_fig;
title('Oz')
%% Split the LogSNR and LogPower plots by attention state - face condition
figure;
plot(faxis,squeeze(nanmean(onprobe_logSNR_ON(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','b')
hold on;
plot(faxis,squeeze(nanmean(onprobe_logSNR_MW(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','r')
xlim([1 30])
hold on;
plot(faxis,squeeze(nanmean(onprobe_logSNR_MB(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','g')
xlim([1 30])
legend({'On-task','Mind-wandering','Mind-blanking'})
title('LogSNR by attention state - both tasks')

figure;
plot(faxis,squeeze(nanmean(onprobe_logPow_ON(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','b')
hold on;
plot(faxis,squeeze(nanmean(onprobe_logPow_MW(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','r')
xlim([1 30])
hold on;
plot(faxis,squeeze(nanmean(onprobe_logPow_MB(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','g')
xlim([1 30])
legend({'On-task','Mind-wandering','Mind-blanking'})
title('Log Power by attention state - both tasks')

%% Comparisons across conditions
myChan=match_str(D.chanlabels,'Oz');

[~,fund_freq(1)]=findclosest(faxis,6);
[~,fund_freq(2)]=findclosest(faxis,7.5);

[~,harm_freq(1)]=findclosest(faxis,6*2);
[~,harm_freq(2)]=findclosest(faxis,7.5*2);

[~,IM_freq(1)]=findclosest(faxis,6+7.5);
% [~,IM_freq(2)]=findclosest(faxis,7.5);

temp_logSNR={onprobe_logSNR_ON,onprobe_logSNR_MW,onprobe_logSNR_MB};
temp_logSNR_FA={onprobe_logSNR_ON_DT,onprobe_logSNR_MW_DT,onprobe_logSNR_MB_DT};
temp_logSNR_DT={onprobe_logSNR_ON_FA,onprobe_logSNR_MW_FA,onprobe_logSNR_MB_FA};

figure;
subplot(1,3,1); format_fig;
hold on;
for nst=1:3
    temp_bar=squeeze(nanmean(temp_logSNR{nst}(:,myChan,fund_freq),3));
    simpleBarPlot(nst,temp_bar,'k',0.9,'r',{1 squeeze(nanmean(temp_logSNR{1}(:,myChan,fund_freq),3)) 0.05},2);
end
title('Fundamental')
set(gca,'xtick',1:3,'xticklabel',{'ON','MW','MB'})

subplot(1,3,2); format_fig;
hold on;
for nst=1:3
    temp_bar=squeeze(nanmean(temp_logSNR{nst}(:,myChan,harm_freq),3));
    simpleBarPlot(nst,temp_bar,'k',0.9,'r',{1 squeeze(nanmean(temp_logSNR{1}(:,myChan,harm_freq),3)) 0.05},2);
end
title('1st Harmonic')
set(gca,'xtick',1:3,'xticklabel',{'ON','MW','MB'})

subplot(1,3,3); format_fig;
hold on;
for nst=1:3
    temp_bar=squeeze(nanmean(temp_logSNR{nst}(:,myChan,IM_freq),3));
    simpleBarPlot(nst,temp_bar,'k',0.9,'r',{1 squeeze(nanmean(temp_logSNR{1}(:,myChan,IM_freq),3)) 0.05},2);
end
title('IM (f1+f2)')
set(gca,'xtick',1:3,'xticklabel',{'ON','MW','MB'})

figure;
subplot(1,3,1); format_fig;
hold on;
for nst=1:3
    temp_bar=squeeze(nanmean(temp_logSNR_FA{nst}(:,myChan,fund_freq),3));
    simpleBarPlot(nst-0.2,temp_bar,'k',0.35,'r',{0 squeeze(nanmean(temp_logSNR_FA{1}(:,myChan,fund_freq),3)) 0.05},2);
    
    temp_bar=squeeze(nanmean(temp_logSNR_DT{nst}(:,myChan,fund_freq),3));
    simpleBarPlot(nst+0.2,temp_bar,[1 1 1]*0.5,0.35,'r',{0 squeeze(nanmean(temp_logSNR_DT{1}(:,myChan,fund_freq),3)) 0.05},2);
end
title('Fundamental')
set(gca,'xtick',1:3,'xticklabel',{'ON','MW','MB'})

subplot(1,3,2); format_fig;
hold on;
for nst=1:3
    temp_bar=squeeze(nanmean(temp_logSNR_FA{nst}(:,myChan,harm_freq),3));
    simpleBarPlot(nst-0.2,temp_bar,'k',0.35,'r',{0 squeeze(nanmean(temp_logSNR_FA{1}(:,myChan,harm_freq),3)) 0.05},2);
    
    temp_bar=squeeze(nanmean(temp_logSNR_DT{nst}(:,myChan,harm_freq),3));
    simpleBarPlot(nst+0.2,temp_bar,[1 1 1]*0.5,0.35,'r',{0 squeeze(nanmean(temp_logSNR_DT{1}(:,myChan,harm_freq),3)) 0.05},2);
end
title('1st Harmonic')
set(gca,'xtick',1:3,'xticklabel',{'ON','MW','MB'})

subplot(1,3,3); format_fig;
hold on;
for nst=1:3
    temp_bar=squeeze(nanmean(temp_logSNR_FA{nst}(:,myChan,IM_freq),3));
    simpleBarPlot(nst-0.2,temp_bar,'k',0.35,'r',{0 squeeze(nanmean(temp_logSNR_FA{1}(:,myChan,IM_freq),3)) 0.05},2);
    
    temp_bar=squeeze(nanmean(temp_logSNR_DT{nst}(:,myChan,IM_freq),3));
    simpleBarPlot(nst+0.2,temp_bar,[1 1 1]*0.5,0.35,'r',{0 squeeze(nanmean(temp_logSNR_DT{1}(:,myChan,IM_freq),3)) 0.05},2);
end
title('IM (f1+f2)')
set(gca,'xtick',1:3,'xticklabel',{'ON','MW','MB'})
