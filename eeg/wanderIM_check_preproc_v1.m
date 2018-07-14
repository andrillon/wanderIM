%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

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
bsl_files=dir([eeg_path filesep 'basel_nfEEG_S3*.mat']);

%% loop across trials for baseline blocks
for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    D=spm_eeg_load([eeg_path filesep filename]);
    fprintf('... processing subject %s\n',D.fname)
    % data epoched for the 8 blokcs of 30s
    % block order: face_IM / face_noIM / digit / digit_noIM (repeated twice)
    % block length: 30s
    
    % load behavioural results
    SubID=D.fname;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    
    left_freq(n)=SubjectInfo.FlickerL;
        right_freq(n)=SubjectInfo.FlickerR;

    % Let's extact the power spectrum and compute the SNR over 10s window
    % 10s windows allows for 0.1Hz bandwith (resolution)
    % we can use a function that I prepare to compute the log SNR and power
    % spectrum
    % % % %    get_logSNR(data,SR,param)
    % % % %    data: EEG/ECoG data (channel x time x trial)
    % % % %    SR: sampling rate in Hz
    % % % %    param: structure with fields:
    % % % %      - method: 'fft' or 'taper'
    % % % %      - mindist: minimal expected distance between peaks (in Hz)
    % % % %      - numTaper: number of tapers if method if taper
    param=[];
    param.method='fft'; % fast fourier transform
    param.mindist=1; % we want to be able to separate peaks separated by at least 1 Hz
    for nC=1:4 % we have 4 conditions in total
        these_trials=[nC nC+4]; % 2 repetitions of the same condition
        
        for nW=1:3 % blocks are 30s long so we can have 3 windows of 10s
            these_times=D.indsample((nW-1)*10):D.indsample(nW*10);
            temp_data=D(1:63,these_times,these_trials); % D contains the data with channels * time * trials
            temp_data=temp_data; %-repmat(mean(temp_data,1),[size(temp_data,1) 1 1]); % re-reference the data to the average

            [logSNR, faxis, logpow]=get_logSNR(temp_data,D.fsample,param);
            
            baseline_logSNR(n,nC,nW,:,:)=squeeze(mean(logSNR,3));
            baseline_logPow(n,nC,nW,:,:)=squeeze(mean(logpow,3));
        end
        
        these_times=D.indsample(0):D.indsample(30);
        temp_data=D(1:63,these_times,these_trials); % D contains the data with channels * time * trials
        temp_data=temp_data; %-repmat(mean(temp_data,1),[size(temp_data,1) 1 1]); % re-reference the data to the average
        
        [logSNR, faxis2, logpow]=get_logSNR(temp_data,D.fsample,param);
        
        baseline_logSNR2(n,nC,1,:,:)=squeeze(mean(logSNR,3));
        baseline_logPow2(n,nC,1,:,:)=squeeze(mean(logpow,3));
        
        % sanity check: computer ERP on block onset
        these_times2=D.indsample(-0.2+0.3):D.indsample(1.5+0.3);
        these_timesbs=D.indsample(-0.2+0.3):D.indsample(0+0.3);
        temp=D(:,these_times2,these_trials)-repmat(mean(D(:,these_timesbs,these_trials),2),[1 length(these_times2) 1]);
        baseline_erp(n,nC,:,:)=squeeze(mean(temp,3));
    end
end

%% Power spectrum and log SNR - do we have signal?
this_ch=match_str(D.chanlabels,'Oz');
figure; 
subplot(1,3,1);
temp=squeeze(mean(mean(mean(baseline_logPow(:,:,:,this_ch,:),1),2),3));
% log Power average across all conditions for channel Oz
plot(faxis,temp,'Color','k','LineWidth',2);
xlim([1 30])
ylabel('log Power')
xlabel('Freq (Hz)')
 format_fig;
 
subplot(1,3,2); 
temp=squeeze(mean(mean(mean(baseline_logSNR(:,:,:,this_ch,:),1),2),3));
plot(faxis,temp,'Color','k','LineWidth',2);
xlim([1 30])
ylabel('log Power')
xlabel('Freq (Hz)')
 format_fig;
 
subplot(1,3,3); format_fig; hold on;
temp=squeeze(mean(mean(mean(baseline_logSNR(:,1,:,this_ch,:),1),2),3));
plot(faxis,temp,'Color','b','LineWidth',2,'LineStyle','-');
temp=squeeze(mean(mean(mean(baseline_logSNR(:,2,:,this_ch,:),1),2),3));
plot(faxis,temp,'Color','b','LineWidth',2,'LineStyle','--');
temp=squeeze(mean(mean(mean(baseline_logSNR(:,3,:,this_ch,:),1),2),3));
plot(faxis,temp,'Color','r','LineWidth',2,'LineStyle','-');
temp=squeeze(mean(mean(mean(baseline_logSNR(:,4,:,this_ch,:),1),2),3));
plot(faxis,temp,'Color','r','LineWidth',2,'LineStyle','--');
xlim([1 30])
ylabel('log Power')
xlabel('Freq (Hz)')

%%
this_ch=match_str(D.chanlabels,'Oz');
figure; 
subplot(1,3,1);
temp=squeeze(mean(mean(mean(baseline_logPow2(:,:,:,this_ch,:),1),2),3));
% log Power average across all conditions for channel Oz
plot(faxis2,temp,'Color','k','LineWidth',2);
xlim([1 30])
ylabel('log Power')
xlabel('Freq (Hz)')
 format_fig;
 
subplot(1,3,2); 
temp=squeeze(mean(mean(mean(baseline_logSNR2(:,:,:,this_ch,:),1),2),3));
plot(faxis2,temp,'Color','k','LineWidth',2);
xlim([1 30])
ylabel('log Power')
xlabel('Freq (Hz)')
 format_fig;
 
subplot(1,3,3); format_fig; hold on;
temp=squeeze(mean(mean(mean(baseline_logSNR2(:,1,:,this_ch,:)-baseline_logSNR2(:,2,:,this_ch,:),1),2),3));
plot(faxis2,temp,'Color','b','LineWidth',2,'LineStyle','-');
% temp=squeeze(mean(mean(mean(baseline_logSNR2(:,2,:,this_ch,:),1),2),3));
% plot(faxis2,temp,'Color','b','LineWidth',2,'LineStyle','--');
temp=squeeze(mean(mean(mean(baseline_logSNR2(:,3,:,this_ch,:)-baseline_logSNR2(:,4,:,this_ch,:),1),2),3));
plot(faxis2,temp,'Color','r','LineWidth',2,'LineStyle','-');
% temp=squeeze(mean(mean(mean(baseline_logSNR2(:,4,:,this_ch,:),1),2),3));
% plot(faxis2,temp,'Color','r','LineWidth',2,'LineStyle','--');
xlim([1 30])
ylabel('log Power')
xlabel('Freq (Hz)')

%% Where are the tags on the scalp ?
% retrieve the channels position
load('../BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
% pos=D.coor2D';
% labels=D.chanlabels(1:63);
figure;
subplot(3,2,1); format_fig; % left fondamental
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,left_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(mean(mean(baseline_logSNR2(n,:,:,:,this_f),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Left F')

subplot(3,2,2); format_fig; % right fondamental
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,right_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(mean(mean(baseline_logSNR2(n,:,:,:,this_f),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Right F')

subplot(3,2,3); format_fig; % left 2nd harmonic
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,left_freq(n));
    temp_topo(n,:)=squeeze(mean(mean(mean(baseline_logSNR2(n,:,:,:,this_f),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Left 2*F')


subplot(3,2,4); format_fig; % right 2nd harmonic
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,right_freq(n));
    temp_topo(n,:)=squeeze(mean(mean(mean(baseline_logSNR2(n,:,:,:,this_f),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Right 2*F')


subplot(3,2,5); format_fig; % IM (10.5 Hz)
[~,this_f]=findclosest(faxis2,10.5);
temp_topo=[];
temp_topo=squeeze(mean(mean(mean(baseline_logSNR2(:,:,:,:,this_f),1),2),3));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]*0.5)
title('IM 10.5')

subplot(3,2,6); format_fig; % IM (13.5 Hz)
[~,this_f]=findclosest(faxis2,13.5);
temp_topo=[];
temp_topo=squeeze(mean(mean(mean(baseline_logSNR2(:,:,:,:,this_f),1),2),3));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]*0.5)
title('IM 13.5')

%% is there a block-modulation?
[~,f1_idx]=findclosest(faxis,6);
[~,f2_idx]=findclosest(faxis,7.5);
[~,harm1_idx]=findclosest(faxis,6*2);
[~,harm2_idx]=findclosest(faxis,7.5*2);
[~,IM1_idx]=findclosest(faxis,10.5);
[~,IM2_idx]=findclosest(faxis,13.5);

logSNR_fond=squeeze(mean(baseline_logSNR(:,:,:,:,[f1_idx f2_idx]),5));
logSNR_2Harm=squeeze(mean(baseline_logSNR(:,:,:,:,[harm1_idx harm2_idx]),5));
logSNR_IM=squeeze(mean(baseline_logSNR(:,:,:,:,[IM1_idx IM2_idx]),5));

figure;
subplot(1,3,1); format_fig
% simpleBarPlot plots bar plot
% INPUT:
% - Pos         : x-position of the center of the bar
% - data        : vector of data to plot
% - colorBar    : color of the bar (string or RGB value or 2 RGB values
% (line and center)
% - widthBar    : width of the bar (e.g. 1)
% - colorError  : color of the error bar (e.g. 'r' for red)
% - sigFlag     : 3-elements cell
%                       0, 1 or 2: no stat, param stat or non-param stat
%                       p-value threshold (e.g. 0.05)
%                       H0 (one value or vector to compare)
% - widthLine   : witdht of the bar line (eg 1)
for nCond=1:4
    temp=squeeze(mean(logSNR_fond(:,nCond,:,this_ch),3)); % average across all windows for a specific condition and channel
    simpleBarPlot(nCond,temp,'k',0.9,'r',[],1);
end
set(gca,'XTick',1:4,'XTickLabel',{'F+IM','F','D+IM','D'})
ylabel('log SNR')
title('Fundamental')

subplot(1,3,2); format_fig
for nCond=1:4
    temp=squeeze(mean(logSNR_2Harm(:,nCond,:,this_ch),3)); % average across all windows for a specific condition and channel
    simpleBarPlot(nCond,temp,'k',0.9,'r',[],1);
end
set(gca,'XTick',1:4,'XTickLabel',{'F+IM','F','D+IM','D'})
ylabel('log SNR')
title('2nd Harmonic')

subplot(1,3,3); format_fig
for nCond=1:4
    temp=squeeze(mean(logSNR_IM(:,nCond,:,this_ch),3))./squeeze(mean(logSNR_fond(:,nCond,:,this_ch),3)); % average across all windows for a specific condition and channel
    simpleBarPlot(nCond,temp,'k',0.9,'r',[],1);
end
set(gca,'XTick',1:4,'XTickLabel',{'F+IM','F','D+IM','D'})
ylabel('log SNR (relative to fond)')
title('IM')

%% is there a block-modulation?
[~,f1_idx]=findclosest(faxis2,6);
[~,f2_idx]=findclosest(faxis2,7.5);
[~,harm1_idx]=findclosest(faxis2,6*2);
[~,harm2_idx]=findclosest(faxis2,7.5*2);
[~,IM1_idx]=findclosest(faxis2,10.5);
[~,IM2_idx]=findclosest(faxis2,13.5);

logSNR_fond=squeeze(mean(baseline_logSNR2(:,:,:,:,[f1_idx f2_idx]),5));
logSNR_2Harm=squeeze(mean(baseline_logSNR2(:,:,:,:,[harm1_idx harm2_idx]),5));
logSNR_IM=squeeze(mean(baseline_logSNR2(:,:,:,:,[IM1_idx IM2_idx]),5));

figure;
subplot(1,3,1); format_fig
for nCond=1:4
    temp=squeeze(mean(logSNR_fond(:,nCond,this_ch),3)); % average across all windows for a specific condition and channel
    simpleBarPlot(nCond,temp,'k',0.9,'r',[],1);
end
set(gca,'XTick',1:4,'XTickLabel',{'F+IM','F','D+IM','D'})
ylabel('log SNR')
title('Fundamental')

subplot(1,3,2); format_fig
for nCond=1:4
    temp=squeeze(mean(logSNR_2Harm(:,nCond,this_ch),3)); % average across all windows for a specific condition and channel
    simpleBarPlot(nCond,temp,'k',0.9,'r',[],1);
end
set(gca,'XTick',1:4,'XTickLabel',{'F+IM','F','D+IM','D'})
ylabel('log SNR')
title('2nd Harmonic')

subplot(1,3,3); format_fig
for nCond=1:4
    temp=squeeze(mean(logSNR_IM(:,nCond,this_ch),3)); % average across all windows for a specific condition and channel
    simpleBarPlot(nCond,temp,'k',0.9,'r',[],1);
end
set(gca,'XTick',1:4,'XTickLabel',{'F+IM','F','D+IM','D'})
ylabel('log SNR (relative to fond)')
title('IM')

%% time-course within block
figure;
subplot(1,3,1); format_fig; hold on;
temp=squeeze(logSNR_fond(:,1,:,this_ch));
errorbar(1:3,mean(temp,1),sem(temp,1),'Color','b','LineStyle','-');
temp=squeeze(logSNR_fond(:,2,:,this_ch));
errorbar((1:3)+0.1,mean(temp,1),sem(temp,1),'Color','b','LineStyle','--');
temp=squeeze(logSNR_fond(:,3,:,this_ch));
errorbar((1:3)+0.2,mean(temp,1),sem(temp,1),'Color','r','LineStyle','-');
temp=squeeze(logSNR_fond(:,4,:,this_ch));
errorbar((1:3)+0.3,mean(temp,1),sem(temp,1),'Color','r','LineStyle','--');
xlim([0.5 3.5])

subplot(1,3,2); format_fig; hold on;
temp=squeeze(logSNR_2Harm(:,1,:,this_ch));
errorbar(1:3,mean(temp,1),sem(temp,1),'Color','b','LineStyle','-');
temp=squeeze(logSNR_2Harm(:,2,:,this_ch));
errorbar((1:3)+0.1,mean(temp,1),sem(temp,1),'Color','b','LineStyle','--');
temp=squeeze(logSNR_2Harm(:,3,:,this_ch));
errorbar((1:3)+0.2,mean(temp,1),sem(temp,1),'Color','r','LineStyle','-');
temp=squeeze(logSNR_2Harm(:,4,:,this_ch));
errorbar((1:3)+0.3,mean(temp,1),sem(temp,1),'Color','r','LineStyle','--');
xlim([0.5 3.5])

subplot(1,3,3); format_fig; hold on;
temp=squeeze(logSNR_IM(:,1,:,this_ch));
errorbar(1:3,mean(temp,1),sem(temp,1),'Color','b','LineStyle','-');
temp=squeeze(logSNR_IM(:,2,:,this_ch));
errorbar((1:3)+0.1,mean(temp,1),sem(temp,1),'Color','b','LineStyle','--');
temp=squeeze(logSNR_IM(:,3,:,this_ch));
errorbar((1:3)+0.2,mean(temp,1),sem(temp,1),'Color','r','LineStyle','-');
temp=squeeze(logSNR_IM(:,4,:,this_ch));
errorbar((1:3)+0.3,mean(temp,1),sem(temp,1),'Color','r','LineStyle','--');
xlim([0.5 3.5])

%% ERP
figure;
% simpleTplot plots time-series
% INPUTS:
% - x: time vector
% - y: data (either vector (averaged-data) or matrix (non-averaged, will plot the mean+SEM) or cell array (will plot mean+SEM for muliple conditions))
% - newF: 1 (plot on a new window) 0 (keep the same window)
% - colorF: color of the plot. if y is vector or matrix 1 color. Or cell
% array of color
% - statsF: [0]: no stat; [1 alpha Ho]: ttest alpha; [2 montecarloalpha
% clusteralpha nperm]: cluster perm ; [3 alpha]: FDR correction
% - lineF: line type of the mean
% - transpF: transparency of the std fill (0->1)
% - jbFlag: use jbfill or not
% - sthFlag: 0 no smooth or smooth with window=sthFlag
% - errFlag: display or not the error estimate (1: yes; 2: dotted line; 0:
% no)
% - lineWidth: width of the lines plotted (normal: 1)
xTime=-0.2:1/D.fsample:1.5;
% average ERP across conditions and use only Oz
this_ch=match_str(D.chanlabels,'Oz');
my_ERP=squeeze(mean(baseline_erp(:,:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'b',0,'-',0.5,1,0,1,1);
