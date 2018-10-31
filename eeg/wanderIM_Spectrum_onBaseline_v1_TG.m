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
bsl_files=dir([eeg_path filesep 'lbasel_nfEEG_S3*.mat']);

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
%             temp_data=lapD(1:63,these_times,these_trials); % D contains the data with channels * time * trials
%             temp_data=temp_data; %-repmat(mean(temp_data,1),[size(temp_data,1) 1 1]); % re-reference the data to the average

%             [logSNR, faxis, logpow]=get_logSNR(temp_data,D.fsample,param);
            
%             baseline_logSNR(n,nC,nW,:,:)=squeeze(mean(logSNR,3));
%             baseline_logPow(n,nC,nW,:,:)=squeeze(mean(logpow,3));
        end
        
        these_times=D.indsample(0):D.indsample(30);
        temp_data=D(1:63,these_times,these_trials); % D contains the data with channels * time * trials
%         temp_data=temp_data; %-repmat(mean(temp_data,1),[size(temp_data,1) 1 1]); % re-reference the data to the average
        
        [logSNR, faxis2, logpow]=get_logSNR(temp_data,D.fsample,param);
        
        baseline_logSNR2(n,nC,1,:,:)=squeeze(mean(logSNR,3));
        baseline_logPow2(n,nC,1,:,:)=squeeze(mean(logpow,3));
        
        % sanity check: computer ERP on block onset
        these_times2=D.indsample(-0.2+0.3):D.indsample(1.5+0.3);
        these_timesbs=D.indsample(-0.2+0.3):D.indsample(0+0.3);
        temp=D(:,these_times2,these_trials)-repmat(mean(D(:,these_timesbs,these_trials),2),[1 length(these_times2) 1]);
        baseline_erp(n,nC,:,:)=squeeze(mean(temp,3));
        
    end
    
    [logSNR, faxis2, logpow]=get_logSNR(D(1:63,D.indsample(0):D.indsample(30),:),D.fsample,param);
    [logSNR3, faxis2, logpow3]=get_logSNR(mean(D(1:63,D.indsample(0):D.indsample(30),:),3),D.fsample,param);
    baseline_logSNR_byTrial(n,:,:)=squeeze(mean(logSNR,3));
    baseline_logSNR_avTrial(n,:,:)=logSNR3;
    
end
%% Power spectrum and log SNR - do we have signal?
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
temp=squeeze(mean(mean(mean(baseline_logSNR2(:,1,:,this_ch,:),1),2),3));
plot(faxis2,temp,'Color','b','LineWidth',2,'LineStyle','-');
temp=squeeze(mean(mean(mean(baseline_logSNR2(:,2,:,this_ch,:),1),2),3));
plot(faxis2,temp,'Color','b','LineWidth',2,'LineStyle','--');
temp=squeeze(mean(mean(mean(baseline_logSNR2(:,3,:,this_ch,:),1),2),3));
plot(faxis2,temp,'Color','r','LineWidth',2,'LineStyle','-');
temp=squeeze(mean(mean(mean(baseline_logSNR2(:,4,:,this_ch,:),1),2),3));
plot(faxis2,temp,'Color','r','LineWidth',2,'LineStyle','--');
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
load('BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
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
[~,f1_idx]=findclosest(faxis2,6);
[~,f2_idx]=findclosest(faxis2,7.5);
[~,harm1_idx]=findclosest(faxis2,6*2);
[~,harm2_idx]=findclosest(faxis2,7.5*2);
[~,IM1_idx]=findclosest(faxis2,11*(7.5-6));
[~,IM2_idx]=findclosest(faxis2,6+7.5);

logSNR_fond=squeeze(mean(baseline_logSNR2(:,:,:,:,[f1_idx f2_idx]),5));
logSNR_2Harm=squeeze(mean(baseline_logSNR2(:,:,:,:,[harm1_idx harm2_idx]),5));
logSNR_IM=squeeze(mean(baseline_logSNR2(:,:,:,:,[IM2_idx]),5));

figure;
typeblocks={'Fa','Fgap','Sq','inSq'};
for l=1:4
subplot(3,4,l); format_fig
temp_topo=squeeze(mean(logSNR_fond(:,l,:),1));
simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1]*2.5)
title(['f1 - ' typeblocks{l}])
colorbar;

subplot(3,4,l+4); format_fig
temp_topo=squeeze(mean(logSNR_2Harm(:,l,:),1));
simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1]*2.5)
title(['2f1 - ' typeblocks{l}])
colorbar;

subplot(3,4,l+8); format_fig
temp_topo=squeeze(mean(logSNR_IM(:,l,:),1));
simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
caxis([-1 1]*2.5)
title(['IM - ' typeblocks{l}])
colorbar;
end


figure
for k=1:3
    subplot(3,2,2*(k-1)+1); format_fig
    switch k
        case 1
            temp_topo=squeeze(mean(logSNR_fond(:,1,:)-logSNR_fond(:,2,:),1));
        case 2
            temp_topo=squeeze(mean(logSNR_2Harm(:,1,:)-logSNR_2Harm(:,2,:),1));
        case 3
            temp_topo=squeeze(mean(logSNR_IM(:,1,:)-logSNR_IM(:,2,:),1));
    end
    simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
    caxis([-1 1]*1)
    
    subplot(3,2,2*k); format_fig
    switch k
        case 1
            temp_topo=squeeze(mean(logSNR_fond(:,3,:)-logSNR_fond(:,4,:),1));
        case 2
            temp_topo=squeeze(mean(logSNR_2Harm(:,3,:)-logSNR_2Harm(:,4,:),1));
        case 3
            temp_topo=squeeze(mean(logSNR_IM(:,3,:)-logSNR_IM(:,4,:),1));
    end
    simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
    caxis([-1 1]*1)
end