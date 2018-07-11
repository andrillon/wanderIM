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
        myfreqs=5:0.1:40;
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
        
        these_times=D.indsample(5):D.indsample(25)-1;
        temp_data=D(1:63,these_times,these_trials); % D contains the data with channels * time * trials
        temp_data=temp_data; %-repmat(mean(temp_data,1),[size(temp_data,1) 1 1]); % re-reference the data to the average
        
        [logSNR, faxis2, logpow]=get_logSNR(temp_data,D.fsample,param);
        
        baseline_logSNR2(n,nC,1,:,:)=squeeze(mean(logSNR,3));
        baseline_logPow2(n,nC,1,:,:)=squeeze(mean(logpow,3));
        
%         % wavelet-decomposition
%         SR=D.fsample;
%         fprintf('%2.0f/%2.0f\n',0,63)
%         for nE=1:63
%             fprintf('\b\b\b\b\b\b%2.0f/%2.0f\n',nE,63)
%         for nTr=1:size(temp_data,3)
%             for nf=1:length(myfreqs)
%                 tpx2=(1:(1/myfreqs(nf)*20*SR))/SR-(1/myfreqs(nf)*20)/2;
%                 sinewave=cos(2*pi*myfreqs(nf)*tpx2);
%                 length_wave(nf)=1/myfreqs(nf)*20;
%                 tpx=linspace(-4, 4, length(sinewave));
%                 gausswave=exp(-(tpx.^2)/2);
%                 mywavelet=sinewave.*gausswave;
%                 temptf=abs(hilbert(conv(temp_data(nE,:,nTr),mywavelet,'same')));
%                 baseline_tf_decomp(n,nC,nE,nf,:,nTr)=temptf(1:10,end);
%             end
%         end
        end
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
% temp=squeeze(mean(mean(mean(baseline_logPow2(:,1,:,this_ch,:)-baseline_logPow2(:,2,:,this_ch,:),1),2),3));
temp=squeeze(mean(mean(mean(baseline_logPow2(:,1,:,this_ch,:),1),2),3));
plot(faxis2,temp,'Color','b','LineWidth',2,'LineStyle','-');
temp=squeeze(mean(mean(mean(baseline_logPow2(:,2,:,this_ch,:),1),2),3));
plot(faxis2,temp,'Color','b','LineWidth',2,'LineStyle','--');
% temp=squeeze(mean(mean(mean(baseline_logPow2(:,3,:,this_ch,:)-baseline_logPow2(:,4,:,this_ch,:),1),2),3));
temp=squeeze(mean(mean(mean(baseline_logPow2(:,3,:,this_ch,:),1),2),3));
plot(faxis2,temp,'Color','r','LineWidth',2,'LineStyle','-');
temp=squeeze(mean(mean(mean(baseline_logPow2(:,4,:,this_ch,:),1),2),3));
plot(faxis2,temp,'Color','r','LineWidth',2,'LineStyle','--');
xlim([1 30])
ylabel('log Power')
xlabel('Freq (Hz)')

%% Where are the tags on the scalp ?
% retrieve the channels position
load('../BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
% pos=D.coor2D';
% labels=D.chanlabels(1:63);
figure;
subplot(3,4,1); format_fig; % left fondamental
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,left_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(mean(mean(baseline_logSNR2(n,1:2,:,:,this_f),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Left F')

subplot(3,4,2); format_fig; % right fondamental
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,right_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(mean(mean(baseline_logSNR2(n,1:2,:,:,this_f),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Right F')

subplot(3,4,5); format_fig; % left 2nd harmonic
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,left_freq(n));
    temp_topo(n,:)=squeeze(mean(mean(mean(baseline_logSNR2(n,1:2,:,:,this_f),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Left 2*F')


subplot(3,4,6); format_fig; % right 2nd harmonic
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,right_freq(n));
    temp_topo(n,:)=squeeze(mean(mean(mean(baseline_logSNR2(n,1:2,:,:,this_f),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Right 2*F')


subplot(3,4,9); format_fig; % IM (10.5 Hz)
[~,this_f]=findclosest(faxis2,10.5);
temp_topo=[];
temp_topo=squeeze(mean(mean(mean(baseline_logSNR2(:,1:2,:,:,this_f),1),2),3));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]*0.5)
title('IM 10.5')

subplot(3,4,10); format_fig; % IM (13.5 Hz)
[~,this_f]=findclosest(faxis2,13.5);
temp_topo=[];
temp_topo=squeeze(mean(mean(mean(baseline_logSNR2(:,1:2,:,:,this_f),1),2),3));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]*0.5)
title('IM 13.5')

subplot(3,4,3); format_fig; % left fondamental
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,left_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(mean(mean(baseline_logSNR2(n,3:4,:,:,this_f),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Left F')

subplot(3,4,4); format_fig; % right fondamental
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,right_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(mean(mean(baseline_logSNR2(n,3:4,:,:,this_f),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Right F')

subplot(3,4,7); format_fig; % left 2nd harmonic
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,left_freq(n));
    temp_topo(n,:)=squeeze(mean(mean(mean(baseline_logSNR2(n,3:4,:,:,this_f),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Left 2*F')


subplot(3,4,8); format_fig; % right 2nd harmonic
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,right_freq(n));
    temp_topo(n,:)=squeeze(mean(mean(mean(baseline_logSNR2(n,3:4,:,:,this_f),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Right 2*F')


subplot(3,4,11); format_fig; % IM (10.5 Hz)
[~,this_f]=findclosest(faxis2,10.5);
temp_topo=[];
temp_topo=squeeze(mean(mean(mean(baseline_logSNR2(:,3:4,:,:,this_f),1),2),3));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]*0.5)
title('IM 10.5')

subplot(3,4,12); format_fig; % IM (13.5 Hz)
[~,this_f]=findclosest(faxis2,13.5);
temp_topo=[];
temp_topo=squeeze(mean(mean(mean(baseline_logSNR2(:,3:4,:,:,this_f),1),2),3));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]*0.5)
title('IM 13.5')


%% Where are the tags on the scalp ?
% retrieve the channels position
load('../BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
% pos=D.coor2D';
% labels=D.chanlabels(1:63);
figure;
subplot(3,4,1); format_fig; % left fondamental
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,left_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(mean(mean(diff(baseline_logSNR2(n,1:2,:,:,this_f),[],2),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Left F')

subplot(3,4,2); format_fig; % right fondamental
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,right_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(mean(mean(diff(baseline_logSNR2(n,1:2,:,:,this_f),[],2),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Right F')

subplot(3,4,5); format_fig; % left 2nd harmonic
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,left_freq(n));
    temp_topo(n,:)=squeeze(mean(mean(mean(diff(baseline_logSNR2(n,1:2,:,:,this_f),[],2),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Left 2*F')


subplot(3,4,6); format_fig; % right 2nd harmonic
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,right_freq(n));
    temp_topo(n,:)=squeeze(mean(mean(mean(diff(baseline_logSNR2(n,1:2,:,:,this_f),[],2),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Right 2*F')


subplot(3,4,9); format_fig; % IM (10.5 Hz)
[~,this_f]=findclosest(faxis2,10.5);
temp_topo=[];
temp_topo=squeeze(mean(mean(mean(diff(baseline_logSNR2(:,1:2,:,:,this_f),[],2),1),2),3));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]*0.5)
title('IM 10.5')

subplot(3,4,10); format_fig; % IM (13.5 Hz)
[~,this_f]=findclosest(faxis2,13.5);
temp_topo=[];
temp_topo=squeeze(mean(mean(mean(diff(baseline_logSNR2(:,1:2,:,:,this_f),[],2),1),2),3));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]*0.5)
title('IM 13.5')

subplot(3,4,3); format_fig; % left fondamental
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,left_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(mean(mean(diff(baseline_logSNR2(n,3:4,:,:,this_f),[],2),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Left F')

subplot(3,4,4); format_fig; % right fondamental
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,right_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(mean(mean(diff(baseline_logSNR2(n,3:4,:,:,this_f),[],2),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Right F')

subplot(3,4,7); format_fig; % left 2nd harmonic
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,left_freq(n));
    temp_topo(n,:)=squeeze(mean(mean(mean(diff(baseline_logSNR2(n,3:4,:,:,this_f),[],2),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Left 2*F')


subplot(3,4,8); format_fig; % right 2nd harmonic
temp_topo=[];
for n=1:length(left_freq)
    [~,this_f]=findclosest(faxis2,right_freq(n));
    temp_topo(n,:)=squeeze(mean(mean(mean(diff(baseline_logSNR2(n,3:4,:,:,this_f),[],2),1),2),3))';
end
simpleTopoPlot2(mean(temp_topo(:,1:63),1), pos', labels,0,[],0,lay,[]);
caxis([-1 1])
title('Right 2*F')


subplot(3,4,11); format_fig; % IM (10.5 Hz)
[~,this_f]=findclosest(faxis2,10.5);
temp_topo=[];
temp_topo=squeeze(mean(mean(mean(diff(baseline_logSNR2(:,3:4,:,:,this_f),[],2),1),2),3));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]*0.5)
title('IM 10.5')

subplot(3,4,12); format_fig; % IM (13.5 Hz)
[~,this_f]=findclosest(faxis2,13.5);
temp_topo=[];
temp_topo=squeeze(mean(mean(mean(diff(baseline_logSNR2(:,3:4,:,:,this_f),[],2),1),2),3));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-1 1]*0.5)
title('IM 13.5')

