% Histogram
% Data: Period of All/Theta10/Theta20 Waves of Fz(Ch2)of each n in Hz
% Bins: 0:0.1:5 Hz
% Figure: 5 by 4, horizontal order

%% Init: Load instead of save
clear all;
close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
% addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
% addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);

%% loop across trials for baseline blocks
all_amp_Waves=[];
histc_record = zeros(20,51); % Documenting purposes; rows = participants, columns = bin/counts
collect_all=[];

for n=1:length(bsl_files)
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    load([eeg_path filesep 'wanderIM_twa2_' SubID])
    
    all_Waves = double(all_Waves); 
    thr_wave = prctile(all_Waves(all_Waves(:,3)==2,4),80);
    Fz_Waves = all_Waves(all_Waves(:,3)==2,:);
    halfw_samples = Fz_Waves(Fz_Waves(:,4)>=thr_wave,6) - Fz_Waves(Fz_Waves(:,4)>=thr_wave,5);
    
    period = 1./(halfw_samples*2/500); % ^ Half-period in samples --> convert to period in Hz
    % Plotting time
    [N edges] = histcounts(period,0:0.1:5.1); % output 'edges' includes rightmost bin
    histc_record(n,:)= N;
    N=(N./(sum(N))*100); % convert to percentage
    % 20 participants, 5 by 4 subplotting
    collect_all = [collect_all ; N];
    subplot(5,4,n)
    xlabel(SubID);
    plot(edges(1:length(edges)-1),N);
    axis([0 5 0 20]);
    xlabel(SubID);
end


suplabel('Period (Hz)','x');
suplabel('Percentage of Waves (%)','y');
suplabel('Period of Theta Waves (80) of Fz (Ch2)','t');





%% Generate Plot for Meanplot
figure;
mean_all = mean(collect_all);
plot(edges(1:length(edges)-1),mean_all);
axis([0 5 0 10])

suplabel('Period (Hz)','x');
suplabel('Percentage of Waves (%)','y');
suplabel('Mean Period of Theta (80) Waves of Fz (Ch2)','t');
%% Corplot for period and p2p amplitude 
clear all;
close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
% addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
% addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);

all_amp_Waves=[];
histc_record = zeros(20,51); % Documenting purposes; rows = participants, columns = bin/counts
all_period=[];
all_p2p=[];

for n=1:length(bsl_files)
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    load([eeg_path filesep 'wanderIM_twa2_' SubID])
    
    all_Waves = double(all_Waves); 
    thr_wave = prctile(all_Waves(all_Waves(:,3)==2,4),90);
    Fz_Waves = all_Waves(all_Waves(:,3)==2,:);
    halfw_samples = Fz_Waves(Fz_Waves(:,4)>=thr_wave,6) - Fz_Waves(Fz_Waves(:,4)>=thr_wave,5);
  
    % ^ Half-period in samples --> convert to period in Hz
    period = 1./(halfw_samples*2/500);
    all_period = [all_period ; period];
    all_p2p = [all_p2p; Fz_Waves(Fz_Waves(:,4)>=thr_wave,4)];
       
end

simpleCorPlot(all_period,all_p2p,[],'pearson');
xlabel('Period (Hz)')
ylabel('Peak-to-peak Amplitude')
title('Correlation between Period of Theta(90) Waves and Peak-to-peak Amplitude in Fz(Ch2) ')



