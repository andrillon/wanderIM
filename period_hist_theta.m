% Histogram
% Data: Period of All/Theta10/Theta20 Waves of Fz(Ch2)of each n in Hz
% Bins: 0:0.1:5 Hz
% Figure: 5 by 4, horizontal order

%% Half-wave period histogram; see relevant slides (latest: 46)
clear all;
% close all;
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

% loop across trials for baseline blocks
all_amp_Waves=[];
histc_record = zeros(20,50); % Documenting purposes; rows = participants, columns = bin/counts
collect_all=[];
figure;

for n=1:length(bsl_files)
    
    if n==5
        continue
    end
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    
    
    all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)
        
    
    Fz_Waves = all_Waves(all_Waves(:,3)==2,:);
    thr_wave = prctile(Fz_Waves(:,9),80);
    halfw_samples = Fz_Waves(Fz_Waves(:,9)>=thr_wave,6) - Fz_Waves(Fz_Waves(:,9)>=thr_wave,5);
%     halfw_samples = Fz_Waves(:,6)-Fz_Waves(:,5);
    
    period = 1./(halfw_samples/500); % ^ Half-period in samples --> convert to Hz
   
    [N edges] = histcounts(period,0:0.2:10); % output 'edges' includes rightmost bin
    histc_record(n,:)= N;
    N=(N./(sum(N))*100); % convert to percentage
    % 20 participants, 5 by 4 subplotting
    collect_all = [collect_all ; N];
    subplot(5,4,n)
    
    xlabel(SubID);
    plot(edges(1:length(edges)-1),N);
    axis([0 10 0 10]);
    xlabel(SubID);
end
% 
% 
[ax1 h1]= suplabel('Half-Wave Period (Hz)','x');
[ax2 h2]= suplabel('Percentage of Waves (%)','y');
[ax3 h3]= suplabel('Half-Wave Period of  Theta Waves of Fz (Ch2) per Participant','t');

set(h1,'Fontsize',24);
set(h2,'Fontsize',24);
set(h3,'Fontsize',24);



%% Generate Plot for Meanplot
figure;
mean_all = mean(collect_all);
plot(edges(1:length(edges)-1),mean_all);
axis([0 10 0 7.2])

xlabel('Half-Wave Period (Hz)');
ylabel('Percentage of Waves (%)');
title('Mean Half-Wave Period of Theta and All Waves of Fz (Ch2)');
%% Corplot for period and maxnegpkamp (slide 46)
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

% for n=1:length(bsl_files)
for n =1:20; 
    if n==5
        continue
    end
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    
    
    all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)
        
    
    Fz_Waves = all_Waves(all_Waves(:,3)==2,:);
    thr_wave = prctile(Fz_Waves(:,9),80);
    Fz_Waves = all_Waves(all_Waves(:,3)==2,:);
    halfw_samples = Fz_Waves(Fz_Waves(:,9)>=thr_wave,6) - Fz_Waves(Fz_Waves(:,9)>=thr_wave,5);
  
    % ^ Half-period in samples --> convert to period in Hz
    
    halfw_samples = Fz_Waves(Fz_Waves(:,9)>=thr_wave,6) - Fz_Waves(Fz_Waves(:,9)>=thr_wave,5);
%     halfw_samples = Fz_Waves(:,6)-Fz_Waves(:,5);
    
    period = 1./(halfw_samples/500); % ^ Half-period in samples --> convert to Hz
    all_period = [all_period ; period];
    all_p2p = [all_p2p; Fz_Waves(Fz_Waves(:,9)>=thr_wave,9)];
       
end

simpleCorPlotsetbin(all_period,all_p2p,0:0.1:10);
xlabel('Period (Hz)')
ylabel('Max Negative Peak Amplitude')
% title('Correlation between Period of Theta(90) Waves and Peak-to-peak Amplitude in Fz(Ch2) ')

%% Threshold Comparison across Participants & P2P


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


% loop across trials for baseline blocks
all_amp_Waves=[];
histc_record = zeros(20,50); % Documenting purposes; rows = participants, columns = bin/counts
collect_all=[];
figure;

for n=1:length(bsl_files)
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    
    nE=2;
    Fz_Waves=all_Waves(all_Waves(:,3)==nE,:);
    thr_Wave=prctile(all_Waves(all_Waves(:,3)==nE,4),80);
    
    
    all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)

    fz_waves = all_Waves(all_Waves(:,3)==2,:);
%     thr_wave=prctile(all_Waves(:,4),80);
    thr_wave=prctile(abs(all_Waves(all_Waves(:,3)==2,9)),80); 

   
    
    all_p2p = Fz_Waves(:,4);
    [N edges] = histcounts(all_p2p, 0:1.5:130);
    
    
    collect_all = [collect_all ; N];
    subplot(5,4,n)
    xlabel(SubID);
    p = plot(edges(1:length(edges)-1),N, '-o','MarkerIndices',round(thr_Wave/1.5));
%     hold on
    p.MarkerFaceColor = 'r';
    p.MarkerSize =7;    
    
    axis([0 130 0 350]);
%     xlabel(SubID);
   
end