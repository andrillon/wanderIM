%% Correlation Plot using LSCP tools
% X: Number of theta waves in Fz for that probe
% Y: Can be anything, sleepiness and mindstate e.g.

%% Init
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
% Vectors for x and y axis
nTheta = [];
sleepiness = [];

for n=1:length(bsl_files)
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'wanderIM_twa2_' SubID])
    
    thr_wave = prctile(all_Waves(all_Waves(:,3)==2,4),80); % Threshold P2P for Theta Wave (Nth percentile)
    Fz_waves = all_Waves(all_Waves(:,3)==2,:);
    nPr_theta = Fz_waves(Fz_waves(:,4)>=thr_wave,2);
    a = hist(nPr_theta,1:60);
    b = find(probe_res(:,5)==1);
    a(b) = NaN;
    nTheta = [nTheta a];
    c = probe_res(:,38);
    c(b) = NaN;
    sleepiness = [sleepiness ; c];
    
end
simpleCorPlotsetbin(nTheta,sleepiness, 1:15);
xlabel('Number of Theta (80) Waves')
ylabel('Alertness (1 (alert) to 4 (sleepy))')




