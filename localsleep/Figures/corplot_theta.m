%% Correlation Plot using LSCP tools
% X: Number of theta waves in Fz for that probe
% Y: Can be anything, sleepiness and mindstate e.g.

%% corplot theta number and sleepiness
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
% Vectors for x and y axis
nTheta = [];
sleepiness = [];

allcorr=[];
nTheta_time=[];
sleepiness_time=[];
for n=1:length(bsl_files)
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'wanderIM_twa2_' SubID])

    a = [2 7 24 29 35 39 56 62 63 3 30 38 57 8 25 12 52 23 13]; %% Good until from 2 to 63
    b=[];
    for x = 1:length(a)
        b = [b ; all_Waves(all_Waves(:,3)==a(x),4)];
    end
    thr_wave = prctile(b,90); % Threshold P2P for Theta Wave (Nth percentile)
    
    Fz_waves = all_Waves(all_Waves(:,3)==2,:);
    nPr_theta = Fz_waves(Fz_waves(:,4)>=thr_wave,2);
    
    nTheta_n = hist(nPr_theta,1:60); % number of Theta waves indexed by probe for the participant in the loop
%     nTheta_n = nTheta_n./mean(nTheta_n);
    nTheta = [nTheta (nTheta_n)]; 
    
    sleepiness_n = (probe_res(:,38)); %sleepiness rating indexed by probe for the participant in the loop
    
    sleepiness = [sleepiness ; sleepiness_n];
    [r,pV]=corr(hist(nPr_theta,1:60)', probe_res(:,38),'type','spearman');
    allcorr=[allcorr ; [r pV]];
    
    for nE=1:63
        nTheta_byCh(n,nE)=sum(all_Waves(all_Waves(:,3)==nE,4)>=thr_wave);
        [r,pV]=corr(hist(nPr_theta,1:60)', probe_res(:,38),'type','spearman');
        
    end
    
    c=0;
    for nprctile=0:10:90
        c=c+1;
        tp_thr_wave = prctile(all_Waves(:,4),nprctile); % Threshold P2P for Theta Wave (Nth percentile)
        nPr_theta = Fz_waves(Fz_waves(:,4)>=tp_thr_wave,2);
        [r_byprctile(n,c),pV]=corr(hist(nPr_theta,1:60)', probe_res(:,38),'type','spearman');
    end
end

%%
figure;
format_fig;
addpath(genpath(path_eeglab));
temp_topo=mean(nTheta_byCh); % vector of 63 values
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
% caxis([0 2.5])
rmpath(genpath(path_eeglab));
colorbar;
caxis([0 max(temp_topo)])
 
%%
figure
simpleCorPlotsetbin(sleepiness,nTheta, 1:4);
% simpleCorPlot(sleepiness,nTheta,[],'pearson');
ylabel('Number of Theta (80) Waves')
xlabel('Alertness (1 (alert) to 4 (sleepy))')


%% barplot theta number and mindstate
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
% Vectors for x and y axis
nTheta = [];
mind_state = [];

for n=1:length(bsl_files)
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'wanderIM_twa2_' SubID])
    
%     thr_wave = prctile(all_Waves(all_Waves(:,3)==2,4),80); % Threshold P2P for Theta Wave (Nth percentile)
    thr_wave = prctile(all_Waves(:,4),90); % Threshold P2P for Theta Wave (Nth percentile)
    Fz_waves = all_Waves(all_Waves(:,3)==2,:);
    nPr_theta = Fz_waves(Fz_waves(:,4)>=thr_wave,2);
    nTheta = [nTheta hist(nPr_theta,1:60)];
    mind_state = [mind_state ; probe_res(:,32)];
end


simpleBarPlot(0.5,nTheta(mind_state==1),'b',1,'r');
simpleBarPlot(2,nTheta(mind_state==2),'b',1,'r');
simpleBarPlot(3.5,nTheta(mind_state==3),'b',1,'r');
% simpleBarPlot(5,nTheta(mind_state==4),'b',1,'r');
xticks([0:0.5:7]);
xticklabels({ [] 'ON' [] [] 'MW' [] [] 'MB' [] [] '?'});
ylabel('Number of Theta Waves Detected');
title('Number of Theta (80) Waves Detected in each Mind-state');



