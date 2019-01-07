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
for n=1
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'fwanderIM_twa2_' SubID])

%     a = [2 7 24 29 35 39 56 62 63 3 30 38 57 8 25 12 52 23 13]; %% Good until from 2 to 63
    a=2;
    b=[];
    for x = 1:length(a)
        b = [b ; all_Waves(all_Waves(:,3)==a(x),4)];
    end
    
%     thr_wave = prctile(b,80); % Threshold P2P for Theta Wave (Nth percentile)
    thr_wave = prctile(all_Waves(all_Waves(:,3)==2,4),80);
    nE_waves = all_Waves(all_Waves(:,3)==2,:);
    nPr_theta = nE_waves(nE_waves(:,4)>=thr_wave,2);
    
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
        nPr_theta = nE_waves(nE_waves(:,4)>=tp_thr_wave,2);
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
    load([eeg_path filesep 'fwanderIM_twa2_' SubID])
    
%     thr_wave = prctile(all_Waves(all_Waves(:,3)==2,4),80); % Threshold P2P for Theta Wave (Nth percentile)
    thr_wave = prctile(all_Waves(all_Waves(:,3)==2,4),90);% Threshold P2P for Theta Wave (Nth percentile)
    nE_waves = all_Waves(all_Waves(:,3)==2,:);
    nPr_theta = nE_waves(nE_waves(:,4)>=thr_wave,2);
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


%% Count Theta Number: Obtain 63 by 20 matrix
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
nTheta = zeros(63,20);

for n=1:length(bsl_files)
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'wanderIM_twa2_' SubID])

for nE = 1:63
    
    thr_wave = prctile(all_Waves(all_Waves(:,3)==nE,4),80);
    nE_waves = all_Waves(all_Waves(:,3)==nE,:);
    nPr_theta = nE_waves(nE_waves(:,4)>=thr_wave,2);
    nTheta(nE,n)= length(nPr_theta);
  
end
end

%%
figure;
format_fig;
addpath(genpath(path_eeglab));
temp_topo=mean(nTheta,2); % vector of 63 values
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
% caxis([0 2.5])
rmpath(genpath(path_eeglab));
colorbar;
caxis([min(temp_topo) max(temp_topo)])


% loop across trials for baseline blocks
% Vectors for x and y axis
nTheta = [];
avg_go=[];

allcorr=[];


for n=1:length(bsl_files)
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'fwanderIM_twa2_' SubID])

    

    thr_wave = prctile(all_Waves(all_Waves(:,3)==2,4),80);
    nE_waves = all_Waves(all_Waves(:,3)==2,:);
    nPr_theta = nE_waves(nE_waves(:,4)>=thr_wave,2);
    
    nTheta_n = hist(nPr_theta,1:60); 
    nTheta = [nTheta (nTheta_n)]; 
   
    numtrials = probe_res(:,6); % index where probe, but split in 6 blocks
    %end result, avg_go is a 1200 length vector with avg go correct per
    %probe (60) per participant (20)
    
    Pr_play = probe_res(:,6);
    for block = 1:6
     block_res = test_res(test_res(:,1)==block,:);
     block_ind = Pr_play(10*block-9:block*10);
     for x = 1:length(block_ind);
         b = find(block_res(:,11)==1 | block_res(:,11)==0);
         b=b(b<block_ind(x));
         b=b(end-1:end);
         b=block_res(b,11);
         avg_go = [avg_go sum(b)];
     end
         
     end
   
    
     
    
    
    
    
    
end

simpleCorPlotsetbin(avg_go,nTheta,0:0.5:2);

%% Corplot theta number and RT of Go Trials

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
avg_go=[];

allcorr=[];


for n=1:length(bsl_files)
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'fwanderIM_twa2_' SubID])

    

    thr_wave = prctile(all_Waves(all_Waves(:,3)==2,4),80);
    nE_waves = all_Waves(all_Waves(:,3)==2,:);
    nPr_theta = nE_waves(nE_waves(:,4)>=thr_wave,2);
    
    nTheta_n = hist(nPr_theta,1:60); 
    nTheta = [nTheta (nTheta_n)]; 
   
    numtrials = probe_res(:,6); % index where probe, but split in 6 blocks
    %end result, avg_go is a 1200 length vector with avg go correct per
    %probe (63) per participant (20)
    
    Pr_play = probe_res(:,6);
    for block = 1:6
     block_res = test_res(test_res(:,1)==block,:);
     block_ind = Pr_play(10*block-9:block*10);
     for x = 1:length(block_ind);
         b = find(block_res(:,12)==1 | block_res(:,12)==0);
         b=b(b<block_ind(x));
         b=b(end-19:end);
         b=block_res(b,10)-block_res(b,8);
         avg_go = [avg_go nanmean(b)];
     end
         
     end
   
    
     
    
    
    
    
    
end

simpleCorPlotsetbin(avg_go,nTheta,0:0.1:1);

%% Plot RT by itself
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
avg_go=[];

allcorr=[];
sleepiness =[];
sleepiness_all=[];

for n=1:length(bsl_files)
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'fwanderIM_twa2_' SubID])

   
   
    numtrials = probe_res(:,6); % index where probe, but split in 6 blocks
    %end result, avg_go is a 1200 length vector with avg go correct per
    %probe (60) per participant (20)
    
    Pr_play = probe_res(:,6);
    for block = 1:6
     block_res = test_res(test_res(:,1)==block,:);
     block_ind = Pr_play(10*block-9:block*10);
     for x = 1:length(block_ind);
         b = find(block_res(:,12)==1 | block_res(:,12)==0);
         b=b(b<block_ind(x));
         b=b(end-19:end);
         b=block_res(b,10)-block_res(b,8);
         avg_go = [avg_go; nanmedian(b)];
     end
         
    end
     
    
    sleepiness_n = (probe_res(:,38)); %sleepiness rating indexed by probe for the participant in the loop
    
    sleepiness = [sleepiness  sleepiness_n];
   
        
    
end
%% By-Trial Analysis for Theta Waves
% For every incorrect Go Trial, look -1s and 1s around Trial Onset to see if there is a Theta Wave



for n=1:length(bsl_files)
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'wanderIM_twa2_' SubID])

    thr_wave = prctile(all_Waves(all_Waves(:,3)==2,4),80);
    nE_waves = all_Waves(all_Waves(:,3)==2,:);
    nPr_theta = nE_waves(nE_waves(:,4)>=thr_wave,2);
    
    nTheta_n = hist(nPr_theta,1:60); 
    nTheta = [nTheta (nTheta_n)]; 
   
    numtrials = probe_res(:,6); % index where probe, but split in 6 blocks
    %end result, avg_go is a 1200 length vector with avg go correct per
    %probe (60) per participant (20)
    
    incor_go = find(test_res(:,12)==0);
    incor_onset = test_res(incor_go,8);
         
    
end



