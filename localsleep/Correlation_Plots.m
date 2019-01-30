%% Correlation Plot using LSCP tools
% This script contains the code used to make all the correlation plots
% For topographies of the correlation coefficients, see topographies.m

% Contents (Ctrl+F: [])
% 1. [001] Corplot: Number of Theta Waves and Sleepiness
% 2. [002] Barplot: Number of Theta Waves in each mindstate
% 3. [003] Corplot: Theta Number and Performance on Go/Nogo Trial
% 4. [004] Corplot: Theta Number and Reaction Time
% 5. [005] Corplot: P2P amplitude and Max positive/negative peak amplitude
% 6. [006] Corplot: P2P or Maxnegpk amplitude and Sleepiness

%% [001] Corplot theta number and sleepiness
clear all;
 close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
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

for n=1:length(bsl_files)
    
    if n==5 %Exclude SubID 306 (Faulty electrode)
        continue
    end    
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'wanderIM_twa2_' SubID]) %Load results from theta wave detection
    filename=bsl_files(n).name;

    nE = 2; %Designate electrode to analyze (Default is Fz/ch2)
    nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
        halfw_samples = nE_Waves(:,6)-nE_Waves(:,5);
    
    period = 1./(halfw_samples/500); % ^ Half-period in samples --> convert to Hz
    d = find(period>=5 & period<=9);
    nE_Waves = nE_Waves(d,:);
    
    
    
    
    
    thr_Wave = prctile(nE_Waves(:,9),20); %Threshold Wave (80th percentile) by Maxnegpkamp    
    nPr_theta = nE_Waves(nE_Waves(:,9)<=thr_Wave,2);
    nTheta_n = hist(nPr_theta,1:60); % number of Theta waves indexed by probe for the participant in the loop
    nTheta = [nTheta (nTheta_n)];     
    sleepiness_n = (probe_res(:,38)); %sleepiness rating indexed by probe for the participant in the loop    
    sleepiness_n=zscore(sleepiness_n);
    sleepiness = [sleepiness ; sleepiness_n];
    

%     subplot(5,4,n) %Plotting each participant in subplots    
%     simpleCorPlotsetbin(sleepiness_n, nTheta_n,-3:3); 

 
end
    sleepiness = round(sleepiness,1);
figure; %mean across participants
[r p] = simpleCorPlotsetbin(sleepiness,nTheta,-4:0.1:4)
ylabel('Number of Theta (80) Waves')
xlabel('Z-score Alertness (Low to High; Alert to Sleepy)')
title('Correlation between Alertness and Number of Theta Waves Detected in Fz')
set(gca,'Fontsize',20);


%% [002] barplot theta number and mindstate
clear all;
close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
% addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
% addpath(path_localsleep)


% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
% addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);
nTheta = [];
mind_state = [];

nE = 2; %Designate electrode

for n=1:length(bsl_files)
    if n==5
        continue
    end
    
      % load file with spm
    filename=bsl_files(n).name;
    D=spm_eeg_load([eeg_path filesep filename]);
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'wanderIM_twa2_' SubID])
    
    
    nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
            halfw_samples = nE_Waves(:,6)-nE_Waves(:,5);
    
    period = 1./(halfw_samples/500); % ^ Half-period in samples --> convert to Hz
    d = find(period>=5 & period<=9);
    nE_Waves = nE_Waves(d,:);
    
    
    thr_wave = prctile(nE_Waves(:,9),20);% Threshold P2P for Theta Wave 
    nPr_theta = nE_Waves(nE_Waves(:,9)<=thr_wave,2);
    nTheta = [nTheta hist(nPr_theta,1:60)];
    mind_state = [mind_state ; probe_res(:,32)];
end


simpleBarPlot(0.5,nTheta(mind_state==1),'b',1,'r');
simpleBarPlot(2,nTheta(mind_state==2),'b',1,'r');
simpleBarPlot(3.5,nTheta(mind_state==3),'b',1,'r');
simpleBarPlot(5,nTheta(mind_state==4),'b',1,'r');
xticks([0:0.5:7]);
xticklabels({ [] 'ON' [] [] 'MW' [] [] 'MB' [] [] '?'});
xlabel('Mind State');
ylabel('Number of Theta Waves Detected');
title('Number of Theta Waves Detected in Fz in each Mind-state');
set(gca,'Fontsize',20)



%% [003] Corplot: Theta Number and Performance on Go/Nogo Trial

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
    
    if n==5
        continue
    end
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'wanderIM_twa2_' SubID])

for nE = 1:63
    
    thr_wave = prctile(all_Waves(all_Waves(:,3)==nE,9),20);
    nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
    nPr_theta = nE_Waves(nE_Waves(:,9)<=thr_wave,2);
    nTheta(nE,n)= length(nPr_theta);
  
end
end

% loop across trials for baseline blocks
% Vectors for x and y axis
nTheta = [];
avg_nogo=[];



for n=1:length(bsl_files)
    
    if n==5
        continue 
    end
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'wanderIM_twa2_' SubID])

 nE=2;
 nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
 
 halfw_samples = nE_Waves(:,6)-nE_Waves(:,5);
    
    period = 1./(halfw_samples/500); % ^ Half-period in samples --> convert to Hz
    d = find(period>=5 & period<=9);
    nE_Waves = nE_Waves(d,:);

    thr_wave = prctile(nE_Waves(:,9),20);
    nPr_theta = nE_Waves(nE_Waves(:,9)<=thr_wave,2);
    
    nTheta_n = hist(nPr_theta,1:60); 
    nTheta = [nTheta (nTheta_n)]; 
   
    numtrials = probe_res(:,6); % index where probe, but split in 6 blocks
    %end result, avg_go is a 1200 length vector with avg go correct per
    %probe (60) per participant (20)
    
    Pr_play = probe_res(:,6); %Numbers here change depending on Go or Nogo Trial
    for block = 1:6
     block_res = test_res(test_res(:,1)==block,:);
     block_ind = Pr_play(10*block-9:block*10);
     for x = 1:length(block_ind);
         b = find(block_res(:,12)==1 | block_res(:,12)==0); % 12 if Go trials, 11 if NoGo trials
         b=b(b<block_ind(x));
         b=b(end-19:end); % end-19 if Go trials, end-1 if NoGo trials
         b=block_res(b,12); % 12 if Go trials, 11 if NoGo trials
         avg_nogo = [avg_nogo sum(b)];
     end
         
    end
   
    
    
     
    
    
    
    
end

% [r p]= simpleCorPlotsetbin(avg_nogo,nTheta,0:1:2)
[r p]= simpleCorPlotsetbin(avg_nogo,nTheta,0:1:20)
% xlabel('Number of Correct NoGo Responses')
xlabel('Number of Correct Go Responses')
ylabel('Number of Theta Waves Detected')
% title('Correlation between NoGo Performance and Theta Waves Detected in Fz')
title('Correlation between Go Performance and Theta Waves Detected in Fz')
set(gca,'Fontsize',24)
%% [004] Corplot theta number and RT of Go Trials

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
rt=[];




for n=1:length(bsl_files)
    if n==5
        continue 
    end
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'wanderIM_twa2_' SubID])
            
    

    
    nE_Waves = all_Waves(all_Waves(:,3)==2,:);
     halfw_samples = nE_Waves(:,6)-nE_Waves(:,5);
    
    period = 1./(halfw_samples/500); % ^ Half-period in samples --> convert to Hz
    d = find(period>=5 & period<=9);
    nE_Waves = nE_Waves(d,:);
    
    thr_wave = prctile(nE_Waves(:,9),20);
    
    nPr_theta = nE_Waves(nE_Waves(:,9)<=thr_wave,2);
    
    nTheta_n = hist(nPr_theta,1:60); 
    nTheta = [nTheta (nTheta_n)]; 
   
    numtrials = probe_res(:,6); % index where probe, but split in 6 blocks
    %end result, rt is a 1200 length vector with reactione time per
    %probe (63) per participant (20)
    
    Pr_play = probe_res(:,6); %Obtaining reactione time
    for block = 1:6
     block_res = test_res(test_res(:,1)==block,:);
     block_ind = Pr_play(10*block-9:block*10);
     for x = 1:length(block_ind);
         b = find(block_res(:,12)==1 | block_res(:,12)==0);
         b=b(b<block_ind(x));
         b=b(end-19:end);
         b=block_res(b,10)-block_res(b,8);
         rt = [rt nanmean(b)];
     end
         
     end
   
    
     
    
    
    
    
    
end

[r p] = simpleCorPlotsetbin(rt,nTheta,0:0.1:1)
title('Correlation between Reaction Time in Go Trials and Theta Waves Detected')
xlabel('Reaction Time (s)')
ylabel('Number of Theta Waves Detected')
set(gca,'Fontsize',24);



%% [005] Corplot p2p amplitude and maxnegpkamp/maxpospkamp

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

all_p2p=[];
all_maxnegpkamp=[];
percent_pos=[];

for n=1:length(bsl_files)
    
    if n==5
        continue
    end
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])

   
    nE_Waves = all_Waves(all_Waves(:,3)==2,:); %Electrode: Fz/Ch2
    
%     nE_waves = all_Waves;
    
    all_p2p=[all_p2p ; nE_Waves(:,4)];
    all_maxnegpkamp=[all_maxnegpkamp;  nE_Waves(:,9)]; %Column 4 = maxpospkamp, 9 = maxnegpkamp
    
end

% simpleCorPlot(all_p2p, abs(all_maxnegpkamp));

simpleCorPlotsetbin(abs(all_maxnegpkamp), all_p2p, [0:10:100]);

%% INCOMPLETE 
% By-Trial Analysis for Theta Waves
% For every incorrect Go Trial, look -1s and 1s around Trial Onset to see if there is a Theta Wave



for n=1:length(bsl_files)
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'wanderIM_twa2_' SubID])

    thr_wave = prctile(all_Waves(all_Waves(:,3)==2,4),80);
    nE_Waves = all_Waves(all_Waves(:,3)==2,:);
    nPr_theta = nE_Waves(nE_Waves(:,4)>=thr_wave,2);
    
    nTheta_n = hist(nPr_theta,1:60); 
    nTheta = [nTheta (nTheta_n)]; 
   
    numtrials = probe_res(:,6); % index where probe, but split in 6 blocks
    %end result, avg_go is a 1200 length vector with avg go correct per
    %probe (60) per participant (20)
    
    incor_go = find(test_res(:,12)==0);
    incor_onset = test_res(incor_go,8);
         
    
end


     
%% [006] corplot maxnegpkamp and sleepiness
clear all;
 close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
% addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);



% loop across trials for baseline blocks
% Vectors for x and y axis


    
maxnegpkamp = [];
sleepiness = [];

all_corrs = [];
% for n=1:length(bsl_files)
    for n = 1:20
    if n==5
        continue
    end
    
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    filename=bsl_files(n).name;
%     D=spm_eeg_load([eeg_path filesep filename]);

    
    nE = 2; %Eletrode: Fz or Ch2   
    nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
      halfw_samples = nE_Waves(:,6)-nE_Waves(:,5);
    
    period = 1./(halfw_samples/500); % ^ Half-period in samples --> convert to Hz
    d = find(period>=5 & period<=9);
    nE_Waves = nE_Waves(d,:);
    
    thr_wave = prctile(nE_Waves(:,9),20);
    
    nE_Waves = nE_Waves(nE_Waves(:,9)<=thr_wave,:);
    
    for npr = 1:60
        maxnegpkamp = [maxnegpkamp mean(nE_Waves(nE_Waves(:,2)==npr,4))]; %Column 4 for P2P, Column 9 for Maxnegpkamp
    end
    
   
    sleepiness_n = (probe_res(:,38)); %sleepiness rating indexed by probe for the participant in the loop
%     sleepiness_n = zscore(sleepiness_n);
    sleepiness = [sleepiness; sleepiness_n];
    

end
    
[r p] = simpleCorPlotsetbin(sleepiness,maxnegpkamp,0:0.5:5)
% simpleCorPlot(sleepiness(:),maxnegpkamp(:));

% title('Correlation between P2P Amplitude and Alertness in Theta Waves Detected')
title('Correlation between MaxNegPkAmp and Alertness in Theta Waves Detected')
xlabel('Alertness Rating (Higher = Sleepy)')     
ylabel('MaxNegPk Amplitude (uV)')
set(gca,'Fontsize',24)
