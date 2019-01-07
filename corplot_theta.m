%% Correlation Plot using LSCP tools
% X: Number of theta waves in Fz for that probe
% Y: Can be anything, sleepiness and mindstate e.g.

%% corplot theta number and sleepiness
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
% nTheta_time=[];
nTheta=[];
% sleepiness_time=[];
all_corrs = [];
for n=1:length(bsl_files)
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    filename=bsl_files(n).name;
%     D=spm_eeg_load([eeg_path filesep filename]);
    
    corr_per=[];

    for nE = 2
    nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
    thr_Wave = prctile(nE_Waves(:,4),80);
    
    nPr_theta = nE_Waves(nE_Waves(:,4)>=thr_Wave,2);
    
    nTheta_n = hist(nPr_theta,1:60); % number of Theta waves indexed by probe for the participant in the loop
%     nTheta_n = nTheta_n./mean(nTheta_n);
    nTheta = [nTheta (nTheta_n)]; 
    
    sleepiness_n = (probe_res(:,38)); %sleepiness rating indexed by probe for the participant in the loop
    
    sleepiness = [sleepiness ; sleepiness_n];
    
    
    
%   
%     [rho pv] = corr(nTheta, sleepiness, 'type', 'spearman'); 
%     corr_per = [corr_per rho];
%     nTheta = [];
%     sleepiness = [];
    end 
    all_corrs = [ all_corrs;corr_per];

    
%     [r,pV]=corr(hist(nPr_theta,1:60)', probe_res(:,38),'type','spearman');
%     allcorr=[allcorr ; [r pV]];
%     
%     for nE=1:63
%         nTheta_byCh(n,nE)=sum(all_Waves(all_Waves(:,3)==nE,4)>=thr_wave);
%         [r,pV]=corr(hist(nPr_theta,1:60)', probe_res(:,38),'type','spearman');
%         
%     end
%     
%     c=0;
%     for nprctile=0:10:90
%         c=c+1;
%         tp_thr_wave = prctile(all_Waves(:,4),nprctile); % Threshold P2P for Theta Wave (Nth percentile)
%         nPr_theta = nE_waves(nE_waves(:,4)>=tp_thr_wave,2);
%         [r_byprctile(n,c),pV]=corr(hist(nPr_theta,1:60)', probe_res(:,38),'type','spearman');
%     end
end
    




%%


    
% figure;
format_fig;
addpath(genpath(path_eeglab));
temp_topo=(mean(vec2mat(all_corrs,63))); % vector of 63 values
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
% caxis([0 2.5])
rmpath(genpath(path_eeglab));
colorbar;
caxis([min(temp_topo) 0.14])

 
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
for nE = 2
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
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)   
    
    thr_wave = prctile(all_Waves(all_Waves(:,3)==nE,9),80);% Threshold P2P for Theta Wave (Nth percentile)
    nE_waves = all_Waves(all_Waves(:,3)==nE,:);
    nPr_theta = nE_waves(nE_waves(:,4)>=thr_wave,2);
    nTheta = [nTheta hist(nPr_theta,1:60)];
    mind_state = [mind_state ; probe_res(:,32)];
end

% subplot(7,9,nE)
simpleBarPlot(0.5,nTheta(mind_state==1),'b',1,'r');
simpleBarPlot(2,nTheta(mind_state==2),'b',1,'r');
simpleBarPlot(3.5,nTheta(mind_state==3),'b',1,'r');
% simpleBarPlot(5,nTheta(mind_state==4),'b',1,'r');
xticks([0:0.5:7]);
xticklabels({ [] 'ON' [] [] 'MW' [] [] 'MB' [] [] '?'});
xlabel(D.chanlabels(nE))
ylabel('Number of Theta Waves Detected');
% title('Number of Theta (80) Waves Detected in each Mind-state');
% nTheta = [];
% mind_state = [];
end


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

for n=[1:20]
    
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])

for nE = 1:63
    
    thr_wave = prctile(all_Waves(all_Waves(:,3)==nE,4),80);
    nE_waves = all_Waves(all_Waves(:,3)==nE,:);
    nPr_theta = nE_waves(nE_waves(:,4)>=thr_wave,2);
    nTheta(nE,n)= length(nPr_theta);
  
end
end

%%

for n=1:20
subplot(5,4,n);
    
% figure;
format_fig;
addpath(genpath(path_eeglab));
temp_topo=mean(nTheta,2); % vector of 63 values
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
% caxis([0 2.5])
rmpath(genpath(path_eeglab));
colorbar;
caxis([min(temp_topo) max(temp_topo)])
end

%%
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
    if n==5
        continue 
    end
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)   
    
    
    

    thr_wave = prctile(all_Waves(all_Waves(:,3)==2,9),80);
    nE_waves = all_Waves(all_Waves(:,3)==2,:);
    nPr_theta = nE_waves(nE_waves(:,9)>=thr_wave,2);
    
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
%  close all;
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
         avg_go = [avg_go; nanmean(b)];
     end
         
    end
     
    
    sleepiness_n = (probe_res(:,38)); %sleepiness rating indexed by probe for the participant in the loop
    
    sleepiness = [sleepiness  sleepiness_n];
   
        
    
end

%% Corplot p2p amplitude and maxnegpkamp/maxpospkamp; Barplot for %P2PA due to Pos/Neg cpomonent; see slide 45

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

    

    
    all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)

    fz_waves = all_Waves(all_Waves(:,3)==2,:);

    negzx_all = fz_waves(:,5);
    npr_all = fz_waves(:,2);

    npr_all(negzx_all>=9250 | negzx_all<=750)=[];
    negzx_all(negzx_all>=9250)=[];
    negzx_all(negzx_all<=750)=[];
    
    
    
    
    
    
    nE_waves = all_Waves(all_Waves(:,3)==2,:);
    
%     nE_waves = all_Waves;
    
    all_p2p=[all_p2p ; nE_waves(:,4)];
    all_maxnegpkamp=[all_maxnegpkamp;  nE_waves(:,9)];
    percent_pos=[percent_pos  mean(abs((nE_waves(:,9)./nE_waves(:,4))).*100)]; 
    
end

% simpleCorPlot(all_p2p, abs(all_maxnegpkamp));

simpleCorPlotsetbin(abs(all_maxnegpkamp), all_p2p, [0:10:100]);

%% Examining how much maxpospkamp attributes to p2p amplitude






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



%% Corplot b/n nTheta and Number of correct "Go"

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

nSleepiness=[];
avg_go=[];

allcorr=[];

all_nTheta=[];
for n=1:length(bsl_files)
    if n==5
        continue 
    end
    tempb=[];
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])

    nSleepiness = [ nSleepiness; probe_res(:,38)];

    numtrials = probe_res(:,6); % index where probe, but split in 6 blocks
    %end result, avg_go is a 1200 length vector with avg go correct per
    %probe (63) per participant (20)
    
    Pr_play = probe_res(:,6);
    for block = 1:6
     block_res = test_res(test_res(:,1)==block,:);
     block_ind = Pr_play(10*block-9:block*10);
     for x = 1:length(block_ind);
         b = find(block_res(:,11)==1 | block_res(:,11)==0);
         b=b(b<block_ind(x));
         b=b(end-1:end);
         b=block_res(b,11);
         avg_go = [avg_go mean(b)]; %Consider median instead of mean
         tempb=[tempb mean(b)];
     end
         
    end
     
     all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)   
    
        
           for nE = 2
        nE_wave = all_Waves(all_Waves(:,3)==nE,:);
        thr_wave = prctile(all_Waves(all_Waves(:,3)==nE,9),80);
        nE_wave = nE_wave(nE_wave(:,9)>=thr_wave,:);
        a = zeros(60,1);
        for nPr = 1:60
%             a(nPr) =  mean(nE_wave(nE_wave(:,2)==nPr,4));
              a(nPr) = length(nE_wave(nE_wave(:,2)==nPr,1));  
%             if isnan(a(nPr))==true
%                 a(nPr)=0;
%             end
        end
        
        
        
           
            all_nTheta = [all_nTheta ; a];
           end
           
%     subplot(5,4,n)
%     simpleCorPlotsetbin(tempb,a,0:1:10);
           
end
% 
simpleCorPlotsetbin(avg_go,all_nTheta,0:0.1:10);
xlabel('Average Number of Correct NoGo Trial')
ylabel('Number of Theta Waves')
%     
%% corplot maxnegpkamp and sleepiness
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
    filename=bsl_files(n).name;
%     D=spm_eeg_load([eeg_path filesep filename]);
    
    corr_per=[];

    all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)


    
nE = 2;
  
    
    nE_waves = all_Waves(all_Waves(:,3)==nE,:);

    thr_wave=prctile(abs(all_Waves(all_Waves(:,3)==nE,4)),80);
    

    
    for npr = 1:60
        maxnegpkamp = [maxnegpkamp mean(nE_waves(nE_waves(:,2)==npr,4))];
    end
    
   
    sleepiness_n = (probe_res(:,38)); %sleepiness rating indexed by probe for the participant in the loop
    
    sleepiness = [sleepiness ; sleepiness_n];
    
    for nE = 1:63
           nE_waves = all_Waves(all_Waves(:,3)==nE,:);
           thr_wave=prctile(abs(all_Waves(all_Waves(:,3)==nE,4)),80);  
           maxnegpkcorr=[];
            for npr = 1:60
                maxnegpkcorr = [maxnegpkcorr mean(nE_waves(nE_waves(:,2)==npr,4))];
            end
        [r,pV]=corr(maxnegpkcorr', probe_res(:,38),'type','spearman');
        all_corrs = [all_corrs r];
    end

end
    
% simpleCorPlotsetbin(sleepiness,abs(maxnegpkamp),1:4);
% xlabel('Alertness Rating')
% ylabel('maxnegpkamp (abs)')
     
