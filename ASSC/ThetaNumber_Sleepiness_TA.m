
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
path_thetawave='/Users/tand0009/WorkGit/projects/inprogress/localsleep/thetawavedetection/Results';



% loop across trials for baseline blocks
% Vectors for x and y axis


    
nTheta = [];
sleepiness = [];

for n=1:20 %length(bsl_files)
    
    if n==5 %Exclude SubID 306 (Faulty electrode)
        continue
    end    
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([path_thetawave filesep 'wanderIM_cont_twa2_' SubID]) %Load results from theta wave detection
    filename=bsl_files(n).name;

    nE = 2; %Designate electrode to analyze (Default is Fz/ch2)
    nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
%         halfw_samples = nE_Waves(:,6)-nE_Waves(:,5);
%     
%     period = 1./(halfw_samples/500); % ^ Half-period in samples --> convert to Hz
%     d = find(period>=5 & period<=9);
%     nE_Waves = nE_Waves(d,:); 
    
    
    
    
    
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
%     sleepiness = round(sleepiness,1);
    
nTheta = mean(reshape(nTheta,[60 19])');
sleepiness = mean(reshape(sleepiness,[60 19])');
% sleepiness=zscore(sleepiness);

%%
figure; %mean across participants
[r p] = simpleCorPlotbin(sleepiness,nTheta,{'o','b','b',72,12});
ylabel('Num. of Frontal Waves')
xlabel({'Subj. Sleepiness (zscore)'})
title('Fz')
set(gca,'Fontsize',20);
% axis([-2 2 4 6])



