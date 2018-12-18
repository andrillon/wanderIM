%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
% close all;
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
    load([eeg_path filesep 'wanderIM_twa2_' SubID],'all_Waves')
    
    events = D.events
    for nPr = 1:60
        temp_pr = events{nPr};
        D_type = {temp_pr.type}; %cell array of type (string)
        D_time = {temp_pr.time}; %cell array of time 
        ref_indpr = D_time{ismember(D_type,'P')};
        D_time = cell2mat(D_time(ismember(D_type,'T'))); % Exclude non-Trial Events from times & convert to vec
               
        
%         D_type = (D_type(ismember(D_type,'T'))); % Exclude non-Trial events from types & convert to vec
        D_time = D_time - ref_indpr;                                
        indsamples(1,:) = D_time-1;
        indsamples(2,:) = D_time+1;      
        samples(1,:) = D.indsample(indsamples(1,:))
        samples(isnan(samples)==true)=1;
        samples(2,:) = D.indsample(indsamples(2,:))
        
        nE = 2;
        nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
        thr_Wave=prctile(all_Waves(all_Waves(:,3)==nE,4),80);
        
        nPr_Waves = nE_Waves(nE_Waves(:,2)==nPr & nE_Waves(:,4)>=thr_Wave,5);
        
        logind = [samples(1,:)]<[nPr_Waves] & [samples(2,:)]>[nPr_Waves];
        nTheta = sum(logind);
        
        clear samples
        clear indsamples
    end
end