%%
clear all
close all

run ../localdef_wanderIM;

addpath(genpath(lscpTools_path));
addpath(genpath(sleepTools_path));
addpath(genpath(spm12_path));
preproc_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
files=dir([preproc_path filesep 'EEG_S*.mat']);

%% General parameters
trial_window=[-0.2 1.3];
hb_window=[-0.1 0.4];
probe_window=[-15 15];

%% Loop on files
redo.hp=1;
redo.notch=0;

for n=16:length(files)
    %%% LOAD
    filename=files(n).name;
    subID=filename(findstr(filename,'_')+2:findstr(filename,'.')-1);
    fprintf('... %s\n',subID)
    
    filebehav=dir([behav_path filesep 'wanderIM_behavres_s' subID '*.mat']);
    Behav=load([behav_path filesep filebehav.name]);
    D=spm_eeg_load([preproc_path filesep filename]);
    
   
    %%% High-pass filter
    if exist([preproc_path filesep 'f' D.fname])==0 || redo.hp==1
        % update channel type for ECG
        D = chantype(D, match_str(D.chanlabels,'ECG'), 'ECG');
        D.save;
        
        
        type = 'butterworth';
        order = 5;
        dirfilt = 'twopass';
        S = [];
        S.D = D;
        S.band = 'high';
        S.type = type;
        S.order = order;
        S.dir = dirfilt;
        S.freq = 0.1;
        S.save=1;
        D = spm_eeg_filter(S);
    else
        D=spm_eeg_load([preproc_path filesep 'f' D.fname]);
    end
    
    %%% NOTCH FILTERS
    if exist([preproc_path filesep 'n' D.fname])==0 || redo.notch==1
        type = 'butterworth';
        order = 4;
        dirfilt = 'twopass';
        S = [];
        S.prefix='n';
        S.D = D;
        S.band = 'stop';
        S.type = type;
        S.order = order;
        S.dir = dirfilt;
        S.freq = [45 55];
        S.save=1;
        D = spm_eeg_filter(S);
    else
        D=spm_eeg_load([preproc_path filesep 'n' D.fname]);
    end
end
%%
%
% figure
% plot(erp_time,mean(erp_try(2,:,nogo1),3),'r')
% hold on
% plot(erp_time,mean(erp_try(2,:,go1),3),'b')
% plot(erp_time,mean(erp_try(2,:,go0),3),'b--')
% plot(erp_time,mean(erp_try(2,:,nogo0),3),'r--')
