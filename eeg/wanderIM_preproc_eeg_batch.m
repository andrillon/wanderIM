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

parambatch=[];
parambatch.preproc_path=preproc_path;
parambatch.behav_path=behav_path;

%% General parameters
trial_window=[-0.2 1.3];
hb_window=[-0.1 0.4];
probe_window=[-15 15];

parambatch.trial_window=trial_window;
parambatch.hb_window=hb_window;
parambatch.probe_window=probe_window;

%% Loop on files
parambatch.redo.hp=0;
parambatch.redo.notch=0;
parambatch.redo.epoch.bsl=0;
parambatch.redo.epoch.pr=0;
parambatch.redo.epoch.trl=0;
parambatch.redo.epoch.hb=0;
parpool(8);
parfor n=15:length(files)
% for n=1:length(files)
    %%% LOAD
    filename=files(n).name;
    subID=filename(findstr(filename,'_')+2:findstr(filename,'.')-1);
    fprintf('... %s\n',subID)
    
    filebehav=dir([behav_path filesep 'wanderIM_behavres_s' subID '*.mat']);
    Behav=load([behav_path filesep filebehav.name]);
    D=spm_eeg_load([preproc_path filesep filename]);
    
    wanderIM_preproc_eeg_subbatch(D,Behav,parambatch);
end
delete(gcp('nocreate'))
%%
%
% figure
% plot(erp_time,mean(erp_try(2,:,nogo1),3),'r')
% hold on
% plot(erp_time,mean(erp_try(2,:,go1),3),'b')
% plot(erp_time,mean(erp_try(2,:,go0),3),'b--')
% plot(erp_time,mean(erp_try(2,:,nogo0),3),'r--')
