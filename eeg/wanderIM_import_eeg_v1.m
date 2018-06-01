%%
clear all
close all

localdef_wanderIM;

addpath(genpath(lscpTools_path));
addpath(genpath(spm12_path));
eeg_path=[root_path filesep 'eeg'];
preproc_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'eeg'];

files=dir([eeg_path filesep 'MWI3*.eeg']);

%% Loop on files
for n=1:length(files)
    subID=files(n).name;
    subID=subID(4:6);
    if exist([eeg_path filesep 'EEG_S' subID '.mat'])==0
        fprintf('... importing in SPM format %s \n',subID)
        %%% Get headers and events
        S = [];
        S.dataset = [eeg_path filesep files(n).name];
        S.outfile = [preproc_path filesep 'EEG_S' subID];
        S.channels = 'all';
        S.timewindow = [];
        S.blocksize = 3276800;
        S.checkboundary = 1;
        S.usetrials = 1;
        S.datatype = 'float32-le';
        S.eventpadding = 0;
        S.saveorigheader = 0;
        S.conditionlabel = {'Undefined'};
        S.inputformat = [];
        S.continuous = true;
        S.autoloc = true;
        S.units='uV'; % it will be lost at montage anyway...
        D = spm_eeg_convert(S);
    else
        fprintf('... loading %s \n',subID)
        D=spm_eeg_load([preproc_path filesep 'EEG_S' subID '.mat']);
    end
end

