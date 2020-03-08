%% Convert EEG ICA data back to spm
clear all
close all
clc;

%%% dependencies
path_eeglab='/Users/tand0009/Work/local/eeglab';
path_SPM12='/Users/tand0009/Work/local/spm12/';
addpath(path_eeglab);
% eeglab;
addpath(genpath(path_SPM12));

path_data_eeglab='/Users/tand0009/Data/WanderIM/preproc_ica/';
path_data_spm=path_data_eeglab;

%%% list files
listfile_eeglab=dir([path_data_eeglab filesep 'probe_infEEG_S*_ica_clean.set']);
listfile_spm=dir([path_data_spm filesep 'probe_infEEG_S*.mat']);

%% Loop across files
for nF=1:length(listfile_eeglab)
    %%% load EEG lab dataset
    filename_eeglab=listfile_eeglab(nF).name;
    EEG = pop_loadset('filename',filename_eeglab,'filepath',path_data_eeglab);
    fprintf('... loading %s (%g/%g)\n',filename_eeglab,nF,length(listfile_eeglab))
    %%% Copy SPM structure
    filename_spm=listfile_spm(nF).name;
    if strcmp(filename_spm(1:17),filename_eeglab(1:17))
        D=spm_eeg_load([path_data_spm filesep filename_spm]);
        S = [];
        S.D = D;
        new_smp_fname = fullfile(path_data_spm,[filename_spm(1:6) 'ica_' filename_spm(7:end)]);
        S.outfile = new_smp_fname;
        D = spm_eeg_copy(S);
        fprintf('... converting to %s (%g/%g)\n',new_smp_fname,nF,length(listfile_eeglab))
        
        % updating meeg-object with new data
        %         chanlabels={EEG.chanlocs.labels};
        %         for nChan=1:length(chanlabels)
        %             D(match_str(D.chanlabels,chanlabels{nChan}),:,:) = EEG.data(nChan,:,:);
        %         end
        D(1:63,:,:) = EEG.data(:,:,:);
        D.save;
        fprintf('... done!\n')
        
    end
end