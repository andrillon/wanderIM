%%
clear all
close all

run ../localdef_wanderIM;

addpath(genpath(lscpTools_path));
addpath(genpath(spm12_path));
addpath((path_eeglab));
addpath(genpath(path_CSD))

preproc_path=[root_path filesep 'preproc_ica'];
behav_path=[root_path filesep 'behav'];
files=dir([preproc_path filesep 'probe_ica_infEEG_S*.mat']);


%% Loop on files
filename=files(1).name;
D=spm_eeg_load([preproc_path filesep filename]);

myLabels=D.chanlabels(1:63);
[M] = ExtractMontage ([path_CSD filesep 'resource' filesep '10-5-System_Mastoids_EGI129.csd'], myLabels');
[G,H] = GetGH(M);

for n=1:length(files)
    % for n=1:length(files)
    %%% LOAD
    filename=files(n).name;
    subID=filename(findstr(filename,'_')+2:findstr(filename,'.')-1);
    fprintf('... %s\n',subID)
    D=spm_eeg_load([preproc_path filesep filename]);
    
    Dnew = clone(D, ['l' fname(D)], [D.nchannels D.nsamples D.ntrials]);
    Dnew(:,:,:)=D(:,:,:);
    for i = 1:D.ntrials
        Dnew(1:63, :, i) =  CSD(D(1:63, :,i),G, H);
    end
    Dnew = conditions(Dnew, ':', conditions(D, 1:D.ntrials));
    Dnew = repl(Dnew, ':', repl(D, 1:D.ntrials));
    Dnew = events(Dnew, ':', events(D, 1:D.ntrials));
    Dnew = trialonset(Dnew, ':', trialonset(D, 1:D.ntrials));
    Dnew = trialtag(Dnew, ':', trialtag(D, 1:D.ntrials));
    save(Dnew);
end