%%
clear all
close all

run ../localdef_wanderIM;

addpath(genpath(lscpTools_path));
% addpath(genpath(sleepTools_path));
addpath(genpath(spm12_path));
preproc_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
files=dir([preproc_path filesep 'nfEEG_S*.mat']);

%% General parameters
hb_window=[-0.1 0.4];

%% Loop on files
for n=1:length(files)
    %%% LOAD
    filename=files(n).name;
    subID=filename(findstr(filename,'_')+2:findstr(filename,'.')-1);
    fprintf('... %s\n',subID)
    
    filebehav=dir([behav_path filesep 'wanderIM_behavres_s' subID '*.mat']);
    D=spm_eeg_load([preproc_path filesep filename]);
    
    %%%%% Epoch HeartBeat
    hb_times=detect_heartbeat(D,2,0);
    all_hb_times{n}=hb_times;
%             pretrig  =  hb_window(1) * D.fsample;
%             posttrig =  hb_window(2) * D.fsample;
%             trllabels=[];
%             trl=[];
%             fprintf('... ... %3.0f\n',0)
%             for ntr=1:length(hb_times)
%                 fprintf('\b\b\b\b\b\b\b\b\b\b\b\b... ... %3.0f\n',round(ntr/length(hb_times)*100))
%                 trlbegin = hb_times(ntr) + pretrig;
%                 trlend   = hb_times(ntr) + posttrig;
%                 offset   = pretrig;
%                 newtrl   = [trlbegin trlend offset];
%                 trl      = [trl; newtrl];
%                 trllabels{ntr}='HB';
%             end
%             S=[];
%             S.prefix = 'hb_';
%             S.D=D;
%             S.bc=1;
%             S.trl=trl;
%             S.conditionlabels=trllabels;
%             S.save=1;
%             D_hb=spm_eeg_epochs(S);
end
