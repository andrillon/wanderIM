%%
clear all
close all

run ../localdef_wanderIM;

addpath(genpath(lscpTools_path));
addpath(genpath(sleepTools_path));
addpath(genpath(spm12_path));
preproc_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
files=dir([preproc_path filesep 'nfEEG_S*.mat']);

%% General parameters
trial_window=[-0.2 1.3];
hb_window=[-0.1 0.4];
probe_window=[-15 15];

%% Loop on files
redo=0;
for n=1:length(files)
    %%% LOAD
    filename=files(n).name;
    subID=filename(findstr(filename,'_')+2:findstr(filename,'.')-1);
    fprintf('... %s\n',subID)
    
    filebehav=dir([behav_path filesep 'wanderIM_behavres_s' subID '*.mat']);
    Behav=load([behav_path filesep filebehav.name]);
    D=spm_eeg_load([preproc_path filesep filename]);
    
   
    %%%%% Epoch HeartBeat
    if exist([preproc_path filesep 'hb_' D.fname])==0 || redo==1
        hb_times=detect_heartbeat(D,0);
        
        pretrig  =  hb_window(1) * D.fsample;
        posttrig =  hb_window(2) * D.fsample;
        trllabels=[];
        trl=[];
        fprintf('... ... %3.0f\n',0)
        for ntr=1:length(hb_times)
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b... ... %3.0f\n',round(ntr/length(hb_times)*100))
            trlbegin = hb_times(ntr) + pretrig;
            trlend   = hb_times(ntr) + posttrig;
            offset   = pretrig;
            newtrl   = [trlbegin trlend offset];
            trl      = [trl; newtrl];
            trllabels{ntr}='HB';
        end
        S=[];
        S.prefix = 'hb_';
        S.D=D;
        S.bc=1;
        S.trl=trl;
        S.conditionlabels=trllabels;
        S.save=1;
        D_hb=spm_eeg_epochs(S);
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
