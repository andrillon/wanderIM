%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
% close all;
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

%% loop across trials for baseline blocks
all_amp_Waves=[];
for n=1:length(bsl_files)
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);

    load([eeg_path filesep 'wanderIM_twa2_' SubID])
    Fz_Waves=all_Waves(all_Waves(:,3)==2,:);
    wave_period=2*(Fz_Waves(:,6)-Fz_Waves(:,5))/500;
    thr_wave=prctile(Fz_Waves(Fz_Waves(:,3)==2,4),80);
    Fz_thWaves=Fz_Waves(Fz_Waves(:,4)>=thr_wave,:);
    thwave_period=2*(Fz_Waves(Fz_Waves(:,4)>=thr_wave,6)-Fz_Waves(Fz_Waves(:,4)>=thr_wave,5))/500;
    
    [nout,bins]=histc(thwave_period,(0:0.1:5));
    all_freq_thetaWave(n,:)=100*nout/sum(nout);
    for nE=1:63
        thr_Wave1(n,nE)=prctile(all_Waves(all_Waves(:,3)==nE,4),80);
        thr_Wave2(n,nE)=prctile(all_Waves(all_Waves(:,3)==nE,4),90);
    end
    
    for npr=1:60
       num_thWaves_perprobe(n,npr)=sum(Fz_thWaves(:,2)==npr);
       
       mindstate_perprobe(n,npr)=probe_res(npr,32);
       sleepiness_perprobe(n,npr)=probe_res(npr,38);
    end
%                 num_Waves(nE,npr)=length(cell2mat(twa_results.channels(nE).maxnegpkamp));
%             amp_Waves(nE,npr)=mean(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))+abs(cell2mat(twa_results.channels(nE).maxpospkamp)));
end

%%
