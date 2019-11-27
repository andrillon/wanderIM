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
addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'nfEEG_S3*.mat']);

%% loop across trials for baseline blocks
all_amp_Waves=[];
redo=0;
for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    D=spm_eeg_load([eeg_path filesep filename]);
    fprintf('... processing subject %s\n',D.fname)
    
    
    % load behavioural results
    SubID=D.fname;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    if exist([eeg_path filesep 'wanderIM_cont_twa2_' SubID '.mat'])==0 || redo==1
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);

    temp_data=D(1:63,:,:); % D contains the data with channels * time * trials
    temp_data=temp_data-repmat(mean(temp_data([10 21],:,:),1),[size(temp_data,1) 1 1]); % D contains the data with channels * time * trials

    [twa_results]=twalldetectnew_TA(temp_data,D.fsample,0);
        all_Waves=[];
        for nE=1:63
            all_Waves=[all_Waves ; [repmat([n 0 nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
                cell2mat(twa_results.channels(nE).negzx)' ...
                cell2mat(twa_results.channels(nE).poszx)' ...
                cell2mat(twa_results.channels(nE).wvend)' ...
                cell2mat(twa_results.channels(nE).maxnegpk)' ...
                cell2mat(twa_results.channels(nE).maxnegpkamp)' ...
                cell2mat(twa_results.channels(nE).maxpospk)' ...
                cell2mat(twa_results.channels(nE).maxpospkamp)' ...
                cell2mat(twa_results.channels(nE).mxdnslp)' ...
                cell2mat(twa_results.channels(nE).mxupslp)' ...
                ]];
        end
    save([eeg_path filesep 'wanderIM_cont_twa2_' SubID],'all_Waves')

    end
%                 num_Waves(nE,npr)=length(cell2mat(twa_results.channels(nE).maxnegpkamp));
%             amp_Waves(nE,npr)=mean(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))+abs(cell2mat(twa_results.channels(nE).maxpospkamp)));
end

%%
    