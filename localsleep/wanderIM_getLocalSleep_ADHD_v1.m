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
eeg_path=[root_path_adhd filesep 'preproc_eeg'];
behav_path=[root_path_adhd filesep 'behav'];
bsl_files=dir([eeg_path filesep 'nfEEG_S4*.mat']);

%% loop across trials for baseline blocks
for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    D=spm_eeg_load([eeg_path filesep filename]);
    fprintf('... processing subject %s\n',D.fname)
    
    
    % load behavioural results
    SubID=D.fname;
    SubID=SubID(findstr(SubID,'_S4')+2:findstr(SubID,'_S4')+4);
    behav_file=dir([behav_path filesep 'MWADHD_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    if exist([eeg_path filesep 'MWADHD_cont_twa2_' SubID '.mat'])==0
        %     left_freq(n)=SubjectInfo.FlickerR; % CAREFUL this is inverted
        %     right_freq(n)=SubjectInfo.FlickerL;
        
        temp_data=D(1:64,:,:); % D contains the data with channels * time * trials
        temp_data=temp_data-repmat(mean(temp_data([17 22],:,:),1),[size(temp_data,1) 1 1]); % D contains the data with channels * time * trials
        
        [twa_results]=twalldetectnew_TA(temp_data,D.fsample,0);
        all_Waves=[];
        for nE=1:64
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
        save([eeg_path filesep 'MWADHD_cont_twa2_' SubID],'all_Waves')
        for nE=1:64
            thr_Wave1(n,nE)=prctile(all_Waves(all_Waves(:,3)==nE,4),80);
            thr_Wave2(n,nE)=prctile(all_Waves(all_Waves(:,3)==nE,4),90);
        end
    end
    %                 num_Waves(nE,npr)=length(cell2mat(twa_results.channels(nE).maxnegpkamp));
    %             amp_Waves(nE,npr)=mean(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))+abs(cell2mat(twa_results.channels(nE).maxpospkamp)));
end

%%
