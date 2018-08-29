%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Teigane's general notes:
%D.conditions - probe conditions
%D.chantype - shows where each channel is coming from - general)
%D.channels - all 66 channels with information about them
%D.fsample - sampling rate
%D.time - time within the epoch - we have from -32 seconds to 32 seconds. - this is the variable we will need to modify to get the right time window for our analysis

%D.indsample - (stands for index of the sample) is a function that will cut
%the data to the time window I need. To be modified when looking at
%baseline vs around probe (baseline cut into 3 10s blocks)

%% Init
clear all;
close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'lprobe_nfEEG_S3*.mat']);

%% loop across trials for baseline blocks
all_probes_mat=[];
for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    D=spm_eeg_load([eeg_path filesep filename]);
    fprintf('... processing subject %s\n',D.fname)
    
    
    % load behavioural results
    SubID=D.fname;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    
    % I probably won't use the left and right flickers individually for now
    % but note: the left and right are showing inverse to logic (needs to
    % be checked)
    left_freq=SubjectInfo.FlickerL;
    right_freq=SubjectInfo.FlickerR;
    
    param=[];
    param.method='fft'; % fast fourier transform
    param.mindist=1; % we want to be able to separate peaks separated by at least 1 Hz
    these_times=D.indsample(-20):D.indsample(0)-1;
    temp_data=D(match_str(D.chanlabels,{'Oz','POz'}),these_times,:); % D contains the data with channels * time * trials
    [logSNR, faxis, logpow]=get_logSNR(temp_data,D.fsample,param);

    %%% aggregate trials between probes
    for nbl=1:6
        these_probes=probe_res(probe_res(:,4)==nbl,:);
        these_trials=test_res(test_res(:,1)==nbl,:);
        
        for npr=1:10
            this_pr_tridx=these_probes(npr,6);
%             if npr==1
%                 last_pr_tridx=0;
%             else
%                 last_pr_tridx=these_probes(npr-1,6);
%             end  
            last_pr_tridx=this_pr_tridx-20;
            
            probe_details(1)=these_probes(npr,31); % look
            probe_details(2)=these_probes(npr,32); % look
            probe_details(3)=these_probes(npr,33); % look
            probe_details(4)=these_probes(npr,34); % look
            probe_details(5)=these_probes(npr,35); % look
            probe_details(6)=these_probes(npr,36); % look
            probe_details(7)=these_probes(npr,37); % look
            probe_details(8)=these_probes(npr,38); % look            
            
            temp_testres=these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,:);
%             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            temp_corr_go=nanmean(temp_testres(:,12));
            temp_corr_nogo=nanmean(temp_testres(:,11));
            temp_rt_go=nanmean(temp_testres(:,10)-temp_testres(:,8));
            % Colum order:
            
            temp_behav=[str2num(SubID) nbl these_probes(npr,5) this_pr_tridx probe_details temp_corr_go temp_corr_nogo temp_rt_go left_freq right_freq];
            
            % select probe
            this_trial_label=sprintf('B%g_P%g_',nbl,npr);
            idx_probe_eeg=find_trials(D.conditions,this_trial_label);
            if length(idx_probe_eeg)~=1
                warning('Problem')
                continue
            end
            
            FOI=[6 7.5 12 15 13.5];
            temp_SNR=[];
            temp_PWR=[];
            for nf=1:length(FOI)
                [~,fidx]=findclosest(faxis,FOI(nf));
                temp_SNR(:,nf)=squeeze(logSNR(:,fidx,idx_probe_eeg));
                temp_PWR(:,nf)=squeeze(logpow(:,fidx,idx_probe_eeg));
            end
            temp_EEGres=[temp_SNR(1,:) temp_SNR(2,:) temp_PWR(1,:) temp_PWR(2,:)];
            
            all_probes_headers={'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','CorrGo','CorrNoGo','RTGo','L_freq','R_freq',...
                'sOZ_f1','sOZ_f2','sOZ_2f1','sOZ_2f2','sOZ_IM',...
                'sPOZ_f1','sPOZ_f2','sPOZ_2f1','sPOZ_2f2','sPOZ_IM',...
                'pOZ_f1','pOZ_f2','pOZ_2f1','pOZ_2f2','pOZ_IM',...
                'pPOZ_f1','pPOZ_f2','pPOZ_2f1','pPOZ_2f2','pPOZ_IM'};
            all_probes_mat=[all_probes_mat ; [temp_behav temp_EEGres]];
        end
    end
end



%% transform into tables and export
% all_headers={'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe'};

tbl_probe=array2table(all_probes_mat,'VariableNames',all_probes_headers);
% 'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe'
tbl_probe.SubID=categorical(tbl_probe.SubID);
tbl_probe.Task=categorical(tbl_probe.Task);
tbl_probe.Look=categorical(tbl_probe.Look);
tbl_probe.State=categorical(tbl_probe.State);
tbl_probe.Orig=categorical(tbl_probe.Orig);
tbl_probe.L_freq=categorical(tbl_probe.L_freq);
tbl_probe.R_freq=categorical(tbl_probe.R_freq);

writetable(tbl_probe,[behav_path filesep 'WanderIM_ProbeResults_EEG.txt']);

