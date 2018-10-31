%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

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
bsl_files=dir([eeg_path filesep 'trial_nfEEG_S3*.mat']);

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
    
    
    left_freq(n)=SubjectInfo.FlickerL;
    right_freq(n)=SubjectInfo.FlickerR;
    
    % sanity check: computer ERP on block onset
    these_timesbs=D.indsample(-0.2):D.indsample(0);
    %     temp_data=D(1:63,these_times2,:)-repmat(mean(D(1:63,these_timesbs,:),2),[1 length(these_times2) 1]);
    temp_data=D(match_str(D.chantype,'EEG'),:,:)-repmat(mean(D(match_str(D.chantype,'EEG'),:,:),1),length(match_str(D.chantype,'EEG')),1);
    
    %%% aggregate trials between probes
    allRTs=test_res(:,10)-test_res(:,8);
    for nbl=1:6
        these_probes=probe_res(probe_res(:,4)==nbl,:);
        these_trials=test_res(test_res(:,1)==nbl,:);
        for npr=1:10
            this_pr_tridx=these_probes(npr,6);
            if npr==1
                last_pr_tridx=0;
            else
                last_pr_tridx=these_probes(npr-1,6);
            end
            these_trialsidx=find(test_res(:,1)==nbl & (test_res(:,4)>last_pr_tridx & test_res(:,4)<this_pr_tridx));
            
            probe_details(1)=these_probes(npr,31); % look
            probe_details(2)=these_probes(npr,32); % look
            probe_details(3)=these_probes(npr,33); % look
            probe_details(4)=these_probes(npr,34); % look
            probe_details(5)=these_probes(npr,35); % look
            probe_details(6)=these_probes(npr,36); % look
            probe_details(7)=these_probes(npr,37); % look
            probe_details(8)=these_probes(npr,38); % look
            
            temp_testresgo=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,12)),:);
            tcorr_go=(temp_testresgo(end-19:end,12))';%/corr_go(n,these_probes(npr,5));
            temp_testresnogo=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,11)),:);
            tcorr_nogo=(temp_testresnogo(end-1:end,11))';%/corr_go(n,these_probes(npr,5));
            
            %             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            temp_corr_go=nanmean(tcorr_go);
            temp_corr_nogo=nanmean(tcorr_nogo);
            temp_rt_go=nanmean(temp_testresgo(:,10)-temp_testresgo(:,8));
            temp_corr_go2=tcorr_go;
            temp_corr_nogo2=tcorr_nogo;
            [dprime, crit]=calc_dprime(temp_corr_go2,temp_corr_nogo2==0);
            % Colum order
            
            EEG_Gr=[];
            EEG_mat=[];
            GoNoGo=~isnan(test_res(these_trialsidx,12));
            CorrUnCorr=nanmin(test_res(these_trialsidx,11:12),[],2);
               Behav_mat=repmat([repmat([str2num(SubID) nbl npr these_probes(npr,5) probe_details],length(GoNoGo),1) GoNoGo CorrUnCorr these_trialsidx these_trialsidx-this_pr_tridx ],size(temp_data,1),1);
         
            for nE=1:size(temp_data,1)
                % P300
                timeW_P300=[0.4 0.6];
                P3_amp=squeeze(mean(temp_data(nE,D.indsample(timeW_P300(1)):D.indsample(timeW_P300(2)),these_trialsidx),2));
                
                % ERN
                timeW_ERN=[-0.1 0];
                temp_RTs=allRTs(these_trialsidx);
                ERN_amp=nan(length(temp_RTs),1);
                for nRT=1:length(temp_RTs)
                    if isnan(temp_RTs(nRT)) %|| temp_RTs(nRT)>1
                        continue;
                    else
                        ERN_amp(nRT,1)=squeeze(mean(temp_data(nE,D.indsample(timeW_ERN(1)+temp_RTs(nRT)):D.indsample(timeW_ERN(2)+temp_RTs(nRT)),these_trialsidx(nRT)),2));
                    end
                end
                
                EEG_Gr=[EEG_Gr ; nE*ones(length(these_trialsidx),1)];
                EEG_mat=[EEG_mat ; [P3_amp ERN_amp]];
                
                                
            end
            
            all_probes_headers={'SubID','nBlock','nProbe','Task','Look','State','Orig','Awa','Int','Eng','Perf','Vig','TrCond','CorrCond','nTrial','nTrial2','P3','ERN','Chan'};
            all_probes_mat=[all_probes_mat ; [Behav_mat EEG_mat EEG_Gr]];
            
        end
    end
end

%%
tbl_probe=array2table(all_probes_mat,'VariableNames',all_probes_headers);
%'SubID','nBlock','nProbe','Task','Look','State','Orig','Awa','Int','Eng','Perf','Vig','TrCond','CorrCond','nTrial','nTrial2','P3','ERN','Chan'
tbl_probe.SubID=categorical(tbl_probe.SubID);
tbl_probe.Task=categorical(tbl_probe.Task);
tbl_probe.Look=categorical(tbl_probe.Look);
tbl_probe.State=categorical(tbl_probe.State);
tbl_probe.Orig=categorical(tbl_probe.Orig);
tbl_probe.TrCond=categorical(tbl_probe.TrCond);
tbl_probe.CorrCond=categorical(tbl_probe.CorrCond);
tbl_probe.nTrial=categorical(tbl_probe.nTrial);
tbl_probe.nTrial2=categorical(tbl_probe.nTrial2);
tbl_probe.Chan=categorical(tbl_probe.Chan);
