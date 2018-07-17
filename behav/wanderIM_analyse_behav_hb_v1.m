%%
clear all
close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))
addpath(genpath(spm12_path))

data_path=[root_path filesep 'behav/'];
preproc_path=[root_path filesep 'preproc_eeg/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

%%
all_test_res=[];
all_probes_mat=[];
for n=1:length(files)
    %% load files
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
    hb_filename=['hb_detection_S' SubID];
    trig_filename=['triggers_S' SubID];
    load([preproc_path filesep hb_filename])
    load([preproc_path filesep trig_filename])

    %% Gather behav
    % load
     temp_perf=min(test_res(:,11:12),[],2);
    temp_cat=(test_res(:,5)==test_res(:,6));
    temp_RT=(test_res(:,10)-test_res(:,8));
    all_test_res=[all_test_res ; [str2num(SubID)*ones(size(test_res,1),1) test_res(:,[1 2 4 5]) temp_perf temp_RT temp_cat]];
    all_test_headers={'SubID','nBlock','Task','nTrial','StimID','Corr','RT','TrCat'};

    % probes
    for nbl=1:6
        these_probes=probe_res(probe_res(:,4)==nbl,:);
        these_trials=test_res(test_res(:,1)==nbl,:);
        for npr=1:10
            this_pr_tridx=these_probes(npr,6);
            this_pr_tridx2=start_probe(npr);
            if npr==1
                last_pr_tridx=0;
                last_pr_tridx2=clean_start_trial(1);
            else
                last_pr_tridx=these_probes(npr-1,6);
                last_pr_tridx2=start_probe(npr-1);
            end            
            probe_details(1)=these_probes(npr,31); % look
            probe_details(2)=these_probes(npr,32); % look
            probe_details(3)=these_probes(npr,33); % look
            probe_details(4)=these_probes(npr,34); % look
            probe_details(5)=these_probes(npr,35); % look
            probe_details(6)=these_probes(npr,36); % look
            probe_details(7)=these_probes(npr,37); % look
            probe_details(8)=these_probes(npr,38); % look            
            
            temp_testres=these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,:);
            tcorr_perf_go=nanmean(temp_testres(:,12));
            tcorr_perf_nogo=nanmean(temp_testres(:,11));
            tcorr_RT=nanmean(temp_testres(:,10)-temp_testres(:,8));
             
            %% Gather HB relative to probes
            this_hbrate=sum(hb_times>last_pr_tridx2 & hb_times<this_pr_tridx2)/((this_pr_tridx2-last_pr_tridx2+1)/500/60);
            this_hbrate2=sum(hb_times>this_pr_tridx2-20*500 & hb_times<this_pr_tridx2)/((this_pr_tridx2-this_pr_tridx2+20*500+1)/500/60);
            
            % Colum order:
            all_probes_headers={'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Go','NoGo','RT','HB','HB2'};
            all_probes_mat=[all_probes_mat ; [str2num(SubID) nbl these_probes(npr,5) this_pr_tridx probe_details tcorr_perf_go tcorr_perf_nogo tcorr_RT this_hbrate this_hbrate2]];
           
        end
    end
end  
    

%%
tbl_probe=array2table(all_probes_mat,'VariableNames',all_probes_headers);
% 'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe'
tbl_probe.SubID=categorical(tbl_probe.SubID);
tbl_probe.Task=categorical(tbl_probe.Task);
tbl_probe.Look=categorical(tbl_probe.Look);
tbl_probe.State=categorical(tbl_probe.State);
tbl_probe.Orig=categorical(tbl_probe.Orig);

writetable(tbl_probe,[data_path filesep 'WanderIM_ProbeResults_HB.txt']);

lme_0= fitlme(tbl_probe,'HB~1+(1|SubID)');
lme_full= fitlme(tbl_probe,'RT~Vig+(1|SubID)');
