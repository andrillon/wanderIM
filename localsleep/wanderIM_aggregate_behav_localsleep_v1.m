%%
clear all
close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))

% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);
eeg_path=[root_path filesep 'preproc_eeg'];

state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];


prticle_Thr=80; % 80 or 90 or 95
LimFrqW=[1 4]; % [1 4] or [4 10]
AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
thisChannel=17;

% load([data_path filesep 'CARS_quest'])
%%
all_test_res=[];
all_probes_mat=[];
RTs=[];
for n=1:20 %length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
%     CARS_flag(n)=CARS_bool(CARS_bool(:,1)==str2num(SubID),2);
    % SART
    %  1: num block
    %  2: block cond (1: Faces / 2: Squares)
    %  3: image set
    %  4: num trial
    %  5: seq trial
    %  6: target
    %  7: resp
    %  8: stim onset
    %  9: stim pre
    % 10: resp onset
    % 11: nogo
    % 12: go
    
    %%% cleaning too fast RTs
    RTs=[RTs; test_res(:,10)-test_res(:,8)];
%     findFalseStarts=find(test_res(:,10)-test_res(:,8)<0.1);
%     test_res(findFalseStarts,[10 11 12])=NaN;
%     warning('correcting for false starts')
    for nbt=1:2
        tp_nogos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,11)),11);
        tp_gos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,12)),12);
        [dprime_test(n,nbt), crit_test(n,nbt)]=calc_dprime((tp_gos==1),(tp_nogos==0));
        corr_go(n,nbt)=nanmean(tp_gos);
        corr_nogo(n,nbt)=nanmean(tp_nogos);
    end
    temp_perf=min(test_res(:,11:12),[],2);    
      temp_cat=(test_res(:,5)==test_res(:,6));
  code_resp=nan(length(temp_perf),1);
    code_resp(temp_perf==1 & temp_cat==0)=1;
    code_resp(temp_perf==1 & temp_cat==1)=0;
    code_resp(temp_perf==0 & temp_cat==0)=0;
    code_resp(temp_perf==0 & temp_cat==1)=1;
    temp_RT=(test_res(:,10)-test_res(:,8));
    all_test_res=[all_test_res ; [str2num(SubID)*ones(size(test_res,1),1) test_res(:,[1 2 4 5]) temp_perf temp_RT temp_cat code_resp]];
    all_test_headers={'SubID','nBlock','Task','nTrial','StimID','Corr','RT','TrCat','Resp'};

    % probe
    %  1: num probe within block
    %  2: probe onset (expected)
    %  3: probe onset (actual)
    %  4: num block
    %  5: block type (1: FACE / 2: SQUARE)
    %  6: num trial
    %  7: resp Q1 Looking (1: Yes / 2: No)
    %  8: resp Q2 Mind-State (1: On / 2: MW / 3: MB / 4: ?)
    %  9: resp Q3 Origine (1: room / 2: personal / 3: task)
    % 10: resp Q4 Awareness (1 (fully) to 4 (not at all))
    % 11: resp Q5 Intention (1 (fully) to 4 (not at all))
    % 12: resp Q6 Engagement (1 (not) to 4 (very))
    % 13: resp Q7 Performance (1 (bad) to 4 (good))
    % 14: resp Q8 Alterness (1 (alert) to 4 (sleepy))
    % 15-22: onset R1-8
    % 23-30: onset Q1-8
    % 31-38: resp number Q1-8
    for nbt=1:2
        temp=probe_res(probe_res(:,5)==nbt,[5 31:38]);
        clear pption pption2
        for nstate=1:4
            pption(nstate)=mean(temp(:,3)==nstate);
            pption2(nstate)=mean(temp(temp(:,3)==2,4)==nstate);
        end
        pption_MS(n,nbt,:)=pption;
        pption_Ori(n,nbt,:)=pption2;
    end
    for nbt=1:2
        tp_probes=probe_res(probe_res(:,5)==nbt,31:38);
        mean_bytask_prf(n,nbt)=nanmean(tp_probes(:,7));
        for nstate=1:3
            mean_count_mindstates(n,nbt,nstate)=sum(tp_probes(:,2)==(nstate));
            mean_count_ori(n,nbt,nstate)=sum(tp_probes(:,3)==(nstate));
            
            if sum(tp_probes(:,2)==(nstate))~=0

                mean_byprobe_awa(n,nbt,nstate)=5-nanmean(tp_probes(tp_probes(:,2)==(nstate),4));
                mean_byprobe_wil(n,nbt,nstate)=5-nanmean(tp_probes(tp_probes(:,2)==(nstate),5));
                mean_byprobe_eng(n,nbt,nstate)=nanmean(tp_probes(tp_probes(:,2)==(nstate),6));
                mean_byprobe_prf(n,nbt,nstate)=nanmean(tp_probes(tp_probes(:,2)==(nstate),7));
                mean_byprobe_vig(n,nbt,nstate)=5-nanmean(tp_probes(tp_probes(:,2)==(nstate),8));
            else
                mean_byprobe_awa(n,nbt,nstate)=nan;
                mean_byprobe_wil(n,nbt,nstate)=nan;
                mean_byprobe_eng(n,nbt,nstate)=nan;
                mean_byprobe_prf(n,nbt,nstate)=nan;
                mean_byprobe_vig(n,nbt,nstate)=nan;
            end
        end
    end
    
    load([eeg_path filesep 'wanderIM_twa3_' SubID])
    
    thisE_Waves=all_Waves(all_Waves(:,3)==thisChannel,:);
    temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./500;
    temp_freq=1./temp_len;
    temp_abs=1./temp_len;
    temp_p2p=thisE_Waves(:,AmpCriterionIdx);
    thr_Wave2(n)=prctile(all_Waves((temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),4),prticle_Thr);
    thisE_Waves_Onsets=thisE_Waves(temp_p2p>thr_Wave2(n),[2 5]);
    
     
    %%% aggregate trials between probes
    countpr=0;
    probe_res=[probe_res (probe_res(:,4)-1)*10+probe_res(:,1)];
    for nbl=1:6
        these_probes=probe_res(probe_res(:,4)==nbl,:);
        these_trials=test_res(test_res(:,1)==nbl,:);
        
        for npr=1:10
            countpr=countpr+1;
            this_pr_tridx=these_probes(npr,6);
            if npr==1
                last_pr_tridx=0;
            else
                last_pr_tridx=these_probes(npr-1,6);
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
%             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            tcorr_perf=min(temp_testres(:,11:12),[],2);
            tcorr_catTr=(temp_testres(:,5)==temp_testres(:,6)); %0 for go, 1 for nogo
            tcorr_RT=temp_testres(:,10)-temp_testres(:,8);
            tcorr_distPr=nan(1,length(tcorr_catTr));
            tcorr_distPr(tcorr_catTr==1)=-sum(tcorr_catTr==1):1:-1;
            tcorr_distPr(tcorr_catTr==0)=-sum(tcorr_catTr==0):1:-1;
            
            % Waves
            thisProbe_WavesOnset=thisE_Waves_Onsets(thisE_Waves_Onsets(:,1)==((nbl-1)*10+npr),2)-20*500;
            theseTrials_Onset=round((temp_testres(:,8)-these_probes(npr,3))*500);
            numAssociatedWaves=zeros(length(theseTrials_Onset),1);
            if ~isempty(thisProbe_WavesOnset)
                for nw=1:length(thisProbe_WavesOnset)
                    thisw=find(theseTrials_Onset<=thisProbe_WavesOnset(nw) & [theseTrials_Onset(2:end) ; 0]>thisProbe_WavesOnset(nw));
                    numAssociatedWaves(thisw)=numAssociatedWaves(thisw)+1;
%                     if numAssociatedWaves(thisw)>2
%                         pause;
%                     end
                end
            end
            % Colum order:
            all_probes_headers={'SubID','nBlock','nProbe','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe','nWave'};
            all_probes_mat=[all_probes_mat ; [repmat([str2num(SubID) nbl countpr these_probes(npr,5) this_pr_tridx probe_details],size(temp_testres,1),1) tcorr_perf tcorr_RT tcorr_catTr tcorr_distPr' numAssociatedWaves]];
        end
    end
end

%% transform into tables and export
tbl_probe=array2table(all_probes_mat,'VariableNames',all_probes_headers);
% 'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe'
tbl_probe.SubID=categorical(tbl_probe.SubID);
tbl_probe.Task=categorical(tbl_probe.Task);
tbl_probe.Look=categorical(tbl_probe.Look);
tbl_probe.State=categorical(tbl_probe.State);
tbl_probe.Orig=categorical(tbl_probe.Orig);
tbl_probe.TrCat=categorical(tbl_probe.TrCat);


tbl_probe.response=nan(length(tbl_probe.TrCat),1);
tbl_probe.response(tbl_probe.TrCat=="1" & tbl_probe.Corr==1)=0;
tbl_probe.response(tbl_probe.TrCat=="1" & tbl_probe.Corr==0)=1;
tbl_probe.response(tbl_probe.TrCat=="0" & tbl_probe.Corr==1)=1;
tbl_probe.response(tbl_probe.TrCat=="0" & tbl_probe.Corr==0)=0;

tbl_probe.stimulus=nan(length(tbl_probe.TrCat),1);
tbl_probe.stimulus=tbl_probe.TrCat=="0";

tbl_probe.vigcat=tbl_probe.Vig>=3;


tbl_probe.cond_v=cell(length(tbl_probe.TrCat),1);
tbl_probe.cond_v(tbl_probe.TrCat=="1" & tbl_probe.State=="1")={'nogo_on'};
tbl_probe.cond_v(tbl_probe.TrCat=="1" & tbl_probe.State=="2")={'nogo_mw'};
tbl_probe.cond_v(tbl_probe.TrCat=="1" & tbl_probe.State=="3")={'nogo_mb'};
tbl_probe.cond_v(tbl_probe.TrCat=="0" & tbl_probe.State=="1")={'go_on'};
tbl_probe.cond_v(tbl_probe.TrCat=="0" & tbl_probe.State=="2")={'go_mw'};
tbl_probe.cond_v(tbl_probe.TrCat=="0" & tbl_probe.State=="3")={'go_mb'};


writetable(tbl_probe,[data_path filesep 'WanderIM_ProbeResults2_thWaves.txt']);

%%
tbl_probe=array2table(all_probes_mat,'VariableNames',all_probes_headers);
% 'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe'
tbl_probe.SubID=categorical(tbl_probe.SubID);
tbl_probe.Task=categorical(tbl_probe.Task);
tbl_probe.Look=categorical(tbl_probe.Look);
tbl_probe.State=categorical(tbl_probe.State);
tbl_probe.Orig=categorical(tbl_probe.Orig);
tbl_probe.TrCat=categorical(tbl_probe.TrCat);


tbl_probe.response=nan(length(tbl_probe.TrCat),1);
tbl_probe.response(tbl_probe.TrCat=="1" & tbl_probe.Corr==1)=0;
tbl_probe.response(tbl_probe.TrCat=="1" & tbl_probe.Corr==0)=1;
tbl_probe.response(tbl_probe.TrCat=="0" & tbl_probe.Corr==1)=1;
tbl_probe.response(tbl_probe.TrCat=="0" & tbl_probe.Corr==0)=0;

tbl_probe.stimulus=nan(length(tbl_probe.TrCat),1);
tbl_probe.stimulus=tbl_probe.TrCat=="0";

tbl_probe.vigcat=tbl_probe.nWave>0;


tbl_probe.cond_v=cell(length(tbl_probe.TrCat),1);
tbl_probe.cond_v(tbl_probe.TrCat=="1" & tbl_probe.nWave==0)={'nogo_w-'};
tbl_probe.cond_v(tbl_probe.TrCat=="1" & tbl_probe.nWave>0)={'nogo_w+'};
tbl_probe.cond_v(tbl_probe.TrCat=="0" & tbl_probe.nWave==0)={'go_w-'};
tbl_probe.cond_v(tbl_probe.TrCat=="0" & tbl_probe.nWave>0)={'go_w+'};


writetable(tbl_probe,[data_path filesep 'WanderIM_ProbeResults2_thWaves_vig.txt']);

mysubs=unique(tbl_probe.SubID);
for ns=1:length(mysubs)
    mean_pp_waves(ns)=mean(tbl_probe.nWave(tbl_probe.SubID==mysubs(ns))>0);
    mean_sleepiness(ns)=mean(tbl_probe.Vig(tbl_probe.SubID==mysubs(ns)));
end

%%
mdl1=fitlme(tbl_probe,'RT~Task*nWave+(1|SubID)');
mdl2=fitlme(tbl_probe,'RT~Task*Vig+(1|SubID)');
mdl3=fitlme(tbl_probe,'RT~Task*nWave*State+Vig+(1|SubID)');