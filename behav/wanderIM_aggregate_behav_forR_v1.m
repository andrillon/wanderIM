%%
clear all
close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))

% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];
% load([pwd filesep '..' filesep 'paper' filesep 'paper_SubID'])

% load([data_path filesep 'CARS_quest'])
%%
all_test_res=[];
all_probes_mat=[];
all_probes_mat2=[];
RTs=[];
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
       probe_res(probe_res(:,32)==4,32)=3;
 SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
     if ~ismember(SubID,GoodSubID)
        continue;
    end
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
    
    %%% aggregate trials between probes
    countpr=0;
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
            % Colum order:
            all_probes_headers={'SubID','nBlock','nProbe','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe'};
            all_probes_mat=[all_probes_mat ; [repmat([str2num(SubID) nbl countpr these_probes(npr,5) this_pr_tridx probe_details],size(temp_testres,1),1) tcorr_perf tcorr_RT tcorr_catTr tcorr_distPr']];
            all_probes_mat2=[all_probes_mat2 ; [str2num(SubID) nbl countpr these_probes(npr,5) this_pr_tridx probe_details nanmean(tcorr_perf) nanmean(tcorr_RT) nanmean(tcorr_catTr) ]];
        end
    end
end

%% transform into tables and export
all_test_res(:,end-1)=all_test_res(:,end-1)==0;
tbl_test=array2table(all_test_res,'VariableNames',all_test_headers);
% 'SubID','nBlock','Task','nTrial','StimID','Corr','RT','TrCat'
tbl_test.SubID=categorical(tbl_test.SubID);
tbl_test.Task=categorical(tbl_test.Task);
tbl_test.StimID=categorical(tbl_test.StimID);
tbl_test.TrCat=categorical(tbl_test.TrCat);

tbl_test.cond_v=cell(length(tbl_test.TrCat),1);
tbl_test.cond_v(tbl_test.TrCat=="1" & tbl_test.Task=="1")={'goF'};
tbl_test.cond_v(tbl_test.TrCat=="0" & tbl_test.Task=="1")={'nogoF'};
tbl_test.cond_v(tbl_test.TrCat=="1" & tbl_test.Task=="2")={'goD'};
tbl_test.cond_v(tbl_test.TrCat=="0" & tbl_test.Task=="2")={'nogoD'};
writetable(tbl_test,[data_path filesep 'WanderIM_TestResults2.txt']);


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


writetable(tbl_probe,[data_path filesep 'WanderIM_ProbeResults2_MW_new.txt']);

uniSub=(unique(tbl_probe.SubID));
for nS=1:length(uniSub)
    for nSt=1:3
        numberState(nS,nSt)=length(unique(tbl_probe.nProbe(tbl_probe.State==num2str(nSt) &tbl_probe.SubID==uniSub(nS))));
    end
    for nVg=1:4
        numberVig(nS,nVg)=length(unique(tbl_probe.nProbe(tbl_probe.Vig==(nVg) & tbl_probe.SubID==uniSub(nS))));
    end
end
%% Examples of models (to do in R)
lme_0= fitlme(tbl_probe,'Corr~1+(1|SubID)');
lme_1= fitlme(tbl_probe,'Corr~1+(1|SubID)');

lme_full= fitlme(tbl_probe,'Corr~Task+nBlock+nTrial+TrCat+DistProbe+State+Task*State+(1|SubID)');
lme_full2= fitlme(tbl_probe,'Corr~Task+nBlock+nTrial+TrCat+DistProbe+State+DistProbe*State+Task*State+(1|SubID)');

[v] = compare(lme_full,lme_full2);
if v.pValue(end)<0.05
    fprintf('More complex model (%s) wins (pV=%1.6f)\n',lme_full2.Formula,v.pValue(end))
else
    fprintf('Less complex model (%s) wins (pV=%1.6f)\n',lme_full.Formula,v.pValue(end))
end