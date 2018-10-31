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
all_probes_mat2=[];
all_probes_STD=[];
all_probes_DEV=[];
all_probes_behav=[];
all_probes_P3amp=[];
load('clusters_P3'); %,'time_W','cluster_pos_ID','cluster_neg_ID')

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
    if size(test_res,1)~=D.ntrials
        warning('pbme trial numbering')
    end
    
    left_freq(n)=SubjectInfo.FlickerL;
    right_freq(n)=SubjectInfo.FlickerR;
    
    these_times2=D.indsample(-0.2):D.indsample(1);
    these_timesbs=D.indsample(-0.2):D.indsample(0);
    %     temp_data=D(1:63,these_times2,:)-repmat(mean(D(1:63,these_timesbs,:),2),[1 length(these_times2) 1]);
    temp_data=D(match_str(D.chantype,'EEG'),:,:);
    temp_data=temp_data-repmat(mean(temp_data(:,these_timesbs,:),2),[1 size(temp_data,2) 1]);
    maxAmp=squeeze(max(abs(temp_data),[],2));
    badCh=mean(maxAmp>300,2)>0.25;
    fprintf('... ... %g bad channels found\n',sum(badCh))
    badTr=mean(maxAmp>300,1)>0.25;
    fprintf('... ... %g (%g %%) bad trials found\n',sum(badTr),mean(badTr))
    temp_data(badCh,:,:)=NaN;
    temp_data(:,:,badTr)=NaN;
    %     temp_data=temp_data-repmat(nanmean(temp_data,1),[size(temp_data,1) 1 1]);
    temp_data=temp_data-repmat(nanmean(temp_data(match_str(D.chanlabels,{'TP9','TP10'}),:,:),1),[size(temp_data,1) 1 1]);
    
    
    %%% aggregate trials between probes
    allRTs=test_res(:,10)-test_res(:,8);
    temp_ERP_STD=cell(1,3);
    temp_ERP_DEV=cell(1,3);
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
            these_trials2=test_res(these_trialsidx,:);
            temp_data2=temp_data(:,:,these_trialsidx);
            
            probe_details(1)=these_probes(npr,31); % look
            probe_details(2)=these_probes(npr,32); % look
            probe_details(3)=these_probes(npr,33); % look
            probe_details(4)=these_probes(npr,34); % look
            probe_details(5)=these_probes(npr,35); % look
            probe_details(6)=these_probes(npr,36); % look
            probe_details(7)=these_probes(npr,37); % look
            probe_details(8)=these_probes(npr,38); % look
            
            temp_testresgo=these_trials2(~isnan(these_trials2(:,12)),:);
            tcorr_go=(temp_testresgo(end-19:end,12))';%/corr_go(n,these_probes(npr,5));
            temp_testresnogo=these_trials2(~isnan(these_trials2(:,11)),:);
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
            GoNoGo=~isnan(these_trials2(:,12));
            CorrUnCorr=nanmin(these_trials2(:,11:12),[],2);
            Behav_mat=[repmat([str2num(SubID) nbl npr these_probes(npr,5) probe_details],length(GoNoGo),1) GoNoGo CorrUnCorr these_trialsidx these_trialsidx-this_pr_tridx ];
            
            % P300
            timeW_P300=[0.3 0.7];
%             P3_amp_pos=squeeze(nanmean(nanmean(temp_data2(elec_cluster_pos,D.indsample(time_W(1)):D.indsample(time_W(2)),:),2),1))';
%             P3_amp_neg=squeeze(nanmean(nanmean(temp_data2(elec_cluster_neg,D.indsample(time_W(1)):D.indsample(time_W(2)),:),2),1))';
%             P3_amp_diff=P3_amp_pos+P3_amp_neg;
            this_ch=match_str(D.chanlabels,'Pz');
            P3_amp=squeeze(nanmean(nanmean(temp_data2(this_ch,D.indsample(time_W(1)):D.indsample(time_W(2)),:),2),1))';
            
            %
            EEG_mat=P3_amp';
            
            all_probes_headers={'SubID','nBlock','nProbe','Task','Look','State','Orig','Awa','Int','Eng','Perf','Vig','TrCond','CorrCond','nTrial','nTrial2','P3amp'};
            all_probes_mat=[all_probes_mat ; [Behav_mat EEG_mat]];
            
            temp_testresgo=find(~isnan(these_trials2(:,12)));
            temp_testresgo=temp_testresgo(end-19:end);%/corr_go(n,these_probes(npr,5));
            this_ch=match_str(D.chanlabels,'Pz');
            
            temp_testresnogo=find(~isnan(these_trials2(:,11)));
            temp_testresnogo=temp_testresnogo(end-1:end);%/corr_go(n,these_probes(npr,5));
            
            all_probes_STD=[all_probes_STD ; squeeze(nanmean(temp_data2(this_ch,:,temp_testresgo),3))];
            all_probes_DEV=[all_probes_DEV ; squeeze(nanmean(temp_data2(this_ch,:,temp_testresnogo),3))];
            all_probes_behav=[all_probes_behav ; [str2num(SubID) nbl npr these_probes(npr,5) probe_details temp_corr_go temp_corr_nogo dprime crit]];
            
            tempP3=squeeze(nanmean(nanmean(temp_data2(this_ch,D.indsample(time_W(1)):D.indsample(time_W(2)),temp_testresgo),2),3))-squeeze(nanmean(nanmean(temp_data2(this_ch,D.indsample(time_W(1)):D.indsample(time_W(2)),temp_testresgo),2),3));
            all_probes_P3amp=[all_probes_P3amp ; tempP3];
            
            if probe_details(2)~=4
                temp_ERP_STD{probe_details(2)}=[temp_ERP_STD{probe_details(2)} ; squeeze(temp_data2(this_ch,:,temp_testresgo))'];
                temp_ERP_DEV{probe_details(2)}=[temp_ERP_DEV{probe_details(2)} ; squeeze(temp_data2(this_ch,:,temp_testresnogo))'];
            end
            
            % P300
            timeW_P300=[0.4 0.6];
            DEV_amp=squeeze(nanmean(temp_data2(:,D.indsample(timeW_P300(1)):D.indsample(timeW_P300(2)),temp_testresnogo),2))';
            STD_amp=squeeze(nanmean(temp_data2(:,D.indsample(timeW_P300(1)):D.indsample(timeW_P300(2)),temp_testresgo),2))';
            P3_amp=nanmean(DEV_amp)-nanmean(STD_amp);
            
            % ERN
            timeW_ERN=[-0.1 0];
            temp_RTs=allRTs(these_trialsidx);
            Resp_amp=nan(length(temp_RTs),size(temp_data2,1));
            for nRT=1:length(temp_RTs)
                if isnan(temp_RTs(nRT)) || temp_RTs(nRT)>0.9 || temp_RTs(nRT)<0.2
                    continue;
                else
                    Resp_amp(nRT,:)=squeeze(mean(temp_data2(:,D.indsample(timeW_ERN(1)+temp_RTs(nRT)):D.indsample(timeW_ERN(2)+temp_RTs(nRT)),nRT),2))';
                end
            end
            ERN_amp=nanmean(Resp_amp(temp_testresnogo,:))-nanmean(Resp_amp(temp_testresgo,:));
            
            EEG_mat2=[P3_amp' ERN_amp'];
            Behav_mat2=repmat([str2num(SubID) nbl npr these_probes(npr,5) probe_details],63,1);
            all_probes_headers2={'SubID','nBlock','nProbe','Task','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Chan','P3','ERN'};
            all_probes_mat2=[all_probes_mat2 ; [Behav_mat2 (1:63)' EEG_mat2]];
            
        end
    end
    for nstate=1:3
        all_probes_STD2(n,nstate,:)=nanmean(temp_ERP_STD{nstate});
        all_probes_DEV2(n,nstate,:)=nanmean(temp_ERP_DEV{nstate});
    end
end

%%
tbl_probe=array2table(all_probes_mat2,'VariableNames',all_probes_headers2);
%'SubID','nBlock','nProbe','Task','Look','State','Orig','Awa','Int','Eng','Perf','Vig','TrCond','CorrCond','nTrial','nTrial2','P3','ERN','Chan'
tbl_probe(tbl_probe.State==4,:)=[];
tbl_probe.SubID=categorical(tbl_probe.SubID);
tbl_probe.Task=categorical(tbl_probe.Task);
tbl_probe.Look=categorical(tbl_probe.Look);
tbl_probe.State=categorical(tbl_probe.State);
tbl_probe.Orig=categorical(tbl_probe.Orig);

tbl_probe.State=reordercats(tbl_probe.State,{'1','2','3'});
tbl_probe.Task=reordercats(tbl_probe.Task,{'1','2'});
F_probe=tbl_probe(tbl_probe.Task=="1",:);
D_probe=tbl_probe(tbl_probe.Task=="2",:);

tbl_probe2=tbl_probe;
tbl_probe2.State=reordercats(tbl_probe2.State,{'2','1','3'});

%%
clear P3_* ERN_*
fprintf('%2.0f/63\n',0)
for nE=1:63
       fprintf('\b\b\b\b\b\b%2.0f/63\n',nE)
       subtbl_probe=tbl_probe(tbl_probe.Chan==nE,:);
       mdl= fitlme(subtbl_probe,'P3~nBlock+Task*State+(1| SubID)');
    P3_tval(1:size(double(mdl.Coefficients(:,4)),1),nE)=double(mdl.Coefficients(:,4));
    P3_pval(1:size(double(mdl.Coefficients(:,4)),1),nE)=double(mdl.Coefficients(:,6));
    P3_varNames=mdl.CoefficientNames;
    
  subtbl_probe=tbl_probe2(tbl_probe2.Chan==nE,:);
       mdl2= fitlme(subtbl_probe,'P3~nBlock+Task*State+(1| SubID)');
       P3_tval(size(double(mdl.Coefficients(:,4)),1)+(1:2),nE)=double(mdl2.Coefficients(match_str(mdl2.CoefficientNames,{'State_3','Task_2:State_3'}),4));
    P3_pval(size(double(mdl.Coefficients(:,4)),1)+(1:2),nE)=double(mdl2.Coefficients(match_str(mdl2.CoefficientNames,{'State_3','Task_2:State_3'}),6));
    P3_varNames=[P3_varNames {'State_3b','Task_2:State_3b'}];
    
  subtbl_probe=tbl_probe(tbl_probe.Chan==nE,:);
       mdl= fitlme(subtbl_probe,'ERN~nBlock+Task*State+(1| SubID)');
       ERN_tval(1:size(double(mdl.Coefficients(:,4)),1),nE)=double(mdl.Coefficients(:,4));
    ERN_pval(1:size(double(mdl.Coefficients(:,4)),1),nE)=double(mdl.Coefficients(:,6));
    ERN_varNames=mdl.CoefficientNames;
    
subtbl_probe=tbl_probe2(tbl_probe2.Chan==nE,:);
       mdl2= fitlme(subtbl_probe,'ERN~nBlock+Task*State+(1| SubID)');
       ERN_tval(size(double(mdl.Coefficients(:,4)),1)+(1:2),nE)=double(mdl2.Coefficients(match_str(mdl2.CoefficientNames,{'State_3','Task_2:State_3'}),4));
    ERN_pval(size(double(mdl.Coefficients(:,4)),1)+(1:2),nE)=double(mdl2.Coefficients(match_str(mdl2.CoefficientNames,{'State_3','Task_2:State_3'}),6));
    ERN_varNames=[ERN_varNames {'State_3b','Task_2:State_3b'}];
end
% 
%%
load(['EasyCap64_layout'])
figure;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
varNames={'Task_2','State_2','State_3','State_3b'};
for nplot=1:length(varNames)
    subplot(1,length(varNames),nplot)
    tempplot=P3_tval(match_str(P3_varNames,varNames{nplot}),:);
    tempplotpV=P3_pval(match_str(P3_varNames,varNames{nplot}),:);
    topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<fdr(tempplotpV,0.05)),'.','k',10,2},'whitebk','on');
    caxis([-1 1]*5)
    title(varNames{nplot})
end
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));



figure;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
for nplot=1:length(varNames)
    subplot(1,length(varNames),nplot)
    tempplot=ERN_tval(match_str(ERN_varNames,varNames{nplot}),:);
    tempplotpV=ERN_pval(match_str(ERN_varNames,varNames{nplot}),:);
    topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<fdr(tempplotpV,0.05)),'.','k',10,2},'whitebk','on');
 caxis([-1 1]*5)
    title(varNames{nplot})
end
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));


