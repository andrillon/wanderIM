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
spec_path='/Volumes/ANDRILLON_HD1/BackUp_WanderIM_data_JULY2020/preproc_eeg/';
eeg_path=spec_path; %[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'trial_nfEEG_S3*.mat']);

%% loop across trials for baseline blocks
all_probes_mat=[];
all_probes_mat2=nan(6*10*length(bsl_files),142);
all_probes_STD=[];
all_probes_DEV=[];
all_probes_behav=[];
all_probes_P3amp=[];
all_probes_ERP=[];
all_probes_mat3=nan(64*6*10*length(bsl_files),19);
load('clusters_P3'); %,'time_W','cluster_pos_ID','cluster_neg_ID')
lastline=0;
lastline2=0;

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
            temp_testresgo=find(~isnan(these_trials2(:,12)));
            go_idx=these_trialsidx(temp_testresgo);
            go_idx=go_idx(end-17:end);
            temp_testresnogo=find(~isnan(these_trials2(:,11)));
            nogo_idx=these_trialsidx(temp_testresnogo);
            nogo_idx=nogo_idx(end-1:end);
            
            these_trialsidx=these_trialsidx(ismember(these_trialsidx,[go_idx ; nogo_idx]));
            temp_data2=temp_data(:,:,these_trialsidx);
            these_trials2=test_res(these_trialsidx,:);
            
            probe_details(1)=these_probes(npr,31); % look
            probe_details(2)=these_probes(npr,32); % look
            probe_details(3)=these_probes(npr,33); % look
            probe_details(4)=these_probes(npr,34); % look
            probe_details(5)=these_probes(npr,35); % look
            probe_details(6)=these_probes(npr,36); % look
            probe_details(7)=these_probes(npr,37); % look
            probe_details(8)=these_probes(npr,38); % look
            
            tcorr_go=these_trials2(~isnan(these_trials2(:,12)),:);
            tcorr_nogo=these_trials2(~isnan(these_trials2(:,11)),:);
            
            %             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            temp_corr_go=nanmean(tcorr_go(:,12));
            temp_corr_nogo=nanmean(tcorr_nogo(:,11));
            temp_rt_go=nanmean(tcorr_go(:,10)-tcorr_go(:,8));
            [dprime, crit]=calc_dprime(tcorr_go(:,12),tcorr_nogo(:,11)==0);
            % Colum order
            
            EEG_Gr=[];
            GoNoGo=~isnan(these_trials2(:,12));
            CorrUnCorr=nanmin(these_trials2(:,11:12),[],2);
            Behav_mat=[str2num(SubID) nbl npr these_probes(npr,5) probe_details mean(GoNoGo==0) mean(CorrUnCorr==0) NaN NaN];
            Behav_mat2=[repmat([str2num(SubID) nbl npr these_probes(npr,5) probe_details],length(GoNoGo),1) GoNoGo CorrUnCorr these_trialsidx these_trialsidx-this_pr_tridx ];
            
            % P300
            timeW_P300=[0.4 0.6];
            this_ch=match_str(D.chanlabels,'Pz');
            P3_amp=squeeze(nanmean(nanmean(temp_data2(this_ch,D.indsample(time_W(1)):D.indsample(time_W(2)),:),2),1))';
            
            %
            EEG_mat=P3_amp';
            
            temp_testresgo=find(~isnan(these_trials2(:,12)));
            this_ch=match_str(D.chanlabels,'Pz');
            
            temp_testresnogo=find(~isnan(these_trials2(:,11)));
            
            all_probes_STD=[all_probes_STD ; squeeze(nanmean(temp_data2(this_ch,:,temp_testresgo),3))];
            all_probes_DEV=[all_probes_DEV ; squeeze(nanmean(temp_data2(this_ch,:,temp_testresnogo),3))];
            all_probes_behav=[all_probes_behav ; [str2num(SubID) nbl npr these_probes(npr,5) probe_details temp_corr_go temp_corr_nogo dprime crit]];
            
            %             myChannelsIdx=match_str(D.chanlabels,{'Fz','Cz','Pz','Oz'});
            %             for nE=1:length(myChannelsIdx)
            %                 all_probes_ERP=[all_probes_ERP ; [Behav_mat2 nE*ones(size(Behav_mat2,1),1) squeeze(temp_data2(myChannelsIdx(nE),:,:))']];
            %             end
            
            tempP3=squeeze(nanmean(nanmean(temp_data2(this_ch,D.indsample(time_W(1)):D.indsample(time_W(2)),temp_testresgo),2),3))-squeeze(nanmean(nanmean(temp_data2(this_ch,D.indsample(time_W(1)):D.indsample(time_W(2)),temp_testresgo),2),3));
            all_probes_P3amp=[all_probes_P3amp ; tempP3];
            
            this_ch2=match_str(D.chanlabels,'Oz');
            if probe_details(2)~=4
                temp_ERP_STD{probe_details(2)}=[temp_ERP_STD{probe_details(2)} ; squeeze(temp_data2(this_ch2,:,temp_testresgo))'];
                temp_ERP_DEV{probe_details(2)}=[temp_ERP_DEV{probe_details(2)} ; squeeze(temp_data2(this_ch2,:,temp_testresnogo))'];
            end
            
            % P300
            timeW_P300=[0.4 0.6];
            P3_amp=squeeze(nanmean(temp_data(:,D.indsample(timeW_P300(1)):D.indsample(timeW_P300(2)),these_trialsidx),2))';
            P3_amp2=nanmean(P3_amp(GoNoGo==0,:))-nanmean(P3_amp(GoNoGo==1,:));
            
            % ERN
            timeW_ERN=[-0.2 0];
            temp_RTs=allRTs(these_trialsidx);
            ERN_amp=nan(length(temp_RTs),size(temp_data,1));
            for nRT=1:length(temp_RTs)
                if isnan(temp_RTs(nRT)) || temp_RTs(nRT)>0.9 || temp_RTs(nRT)<0.2
                    continue;
                else
                    ERN_amp(nRT,:)=squeeze(mean(temp_data(:,D.indsample(timeW_ERN(1)+temp_RTs(nRT)):D.indsample(timeW_ERN(2)+temp_RTs(nRT)),these_trialsidx(nRT)),2))';
                end
            end
            ERN_amp2=nanmean(ERN_amp(CorrUnCorr==0,:))-nanmean(ERN_amp(CorrUnCorr==1,:));
            
            EEG_mat2=[P3_amp2 ERN_amp2];
            Chan_headers=[];
            Chan_headers2=[];
            for nE=1:63
                Chan_headers{nE}=sprintf('P3_Ch%g',nE);
                Chan_headers2{nE}=sprintf('RN_Ch%g',nE);
            end
            all_probes_headers2=[{'SubID','nBlock','nProbe','Task','Look','State','Orig','Awa','Int','Eng','Perf','Vig','TrCond','CorrCond','nTrial','nTrial2'}, Chan_headers, Chan_headers2];
            all_probes_mat2(lastline+(1:size(Behav_mat,1)),:)=[Behav_mat EEG_mat2];
            lastline=lastline+size(Behav_mat,1);
            
            EEG_mat3=[0 0];
            EEGGR_mat3=[Behav_mat 0];
            EEG_mat3=[EEG_mat3; [P3_amp2' ERN_amp2']];
            EEGGR_mat3=[EEGGR_mat3 ; [repmat(Behav_mat,63,1) (1:63)']];
            all_probes_mat3(lastline2+(1:size(EEGGR_mat3,1)),:)=[EEGGR_mat3 EEG_mat3];
            lastline2=lastline2+size(EEGGR_mat3,1);
            all_probes_headers3={'SubID','nBlock','nProbe','Task','Look','State','Orig','Awa','Int','Eng','Perf','Vig','TrCond','CorrCond','nTrial','nTrial2','Chan','P3','ERN'};
        end
    end
    for nstate=1:3
        all_probes_STD2(n,nstate,:)=nanmean(temp_ERP_STD{nstate});
        all_probes_DEV2(n,nstate,:)=nanmean(temp_ERP_DEV{nstate});
        nall_probes_STD2(n,nstate)=size(temp_ERP_STD{nstate},1);
        nall_probes_DEV2(n,nstate)=size(temp_ERP_DEV{nstate},1);
    end
end

%%
% Behav_mat2=[repmat([str2num(SubID) nbl npr these_probes(npr,5) probe_details],length(GoNoGo),1) GoNoGo CorrUnCorr these_trialsidx these_trialsidx-this_pr_tridx ];
figure; format_fig; hold on;
% for ntask=1:2
%     subplot(2,1,ntask)
for nstate=1:3
    %         tempplot=all_probes_ERP(all_probes_ERP(:,size(Behav_mat2,2)+1)==2 & all_probes_ERP(:,4)==ntask & all_probes_ERP(:,7)==nstate & all_probes_ERP(:,13)==1,size(Behav_mat2,2)+2:end);
    templot=squeeze(all_probes_DEV2(all_probes_DEV2(:,nstate)>0,nstate,:));
    plot(D.time,nanmean(templot),'Color',Colors(nstate,:),'LineStyle','--')
    templot=squeeze(all_probes_STD2(all_probes_DEV2(:,nstate)>0,nstate,:));
    plot(D.time,nanmean(templot),'Color',Colors(nstate,:),'LineWidth',2)
end
% end
%%
tbl_probe=array2table(all_probes_mat2,'VariableNames',all_probes_headers2);
%'SubID','nBlock','nProbe','Task','Look','State','Orig','Awa','Int','Eng','Perf','Vig','TrCond','CorrCond','nTrial','nTrial2','P3','ERN','Chan'
tbl_probe(tbl_probe.State==4,:)=[];
tbl_probe.SubID=categorical(tbl_probe.SubID);
tbl_probe.Task=categorical(tbl_probe.Task);
tbl_probe.Look=categorical(tbl_probe.Look);
tbl_probe.State=categorical(tbl_probe.State);
tbl_probe.Orig=categorical(tbl_probe.Orig);
% tbl_probe.TrCond=categorical(tbl_probe.TrCond);
% tbl_probe.CorrCond=categorical(tbl_probe.CorrCond);
% tbl_probe.nTrial=categorical(tbl_probe.nTrial);
% tbl_probe.nTrial2=categorical(tbl_probe.nTrial2);

% tbl_probe.TrCond=reordercats(tbl_probe.TrCond,{'1','0'});
% tbl_probe.CorrCond=reordercats(tbl_probe.CorrCond,{'1','0'});
tbl_probe.State=reordercats(tbl_probe.State,{'1','2','3'});
tbl_probe.Task=reordercats(tbl_probe.Task,{'1','2'});
F_probe=tbl_probe(tbl_probe.Task=="1",:);
D_probe=tbl_probe(tbl_probe.Task=="2",:);

tbl_probe2=tbl_probe;
tbl_probe2.State=reordercats(tbl_probe2.State,{'2','1','3'});

%% Models with all E
tbl_allE=array2table(all_probes_mat3,'VariableNames',all_probes_headers3);
tbl_allE.P3(abs(tbl_allE.P3)>100)=nan;
tbl_allE.ERN(abs(tbl_allE.ERN)>100)=nan;
tbl_allE(tbl_allE.State==4,:)=[];

tbl_allE.SubID=categorical(tbl_allE.SubID);
tbl_allE.Task=categorical(tbl_allE.Task);
tbl_allE.Look=categorical(tbl_allE.Look);
tbl_allE.State=categorical(tbl_allE.State);
tbl_allE.Orig=categorical(tbl_allE.Orig);
% tbl_allE.TrCond=categorical(tbl_allE.TrCond);
% tbl_allE.CorrCond=categorical(tbl_allE.CorrCond);
% tbl_allE.nTrial=categorical(tbl_allE.nTrial);
% tbl_allE.nTrial2=categorical(tbl_allE.nTrial2);
tbl_allE.Chan=categorical(tbl_allE.Chan);

% tbl_allE.TrCond=reordercats(tbl_allE.TrCond,{'1','0'});
% tbl_allE.CorrCond=reordercats(tbl_allE.CorrCond,{'1','0'});
tbl_allE.State=reordercats(tbl_allE.State,{'1','2','3'});

%  mdl_P3_allEnotask= fitlme(tbl_allE,'P3~nBlock+Chan*State*TrCond+(1| SubID)');
%  mdl_P3_allEnoblock= fitlme(tbl_allE,'P3~Task*Chan*State*TrCond+(1| SubID)');

mdl_P3_allE= fitlme(tbl_allE,'P3~nBlock+Task*Chan*State+(1| SubID)');
mdl_ERN_allE= fitlme(tbl_allE,'ERN~nBlock+Task*Chan*State+(1| SubID)');

tbl_allE2=tbl_allE;
tbl_allE2.State=reordercats(tbl_allE2.State,{'2','1','3'});
mdl_P3_allE2= fitlme(tbl_allE2,'P3~nBlock+Task*Chan*State+(1| SubID)');
mdl_ERN_allE2= fitlme(tbl_allE2,'ERN~nBlock+Task*Chan*State+(1| SubID)');

%% Check main Effects of trial category
figure;
load(['EasyCap64_layout']);
stat_thr=0.05;
scaleVal=10;

this_model=mdl_P3_allE;
pVal=double(this_model.Coefficients(:,6));
beta=double(this_model.Coefficients(:,2));

subplot(2,2,1); format_fig;
ChanIdx=find_trials(this_model.CoefficientNames,'^Chan');
tempplot=beta(ChanIdx);
tempplotpV=pVal(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('P3')
colorbar;

subplot(2,2,3); format_fig;
ChanIdx=find_trials(this_model.CoefficientNames,'^Task_2:Chan');
tempplot=beta(ChanIdx);
tempplotpV=pVal(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('P3*Task')
colorbar;

this_model=mdl_ERN_allE;
pVal=double(this_model.Coefficients(:,6));
beta=double(this_model.Coefficients(:,2));

subplot(2,2,2); format_fig;
ChanIdx=find_trials(this_model.CoefficientNames,'^Chan');
tempplot=beta(ChanIdx);
tempplotpV=pVal(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('ERN')
colorbar;

subplot(2,2,4); format_fig;
ChanIdx=find_trials(this_model.CoefficientNames,'^Task_2:Chan');
tempplot=beta(ChanIdx);
tempplotpV=pVal(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.','k',10,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('ERN*Task')
colorbar;

%% Plot effect of mind-states
figure;
load(['EasyCap64_layout']);
stat_thr=0.05;
scaleVal=5;

this_model=mdl_P3_allE;
pVal=double(this_model.Coefficients(:,6));
beta=double(this_model.Coefficients(:,2));
this_model2=mdl_P3_allE2;
pVal2=double(this_model2.Coefficients(:,6));
beta2=double(this_model2.Coefficients(:,2));


subplot(2,3,1); format_fig;
ChanIdx=find_trials(this_model.CoefficientNames,'^State_2:Chan');
tempplot=beta(ChanIdx);
tempplotpV=pVal(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.',[1 1 1]*0.5,14,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('P3 MW vs ON')
colorbar;

subplot(2,3,2); format_fig;
ChanIdx=find_trials(this_model.CoefficientNames,'^State_3:Chan');
tempplot=beta(ChanIdx);
tempplotpV=pVal(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.',[1 1 1]*0.5,14,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('P3 MB vs ON')
colorbar;

subplot(2,3,3); format_fig;
ChanIdx=find_trials(this_model2.CoefficientNames,'^State_3:Chan');
tempplot=beta2(ChanIdx);
tempplotpV=pVal2(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.',[1 1 1]*0.5,14,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('P3 MB vs MW')
colorbar;


this_model=mdl_ERN_allE;
pVal=double(this_model.Coefficients(:,6));
beta=double(this_model.Coefficients(:,2));
this_model2=mdl_ERN_allE2;
pVal2=double(this_model2.Coefficients(:,6));
beta2=double(this_model2.Coefficients(:,2));


subplot(2,3,4); format_fig;
ChanIdx=find_trials(this_model.CoefficientNames,'^State_2:Chan');
tempplot=beta(ChanIdx);
tempplotpV=pVal(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.',[1 1 1]*0.5,14,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('ERN MW vs ON')
colorbar;

subplot(2,3,5); format_fig;
ChanIdx=find_trials(this_model.CoefficientNames,'^State_3:Chan');
tempplot=beta(ChanIdx);
tempplotpV=pVal(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.',[1 1 1]*0.5,14,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('ERN MB vs ON')
colorbar;

subplot(2,3,6); format_fig;
ChanIdx=find_trials(this_model2.CoefficientNames,'^State_3:Chan');
tempplot=beta2(ChanIdx);
tempplotpV=pVal2(ChanIdx);
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<stat_thr),'.',[1 1 1]*0.5,14,2},'whitebk','on');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
caxis([-1 1]*scaleVal)
title('ERN MB vs MW')
colorbar;

%%
clear P3_* ERN_*
fprintf('%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b%2.0f/63\n',nE)
    mdl= fitlme(tbl_probe,sprintf('P3_Ch%g~nBlock+Task*State*TrCond+(1| SubID)',nE));
    P3_tval(1:size(double(mdl.Coefficients(:,4)),1),nE)=double(mdl.Coefficients(:,4));
    P3_pval(1:size(double(mdl.Coefficients(:,4)),1),nE)=double(mdl.Coefficients(:,6));
    P3_varNames=mdl.CoefficientNames;
    
    mdl2= fitlme(tbl_probe2,sprintf('P3_Ch%g~nBlock+Task*State*TrCond+(1| SubID)',nE));
    P3_tval(size(double(mdl.Coefficients(:,4)),1)+(1:2),nE)=double(mdl2.Coefficients(match_str(mdl2.CoefficientNames,{'State_3:TrCond_0','Task_2:State_3:TrCond_0'}),4));
    P3_pval(size(double(mdl.Coefficients(:,4)),1)+(1:2),nE)=double(mdl2.Coefficients(match_str(mdl2.CoefficientNames,{'State_3:TrCond_0','Task_2:State_3:TrCond_0'}),6));
    P3_varNames=[P3_varNames {'State_3b:TrCond_0','Task_2:State_3b:TrCond_0'}];
    
    mdl= fitlme(tbl_probe,sprintf('RN_Ch%g~nBlock+Task*State*CorrCond+(1| SubID)',nE));
    ERN_tval(1:size(double(mdl.Coefficients(:,4)),1),nE)=double(mdl.Coefficients(:,4));
    ERN_pval(1:size(double(mdl.Coefficients(:,4)),1),nE)=double(mdl.Coefficients(:,6));
    ERN_varNames=mdl.CoefficientNames;
    
    mdl2= fitlme(tbl_probe2,sprintf('RN_Ch%g~nBlock+Task*State*CorrCond+(1| SubID)',nE));
    ERN_tval(size(double(mdl.Coefficients(:,4)),1)+(1:2),nE)=double(mdl2.Coefficients(match_str(mdl2.CoefficientNames,{'State_3:CorrCond_0','Task_2:State_3:CorrCond_0'}),4));
    ERN_pval(size(double(mdl.Coefficients(:,4)),1)+(1:2),nE)=double(mdl2.Coefficients(match_str(mdl2.CoefficientNames,{'State_3:CorrCond_0','Task_2:State_3:CorrCond_0'}),6));
    ERN_varNames=[ERN_varNames {'State_3b:CorrCond_0','Task_2:State_3b:CorrCond_0'}];
end

%%
load(['EasyCap64_layout'])
figure;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
varNames={'TrCond_0','State_2:TrCond_0','State_3:TrCond_0','State_3b:TrCond_0'};
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
varNames={'CorrCond_0','State_2:CorrCond_0','State_3:CorrCond_0','State_3b:CorrCond_0'};
for nplot=1:length(varNames)
    subplot(1,length(varNames),nplot)
    tempplot=ERN_tval(match_str(ERN_varNames,varNames{nplot}),:);
    tempplotpV=ERN_pval(match_str(ERN_varNames,varNames{nplot}),:);
    topoplot(tempplot, layout.chaninfo,'style','map','electrodes','off','emarker2',{find(tempplotpV<fdr(tempplotpV,0.05)),'.','k',10,2},'whitebk','on');
    caxis([-1 1]*5)
    title(varNames{nplot})
end
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
