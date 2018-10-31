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
            Behav_mat=[repmat([str2num(SubID) nbl npr these_probes(npr,5) probe_details],length(GoNoGo),1) GoNoGo CorrUnCorr these_trialsidx these_trialsidx-this_pr_tridx ];
            
            % P300
            timeW_P300=[0.4 0.7];
            P3_amp=squeeze(nanmean(temp_data(:,D.indsample(timeW_P300(1)):D.indsample(timeW_P300(2)),these_trialsidx),2))';
            
            % ERN
            timeW_ERN=[-0.1 0];
            temp_RTs=allRTs(these_trialsidx);
            ERN_amp=nan(length(temp_RTs),size(temp_data,1));
%             for nRT=1:length(temp_RTs)
%                 if isnan(temp_RTs(nRT)) %|| temp_RTs(nRT)>1
%                     continue;
%                 else
%                     ERN_amp(nRT,:)=squeeze(mean(temp_data(:,D.indsample(timeW_ERN(1)+temp_RTs(nRT)):D.indsample(timeW_ERN(2)+temp_RTs(nRT)),these_trialsidx(nRT)),2))';
%                 end
%             end
%             
            EEG_mat=[EEG_mat ; [P3_amp ERN_amp]];
            Chan_headers=[];
            Chan_headers2=[];
            for nE=1:63
                Chan_headers{nE}=sprintf('P3_Ch%g',nE);
                Chan_headers2{nE}=sprintf('RN_Ch%g',nE);
            end
            all_probes_headers=[{'SubID','nBlock','nProbe','Task','Look','State','Orig','Awa','Int','Eng','Perf','Vig','TrCond','CorrCond','nTrial','nTrial2'}, Chan_headers, Chan_headers2];
            all_probes_mat=[all_probes_mat ; [Behav_mat EEG_mat]];
            
        end
    end
end

%%
tbl_probe=array2table(all_probes_mat,'VariableNames',all_probes_headers);
%'SubID','nBlock','nProbe','Task','Look','State','Orig','Awa','Int','Eng','Perf','Vig','TrCond','CorrCond','nTrial','nTrial2','P3','ERN','Chan'
tbl_probe(tbl_probe.State==4,:)=[];
tbl_probe.SubID=categorical(tbl_probe.SubID);
tbl_probe.Task=categorical(tbl_probe.Task);
tbl_probe.Look=categorical(tbl_probe.Look);
tbl_probe.State=categorical(tbl_probe.State);
tbl_probe.Orig=categorical(tbl_probe.Orig);
tbl_probe.TrCond=categorical(tbl_probe.TrCond);
tbl_probe.CorrCond=categorical(tbl_probe.CorrCond);
tbl_probe.nTrial=categorical(tbl_probe.nTrial);
tbl_probe.nTrial2=categorical(tbl_probe.nTrial2);

tbl_probe.State=reordercats(tbl_probe.State,{'1','2','3'});
tbl_probe.Task=reordercats(tbl_probe.Task,{'1','2'});
F_probe=tbl_probe(tbl_probe.Task=="1",:);
D_probe=tbl_probe(tbl_probe.Task=="2",:);

fprintf('%2.0f/63\n',0)
clear pval_P3_F beta_P3_F pval_P3_D beta_P3_D pval_P3 beta_P3
for nE=1:63
    fprintf('\b\b\b\b\b\b%2.0f/63\n',nE)
    temp_mdl= fitlme(F_probe,sprintf('P3_Ch%g~State*TrCond + (1| SubID)',nE));
    pval_P3_F(nE,:)=double(temp_mdl.Coefficients(:,6));
    beta_P3_F(nE,:)=double(temp_mdl.Coefficients(:,2));
    mod_varnanmesF=temp_mdl.CoefficientNames;

    temp_mdl= fitlme(D_probe,sprintf('P3_Ch%g~State*TrCond + (1| SubID)',nE));
    pval_P3_D(nE,:)=double(temp_mdl.Coefficients(:,6));
    beta_P3_D(nE,:)=double(temp_mdl.Coefficients(:,2));
    mod_varnanmesD=temp_mdl.CoefficientNames;
    
    temp_mdl= fitlme(tbl_probe,sprintf('P3_Ch%g~Task*State*TrCond + (1| SubID)',nE));
    pval_P3(nE,:)=double(temp_mdl.Coefficients(:,6));
    beta_P3(nE,:)=double(temp_mdl.Coefficients(:,2));
    mod_varnanmes=temp_mdl.CoefficientNames;
%     temp_mdl= fitlme(tbl_probe,sprintf('RN_Ch%g~Task*State*CorrCond + (1| SubID)',nE));
%     pval_ERN(nE,:)=double(temp_mdl.Coefficients(:,6));
%     beta_ERN(nE,:)=double(temp_mdl.Coefficients(:,2));
%         mod_varnanmes2=temp_mdl.CoefficientNames;
end

%%
load('../BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
figure;
names_cont={'TrCond_1','State_2:TrCond_1','State_3:TrCond_1'};
for nsubplot=1:3
    subplot(2,3,nsubplot)
    nrow=match_str(mod_varnanmesF,names_cont{1,nsubplot});
    temp_topo=beta_P3_F(:,nrow)';
    temp_topo(pval_P3_F(:,nrow)>fdr(pval_P3_F(:,2:end),0.05))=0;
    simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
        caxis([-1 1]*3)
    title(mod_varnanmesF{nrow})
    colorbar;
    
     subplot(2,3,nsubplot+3)
    nrow=match_str(mod_varnanmesD,names_cont{1,nsubplot});
    temp_topo=beta_P3_D(:,nrow)';
    temp_topo(pval_P3_D(:,nrow)>fdr(pval_P3_D(:,2:end),0.05))=0;
    simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
        caxis([-1 1]*3)
    title(mod_varnanmesD{nrow})
    colorbar;
end

%%
figure;
names_cont={'TrCond_1','State_2:TrCond_1','State_3:TrCond_1','Task_2:State_2:TrCond_1','Task_2:State_3:TrCond_1'};
for nsubplot=1:5
    subplot(1,5,nsubplot)
    nrow=match_str(mod_varnanmes,names_cont{1,nsubplot});
    temp_topo=beta_P3(:,nrow)';
    temp_topo(pval_P3(:,nrow)>fdr(pval_P3(:,2:end),0.05))=0;
    simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
        caxis([-1 1]*3)
    title(mod_varnanmes{nrow})
    colorbar;
end
% figure;
% names_cont={'CorrCond_1','State_1:CorrCond_1','State_3:CorrCond_1';'Task_1:CorrCond_1','Task_1:State_1:CorrCond_1','Task_1:State_3:CorrCond_1'};
% for nsubplot=1:3
%     subplot(2,3,nsubplot)
%     nrow=match_str(mod_varnanmes2,names_cont{1,nsubplot});
%     temp_topo=beta_ERN(:,nrow)';
%     temp_topo(pval_ERN(:,nrow)>0.1)=0;
%     simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
%         caxis([-1 1]*3)
%     title(mod_varnanmes2{nrow})
%     colorbar;
%     
%      subplot(2,3,nsubplot+3)
%     nrow=match_str(mod_varnanmes2,names_cont{2,nsubplot});
%     temp_topo=beta_ERN(:,nrow)';
%     temp_topo(pval_ERN(:,nrow)>0.1)=0;
%     simpleTopoPlot2(temp_topo, pos', labels,0,'parula',0,lay,[]);
%         caxis([-1 1]*3)
%     title(mod_varnanmes2{nrow})
%     colorbar;
% end
