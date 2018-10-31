%%%%% Preprocessing check
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
    these_times2=D.indsample(-0.2):D.indsample(1.3);
    these_timesbs=D.indsample(-0.2):D.indsample(0);
    %     temp_data=D(1:63,these_times2,:)-repmat(mean(D(1:63,these_timesbs,:),2),[1 length(these_times2) 1]);
    temp_data=D(match_str(D.chantype,'EEG'),these_times2,:)-mean(D(match_str(D.chanlabels,{'TP7','TP8'}),these_times2,:),1);
    
    %%% aggregate trials between probes
    allRTs=test_res(:,10)-test_res(:,8);
    all_standardTr=[];
    all_deviantTr=[];
    all_standard_ProbeCat=[];
    all_deviant_ProbeCat=[];
    all_corrTr=[];
    all_uncorrTr=[];
    all_corr_ProbeCat=[];
    all_uncorr_ProbeCat=[];
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
            %
            standardTr=intersect(these_trialsidx,find(test_res(:,5)~=test_res(:,6)));
            standardTr=standardTr(end-17:end);
            all_standardTr=[all_standardTr ; standardTr];
            
            deviantTr=intersect(these_trialsidx,find(test_res(:,5)==test_res(:,6)));
            deviantTr=deviantTr(end-1:end);
            all_deviantTr=[all_deviantTr ; deviantTr];
            
            all_standard_ProbeCat=[all_standard_ProbeCat ; repmat([nbl npr these_probes(npr,5) these_probes(npr,32)],length(standardTr),1)];
            all_deviant_ProbeCat=[all_deviant_ProbeCat ; repmat([nbl npr these_probes(npr,5) these_probes(npr,32)],length(deviantTr),1)];
            
            corrTr=intersect(these_trialsidx,find((~isnan(test_res(:,11)) & test_res(:,11)==1) | (~isnan(test_res(:,12)) & test_res(:,12)==1)));
            corrTr=corrTr(end-17:end);
            all_corrTr=[all_corrTr ; corrTr];
            
            uncorrTr=intersect(these_trialsidx,find((~isnan(test_res(:,11)) & test_res(:,11)==0) | (~isnan(test_res(:,12)) & test_res(:,12)==0)));
            uncorrTr=uncorrTr;
            all_uncorrTr=[all_uncorrTr ; uncorrTr];
            
            all_corr_ProbeCat=[all_corr_ProbeCat ; repmat([nbl npr these_probes(npr,5) these_probes(npr,32)],length(corrTr),1)];
            all_uncorr_ProbeCat=[all_uncorr_ProbeCat ; repmat([nbl npr these_probes(npr,5) these_probes(npr,32)],length(uncorrTr),1)];
            %                         corrTr=find((~isnan(test_res(:,11)) & test_res(:,11)==1) | (~isnan(test_res(:,12)) & test_res(:,12)==1));
            
        end
    end
    for ntask=1:2
        for nstate=1:3
            these_tridx=all_standardTr(all_standard_ProbeCat(:,3)==ntask & all_standard_ProbeCat(:,4)==nstate);
            these_tridx2=all_deviantTr(all_deviant_ProbeCat(:,3)==ntask & all_deviant_ProbeCat(:,4)==nstate);
            STD_erp(n,ntask,nstate,:,:)=(mean(temp_data(:,:,these_tridx),3));
            DEV_erp(n,ntask,nstate,:,:)=(mean(temp_data(:,:,these_tridx2),3));
            numSTD_erp(n,ntask,nstate)=length(these_tridx);
            numDEV_erp(n,ntask,nstate)=length(these_tridx2);
            
            
            these_tridx=all_corrTr(all_corr_ProbeCat(:,3)==ntask & all_corr_ProbeCat(:,4)==nstate);
            these_tridx2=all_corrTr(all_uncorr_ProbeCat(:,3)==ntask & all_uncorr_ProbeCat(:,4)==nstate);
            
            corrTr=these_tridx;
            corrRTs=allRTs(these_tridx);
            corrTr(isnan(corrRTs) | corrRTs>1)=[];
            corrRTs(isnan(corrRTs) | corrRTs>1)=[];
            timeW=D.indsample(-0.2):D.indsample(0.3)-1;
            tempERP=nan(63,length(timeW),length(corrRTs));
            for nRT=1:length(corrRTs)
                tempERP(:,:,nRT)=temp_data(:,timeW+round(corrRTs(nRT)*D.fsample),corrTr(nRT));
            end
            Corr_erp(n,ntask,nstate,:,:)=(mean(tempERP,3));
            
            uncorrTr=these_tridx;
            uncorrRTs=allRTs(these_tridx2);
            uncorrTr(isnan(uncorrRTs) | uncorrRTs>1)=[];
            uncorrRTs(isnan(uncorrRTs) | uncorrRTs>1)=[];
            timeW=D.indsample(-0.2):D.indsample(0.3)-1;
            tempERP=nan(63,length(timeW),length(uncorrRTs));
            for nRT=1:length(uncorrRTs)
                tempERP(:,:,nRT)=temp_data(:,timeW+round(uncorrRTs(nRT)*D.fsample),uncorrTr(nRT));
            end
            Uncorr_erp(n,ntask,nstate,:,:)=(mean(tempERP,3));
        end
    end
    
    
    
end


%%
timeP3_W=[0.4 0.7];
xTime=-0.2:1/D.fsample:1.3;
P3_Amp=mean(DEV_erp(:,:,:,13,xTime>timeP3_W(1) & xTime<timeP3_W(2))-STD_erp(:,:,:,13,xTime>timeP3_W(1) & xTime<timeP3_W(2)),5);

f1=figure; format_fig;
f2=figure; format_fig;

for ntask=1:2
    for nstate=1:3
        figure(f1);
        subplot(3,3,3*(ntask-1)+nstate)
        if ntask==3
            sub_idx=squeeze(sum(numDEV_erp(:,:,nstate),2))>4;
            std=squeeze(nanmean(STD_erp(sub_idx,:,nstate,13,:),2));
            dev=squeeze(nanmean(DEV_erp(sub_idx,:,nstate,13,:),2));
            simpleTplot(xTime,std,0,'b',0,'-',0.5,0,0,0,3);
            simpleTplot(xTime,dev,0,'r',0,'-',0.5,0,0,0,3);
        else
            sub_idx=squeeze(sum(numDEV_erp(:,ntask,nstate),2))>4;
            std=squeeze(nanmean(STD_erp(sub_idx,ntask,nstate,13,:),2));
            dev=squeeze(nanmean(DEV_erp(sub_idx,ntask,nstate,13,:),2));
            simpleTplot(xTime,std,0,'b',0,'-',0.5,0,0,0,3);
            simpleTplot(xTime,dev,0,'r',0,'-',0.5,0,0,0,3);
        end
        ylim([-2 10])
        xlim([-0.1 1])
        
        figure(f2);
        if ntask==3
            subplot(1,3,3)
            sub_idx=squeeze(sum(numDEV_erp(:,:,nstate),2))>4;
            simpleBarPlot(nstate,squeeze(mean(P3_Amp(sub_idx,:,nstate),2)),'k',0.9,'r',{1 0 0.05},3);
        else
            subplot(1,3,ntask)
            sub_idx=squeeze(sum(numDEV_erp(:,ntask,nstate),2))>4;
            simpleBarPlot(nstate,squeeze(mean(P3_Amp(sub_idx,ntask,nstate),2)),'k',0.9,'r',{1 0 0.05},3);
        end
    end
end

%%
xTime=-0.2:1/D.fsample:0.3-1/D.fsample;
timeERN_W=[-0.1 0.1];
ERN_Amp=mean(Uncorr_erp(:,:,:,5,xTime>timeERN_W(1) & xTime<timeERN_W(2))-Corr_erp(:,:,:,5,xTime>timeERN_W(1) & xTime<timeERN_W(2)),5);

f1=figure; format_fig;
f2=figure; format_fig;
for ntask=1:2
    for nstate=1:3
        figure(f1);
        subplot(3,3,3*(ntask-1)+nstate)
        if ntask==3
            sub_idx=squeeze(sum(numDEV_erp(:,:,nstate),2))>4;
            std=squeeze(nanmean(Corr_erp(sub_idx,:,nstate,2,:),2));
            dev=squeeze(nanmean(Uncorr_erp(sub_idx,:,nstate,2,:),2));
            simpleTplot(xTime,std,0,'b',0,'-',0.5,0,0,0,3);
            simpleTplot(xTime,dev,0,'r',0,'-',0.5,0,0,0,3);
        else
            sub_idx=squeeze(sum(numDEV_erp(:,:,nstate),2))>4;
            std=squeeze(nanmean(Corr_erp(sub_idx,ntask,nstate,2,:),2));
            dev=squeeze(nanmean(Uncorr_erp(sub_idx,ntask,nstate,2,:),2));
            simpleTplot(xTime,std,0,'b',0,'-',0.5,0,0,0,3);
            simpleTplot(xTime,dev,0,'r',0,'-',0.5,0,0,0,3);
        end
        ylim([-1 10])
        xlim([-0.2 0.3])
        
        figure(f2);
        if ntask==3
            subplot(1,3,3)
            sub_idx=squeeze(sum(numDEV_erp(:,:,nstate),2))>4;
            simpleBarPlot(nstate,squeeze(mean(ERN_Amp(sub_idx,:,nstate),2)),'k',0.9,'r',{1 0 0.05},3);
        else
            subplot(1,3,ntask)
            sub_idx=squeeze(sum(numDEV_erp(:,ntask,nstate),2))>4;
            simpleBarPlot(nstate,squeeze(mean(ERN_Amp(sub_idx,ntask,nstate),2)),'k',0.9,'r',{1 0 0.05},3);
        end
    end
end

%%
timeP3_W=[0.4 0.7];
xTime=-0.2:1/D.fsample:1.3;
xTime2=-0.2:1/D.fsample:0.3-1/D.fsample;

f1=figure; format_fig;
f2=figure; format_fig;
Colors={'k','r','g'};
for ntask=1:2
    for nstate=1:3
        figure(f1); hold on;
        subplot(1,2,ntask)
        sub_idx=squeeze(sum(numDEV_erp(:,ntask,nstate),2))>4;
        std=squeeze(nanmean(STD_erp(sub_idx,ntask,nstate,13,:),2));
        dev=squeeze(nanmean(DEV_erp(sub_idx,ntask,nstate,13,:),2));
        simpleTplot(xTime,dev-std,0,Colors{nstate},0,'-',0.5,0,0,0,3);
        
        ylim([-2 10])
        xlim([-0.1 1])
        
        
        figure(f2); hold on;
        subplot(1,2,ntask)
        sub_idx=squeeze(sum(numDEV_erp(:,ntask,nstate),2))>4;
        std=squeeze(nanmean(Corr_erp(sub_idx,ntask,nstate,5,:),2));
        dev=squeeze(nanmean(Uncorr_erp(sub_idx,ntask,nstate,5,:),2));
        simpleTplot(xTime2,dev-std,0,Colors{nstate},0,'-',0.5,0,0,0,3);
        ylim([-2 10])
        xlim([-0.2 0.3])
    end
end
