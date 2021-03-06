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

fixedSW_Thr=75;

%% Align electrodes
CanonicalTemplate=load(['/Users/tand0009/WorkGit/projects/inprogress/wanderIM' filesep 'eeg' filesep 'EasyCap64_layout']);

MBIelecs=load('/Users/tand0009/Work/PostDoc/Monash/Wanderlust/Analyses/BrainVision/myLayout_BV64.mat');
MBIlabels=MBIelecs.lay_BV64.labels;
MICCNelecs=load('/Users/tand0009/WorkGit/projects/inprogress/wanderIM/BrainVision_63ChLayout.mat');
MICCNelecslabels=MICCNelecs.labels;
for nE=1:length(MICCNelecslabels)
    if isempty(match_str(CanonicalTemplate.layout.label, MICCNelecslabels{nE}))
        MICCNidx(nE)=nan;
    else
        MICCNidx(nE)=match_str(CanonicalTemplate.layout.label, MICCNelecslabels{nE});
    end
end
for nE=1:length(MBIlabels)
    if isempty(match_str(CanonicalTemplate.layout.label, MBIlabels{nE}))
        MBIidx(nE)=nan;
    else
        MBIidx(nE)=match_str(CanonicalTemplate.layout.label, MBIlabels{nE});
    end
end

%% loop across trials for baseline blocks
eeg_path=[root_path_adhd filesep 'preproc_eeg'];
behav_path=[root_path_adhd filesep 'behav'];
localsleep_files=dir([eeg_path filesep 'MWADHD_cont_twa2_*.mat']);
for n=1:5%length(localsleep_files)
    % load file with spm
    filename=localsleep_files(n).name;
    
    
    
    % load behavioural results
    SubID=filename;
    SubID=SubID(findstr(SubID,'_4')+1:findstr(SubID,'_4')+3);
    D=spm_eeg_load([eeg_path filesep 'nfEEG_S' SubID '.mat']);
    data=D(1:63,:,:); % D contains the data with channels * time * trials
    data=data-repmat(mean(data([17 22],:,:),1),[size(data,1) 1 1]); % D contains the data with channels * time * trials
    fprintf('... loading subject %s\n',SubID)
    load([eeg_path filesep 'MWADHD_cont_twa2_' SubID])
    
    % load behaviour and get vig rating
    behav_file=dir([behav_path filesep 'MWADHD_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    mean_vig(n)=nanmean(probe_res(:,end));
    mean_mb(n)=nanmean(probe_res(:,32)>=3);
    mean_rt(n)=nanmean(test_res(test_res(:,2)==2 & ~isnan(test_res(:,12)),10)-test_res(test_res(:,2)==2 & ~isnan(test_res(:,12)),8));
    mean_stdrt(n)=nanstd(test_res(test_res(:,2)==2 & ~isnan(test_res(:,12)),10)-test_res(test_res(:,2)==2 & ~isnan(test_res(:,12)),8));
    mean_go(n)=nanmean(test_res(test_res(:,2)==2 & ~isnan(test_res(:,12)),12)==0);
    mean_nogo(n)=nanmean(test_res(test_res(:,2)==2 & ~isnan(test_res(:,11)),11)==0);
    
    all_Waves(abs(all_Waves(:,11))>200 | abs(all_Waves(:,9))>200,:)=[];
    for nE=1:64
        nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
        %         thr_Wave(n,nE) = prctile(all_Waves(:,9),10); %Threshold Wave (80th percentile) by P2P
        %         n_Wave(n,nE) = sum((nE_Waves(:,9))<=thr_Wave(n,nE))/(D.nsamples/D.fsample);
        thr_Wave(n,nE) = prctile(all_Waves(:,11)-all_Waves(:,9),90); %Threshold Wave (80th percentile) by P2P
        n_Wave(n,nE) = sum((nE_Waves(:,11)-nE_Waves(:,9))>=fixedSW_Thr)/(D.nsamples/D.fsample);
    end
    
    % MBI: Cz is channel 14
    thisch=25; %Pz
    nE_Waves = all_Waves(all_Waves(:,3)==thisch,:);
    wvp2p=double(nE_Waves(:,11)-nE_Waves(:,9));
    negmap=nE_Waves(:,9);
    frewv=1./((nE_Waves(:,7)-nE_Waves(:,5))/D.fsample);
    %     selectedwves=(negmap<=thr_Wave(n,thisch)) & (frewv<4);
    %     selectedwves=(wvp2p>=thr_Wave(n,thisch)) & (frewv<4);
    selectedwves=(wvp2p>=fixedSW_Thr) & (frewv<4);
    wvidx=double(nE_Waves(selectedwves,5));
    Cz_wavedens(n)=length(wvidx)/(D.nsamples/D.fsample);
    Cz_waveamp(n)=mean(wvp2p(selectedwves));
    if length(wvidx)>20
        tempERP=nan(length(wvidx),length((-0.5*D.fsample:1*D.fsample)));
        maxabsamp=nan(length(wvidx),1);
        for nw=1:length(wvidx)
            if min(wvidx(nw)+(-0.5*D.fsample:1*D.fsample))>0 && max(wvidx(nw)+(-0.5*D.fsample:1*D.fsample))<=size(D,2)
                maxabsamp(nw)=max(abs(data(thisch,wvidx(nw)+(-0.5*D.fsample:1*D.fsample),1)));
                if maxabsamp(nw)<150
                    tempERP(nw,:)=data(thisch,wvidx(nw)+(-0.5*D.fsample:1*D.fsample),1);
                end
            end
        end
        tempERP=tempERP-repmat(nanmean(tempERP(:,1:0.2*D.fsample),2),[1 size(tempERP,2)]);
        mean_ERPwave(n,:)=nanmean(tempERP);
    else
        mean_ERPwave(n,:)=nan(1,length((-0.5*D.fsample:1*D.fsample)));
    end
    %                 num_Waves(nE,npr)=length(cell2mat(twa_results.channels(nE).maxnegpkamp));
    %             amp_Waves(nE,npr)=mean(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))+abs(cell2mat(twa_results.channels(nE).maxpospkamp)));
end

%%
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
localsleep_files=dir([eeg_path filesep 'MWCTR_cont_twa2_*.mat']);
nc=0;
for n=1:length(localsleep_files)
    % load file with spm
    filename=localsleep_files(n).name;
    
    
    % load behavioural results
    SubID=filename;
    SubID=SubID(findstr(SubID,'_3')+1:findstr(SubID,'_3')+3);
    D=spm_eeg_load([eeg_path filesep 'nfEEG_S' SubID '.mat']);
    data=D(1:63,:,:); % D contains the data with channels * time * trials
    data=data-repmat(mean(data([10 21],:,:),1),[size(data,1) 1 1]); % D contains the data with channels * time * trials
    fprintf('... loading subject %s\n',SubID)
    load([eeg_path filesep 'MWCTR_cont_twa2_' SubID])
    
    % load behaviour and get vig rating
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    
    all_Waves(abs(all_Waves(:,11))>200 | abs(all_Waves(:,9))>200,:)=[];
    if prctile(all_Waves(:,11)-all_Waves(:,9),90)>100
        continue;
    end
    nc=nc+1;
    for nE=1:63
        nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
        %         thr_Wave_ctr(n,nE) = prctile(all_Waves(:,9),10); %Threshold Wave (80th percentile) by P2P
        %         n_Wave_ctr(n,nE) = sum((nE_Waves(:,9))<=thr_Wave_ctr(n,nE))/(D.nsamples/D.fsample);
        thr_Wave_ctr(nc,nE) = prctile(all_Waves(:,11)-all_Waves(:,9),90); %Threshold Wave (80th percentile) by P2P
        n_Wave_ctr(nc,nE) = sum((nE_Waves(:,11)-nE_Waves(:,9))>=fixedSW_Thr)/(D.nsamples/D.fsample);
    end
    mean_vig_ctr(nc)=nanmean(probe_res(:,end));
    mean_mb_ctr(nc)=nanmean(probe_res(:,32)>=3);
    mean_rt_ctr(nc)=nanmean(test_res(test_res(:,2)==2 & ~isnan(test_res(:,12)),10)-test_res(test_res(:,2)==2 & ~isnan(test_res(:,12)),8));
    mean_stdrt_ctr(nc)=nanstd(test_res(test_res(:,2)==2 & ~isnan(test_res(:,12)),10)-test_res(test_res(:,2)==2 & ~isnan(test_res(:,12)),8));
    mean_go_ctr(nc)=nanmean(test_res(test_res(:,2)==2 & ~isnan(test_res(:,12)),12)==0);
    mean_nogo_ctr(nc)=nanmean(test_res(test_res(:,2)==2 & ~isnan(test_res(:,11)),11)==0);
    
    % MICCN: Cz is channel 24
    thisch=13; %Pz
    nE_Waves = all_Waves(all_Waves(:,3)==thisch,:);
    wvp2p=double(nE_Waves(:,11)-nE_Waves(:,9));
    negmap=nE_Waves(:,9);
    frewv=1./((nE_Waves(:,7)-nE_Waves(:,5))/D.fsample);
    %     selectedwves=(negmap<=thr_Wave_ctr(n,thisch)) & (frewv<4);
    %     selectedwves=(wvp2p>=thr_Wave_ctr(nc,thisch)) & (frewv<4);
    selectedwves=(wvp2p>=fixedSW_Thr) & (frewv<4);
    wvidx=double(nE_Waves(selectedwves,5));
    Cz_wavedens_ctr(nc)=length(wvidx)/(D.nsamples/D.fsample);
    Cz_waveamp_ctr(nc)=mean(wvp2p(selectedwves));
    if length(wvidx)>20
        tempERP=nan(length(wvidx),length((-0.5*D.fsample:1*D.fsample)));
        maxabsamp=nan(length(wvidx),1);
        for nw=1:length(wvidx)
            if min(wvidx(nw)+(-0.5*D.fsample:1*D.fsample))>0 && max(wvidx(nw)+(-0.5*D.fsample:1*D.fsample))<=size(D,2)
                maxabsamp(nw)=max(abs(data(thisch,wvidx(nw)+(-0.5*D.fsample:1*D.fsample),1)));
                if maxabsamp(nw)<150
                    tempERP(nw,:)=data(thisch,wvidx(nw)+(-0.5*D.fsample:1*D.fsample),1);
                end
            end
        end
        tempERP=tempERP-repmat(nanmean(tempERP(:,1:0.2*D.fsample),2),[1 size(tempERP,2)]);
        mean_ERPwave_ctr(nc,:)=nanmean(tempERP);
    else
        mean_ERPwave_ctr(nc,:)=nan(1,length((-0.5*D.fsample:1*D.fsample)));
    end
    %                 num_Waves(nE,npr)=length(cell2mat(twa_results.channels(nE).maxnegpkamp));
    %             amp_Waves(nE,npr)=mean(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))+abs(cell2mat(twa_results.channels(nE).maxpospkamp)));
end

%%
figure;
scatter(MBIidx(~isnan(MBIidx)),mean(n_Wave(:,MBIidx(~isnan(MBIidx)))),'or')
hold on
scatter(MICCNidx,mean(n_Wave_ctr(:,MICCNidx)),'ob')

%%
figure; format_fig
% plot((mean_ERPwave)','k')
hold on
% plot((mean_ERPwave_ctr)','r')
simpleTplot((-0.5*D.fsample:1*D.fsample)/D.fsample,(mean_ERPwave),0,'r',0,'-',0.5,1,10,0,2);
simpleTplot((-0.5*D.fsample:1*D.fsample)/D.fsample,(mean_ERPwave_ctr),0,'b',0,'-',0.5,1,10,0,2);


%%
figure; set(gcf,'Position',[775   293   465   685])
nsubc=0;
mysubs=[2 3];
for nsubplot=mysubs
    nsubc=nsubc+1;
    subplot(length(mysubs),1,nsubc)
    hold on; format_fig;
    switch nsubplot
        case 1
            temp_ctr=mean_vig_ctr;
            temp_adhd=mean_vig;
            templabel='Vigilance';
        
                    temp_wavedens_adhd=n_Wave(:,match_str(MBIlabels,'Pz'));
            temp_wavedens_ctr=n_Wave_ctr(:,match_str(MICCNelecslabels,'Pz'));
        case 2
            temp_ctr=mean_go_ctr;
            temp_adhd=mean_go;
            templabel='Go err';
        
                    temp_wavedens_adhd=n_Wave(:,match_str(MBIlabels,'Pz'));
            temp_wavedens_ctr=n_Wave_ctr(:,match_str(MICCNelecslabels,'Pz'));
        case 3
            temp_ctr=mean_nogo_ctr;
            temp_adhd=mean_nogo;
            templabel='NoGo err';
        
            temp_wavedens_adhd=nanmean(n_Wave(:,match_str(MBIlabels,{'Fp1','Fp2'})),2);
            temp_wavedens_ctr=nanmean(n_Wave_ctr(:,match_str(MICCNelecslabels,{'Fp1','Fp2'})),2);
        case 4
            temp_ctr=mean_rt_ctr;
            temp_adhd=mean_rt;
            templabel='RT';
            
            temp_wavedens_adhd=n_Wave(:,match_str(MBIlabels,'Pz'));
            temp_wavedens_ctr=n_Wave_ctr(:,match_str(MICCNelecslabels,'Pz'));
            case 5
            temp_ctr=mean_mb_ctr;
            temp_adhd=mean_mb;
            templabel='MB';
            
            temp_wavedens_adhd=n_Wave(:,match_str(MBIlabels,'Pz'));
            temp_wavedens_ctr=n_Wave_ctr(:,match_str(MICCNelecslabels,'Pz'));
    end
%     scatter(temp_wavedens_ctr,(temp_ctr),'Marker','o','MarkerFaceColor','w','MarkerEdgeColor','b','SizeData',72,'LineWidth',3);
%     scatter(temp_wavedens_adhd,temp_adhd,'Marker','o','MarkerFaceColor','w','MarkerEdgeColor','r','SizeData',72,'LineWidth',3);
     scatter(nanmean(temp_wavedens_ctr),nanmean(temp_ctr),'Marker','s','MarkerFaceColor','b','MarkerEdgeColor','k','SizeData',288,'LineWidth',2);
     scatter(nanmean(temp_wavedens_adhd),nanmean(temp_adhd),'Marker','s','MarkerFaceColor','r','MarkerEdgeColor','k','SizeData',288,'LineWidth',2);
     line([1 1]*nanmean(temp_wavedens_adhd),[1 1]*nanmean(temp_adhd)+[-1 1]*sem(temp_adhd),'Color','k','LineWidth',2);
     line([1 1]*nanmean(temp_wavedens_adhd)+[-1 1]*sem(temp_wavedens_adhd),[1 1]*nanmean(temp_adhd),'Color','k','LineWidth',2);

line([1 1]*nanmean(temp_wavedens_ctr),[1 1]*nanmean(temp_ctr)+[-1 1]*sem(temp_ctr),'Color','k','LineWidth',2);
     line([1 1]*nanmean(temp_wavedens_ctr)+[-1 1]*sem(temp_wavedens_ctr),[1 1]*nanmean(temp_ctr),'Color','k','LineWidth',2);

       scatter(temp_wavedens_ctr,(temp_ctr),'Marker','o','MarkerEdgeColor','b','SizeData',72,'LineWidth',3);
    scatter(temp_wavedens_adhd,temp_adhd,'Marker','o','MarkerEdgeColor','r','SizeData',72,'LineWidth',3);
    
     xlabel('Slow-Wave density (.s^{-1})')
    ylabel(templabel)
    [r pV]=corr([temp_wavedens_ctr ; temp_wavedens_adhd],[temp_ctr temp_adhd]' ,'type','spearman');
    fprintf('Wave Density X %s: r=%g p=%g\n',templabel,r,pV)
end
%%
figure; hold on; format_fig;
scatter(Cz_waveamp_ctr,mean_vig_ctr,'Marker','s','MarkerFaceColor','b','MarkerEdgeColor','b','SizeData',72);
scatter(Cz_waveamp,mean_vig,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','SizeData',72);
xlabel('Local slow wave amplitude (\muV)')
ylabel('Subjective Fatigue')
[r pV]=corr([Cz_waveamp_ctr Cz_waveamp]',[mean_vig_ctr mean_vig]','type','spearman')

%%
figure
for nsubplot=1:5
    subplot(1,5,nsubplot); format_fig;
    switch nsubplot
        case 1
            temp_ctr=mean_vig_ctr;
            temp_adhd=mean_vig;
            templabel='Vigilance';
        case 2
            temp_ctr=mean_go_ctr;
            temp_adhd=mean_go;
            templabel='Go err';
        case 3
            temp_ctr=mean_nogo_ctr;
            temp_adhd=mean_nogo;
            templabel='NoGo err';
        case 4
            temp_ctr=mean_rt_ctr;
            temp_adhd=mean_rt;
            templabel='RT';
            case 5
            temp_ctr=mean_mb_ctr;
            temp_adhd=mean_mb;
            templabel='MB';
    end
    reordern_Wave=nan(size(n_Wave,1),size(n_Wave_ctr,2));
    for nE=1:length(MICCNidx)
        if ~isempty(match_str(MBIlabels,CanonicalTemplate.layout.label{nE}))
            reordern_Wave(:,nE)=n_Wave(:,match_str(MBIlabels,CanonicalTemplate.layout.label{nE}));
        end
    end
    [r, pV]=corr([n_Wave_ctr ; reordern_Wave],[temp_ctr temp_adhd]','type','spearman');
    % plot(1:61,r);
   
    addpath(genpath(path_eeglab));
    temp_topo=r;
    tempplotpV=pV;
    stat_thr=0.05;
    topoplot(temp_topo, CanonicalTemplate.layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{find(tempplotpV<=stat_thr),'.','k',24,2});
    rmpath(genpath(path_eeglab));
    title(templabel)
end
set(gca,'XTick',1:61,'XTickLabel',MBIlabels)
