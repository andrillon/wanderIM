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

    %     temp_data=temp_data-repmat(mean(temp_data,1),[size(temp_data,1) 1 1]); % re-reference the data to the average
    ontrial_erp(n,:,:)=(nanmean(temp_data,3));
    
    standardTr=find(test_res(:,5)~=test_res(:,6));
    STD_erp(n,:,:)=(nanmean(temp_data(:,:,standardTr),3));
    standardTr_F=find(test_res(:,2)==1 & test_res(:,5)~=test_res(:,6));
    STD_erp_F(n,:,:)=(nanmean(temp_data(:,:,standardTr_F),3));
    standardTr_D=find(test_res(:,2)==2 & test_res(:,5)~=test_res(:,6));
    STD_erp_D(n,:,:)=(nanmean(temp_data(:,:,standardTr_D),3));
    
    deviantTr=find(test_res(:,5)==test_res(:,6));
    DEV_erp(n,:,:)=(nanmean(temp_data(:,:,deviantTr),3));
    deviantTr_F=find(test_res(:,2)==1 & test_res(:,5)==test_res(:,6));
    DEV_erp_F(n,:,:)=(nanmean(temp_data(:,:,deviantTr_F),3));
    deviantTr_D=find(test_res(:,2)==2 & test_res(:,5)==test_res(:,6));
    DEV_erp_D(n,:,:)=(nanmean(temp_data(:,:,deviantTr_D),3));
    fprintf('... ... %g standard and %g deviant trials\n',length(standardTr),length(deviantTr))
    
%     temp_data=D(match_str(D.chantype,'EEG'),these_times2,:)-mean(D(match_str(D.chanlabels,{'TP7','TP8'}),these_times2,:),1);
    
    these_trials=find(test_res(:,5)==test_res(:,6) & test_res(:,11)==1);
    corrDEV_erp(n,:,:)=(nanmean(temp_data(:,:,these_trials),3));
    
    these_trials=find(test_res(:,5)==test_res(:,6) & test_res(:,11)==0);
    uncorrDEV_erp(n,:,:)=(nanmean(temp_data(:,:,these_trials),3));
    
    these_trials=find(test_res(:,5)~=test_res(:,6) & test_res(:,12)==1);
    corrSTD_erp(n,:,:)=(nanmean(temp_data(:,:,these_trials),3));
    
    these_trials=find(test_res(:,5)~=test_res(:,6) & test_res(:,12)==0);
    uncorrSTD_erp(n,:,:)=(nanmean(temp_data(:,:,these_trials),3));
    
    allRTs=test_res(:,10)-test_res(:,8);
    corrTr=find((~isnan(test_res(:,11)) & test_res(:,11)==1) | (~isnan(test_res(:,12)) & test_res(:,12)==1));
    corrRTs=allRTs(corrTr);
    corrTr(isnan(corrRTs) | corrRTs>0.9 | corrRTs<0.2)=[];
    corrRTs(isnan(corrRTs) | corrRTs>0.9 | corrRTs<0.2)=[];
    timeW=-0.4*D.fsample:0.3*D.fsample-1;
    tempERP=nan(63,length(timeW),length(corrRTs));
    for nRT=1:length(corrRTs)
        tempERP(:,:,nRT)=temp_data(:,timeW+round(corrRTs(nRT)*D.fsample)+D.indsample(0),corrTr(nRT));
    end
    tempERP=tempERP-repmat(nanmean(tempERP(:,1:0.2*D.fsample,:),2),[1 size(tempERP,2) 1]);
    Corr_erp(n,:,:)=(nanmean(tempERP,3));
    
    uncorrTr=find((~isnan(test_res(:,11)) & test_res(:,11)==0) | (~isnan(test_res(:,12)) & test_res(:,12)==0));
    uncorrRTs=allRTs(uncorrTr);
    corrRTs=allRTs(uncorrTr);
    corrTr(isnan(corrRTs) | corrRTs>0.9 | corrRTs<0.2)=[];
    corrRTs(isnan(corrRTs) | corrRTs>0.9 | corrRTs<0.2)=[];
    timeW=-0.4*D.fsample:0.3*D.fsample-1;
    tempERP=nan(63,length(timeW),length(corrRTs));
    for nRT=1:length(corrRTs)
        tempERP(:,:,nRT)=temp_data(:,timeW+round(corrRTs(nRT)*D.fsample)+D.indsample(0),corrTr(nRT));
    end
    tempERP=tempERP-repmat(nanmean(tempERP(:,1:0.2*D.fsample,:),2),[1 size(tempERP,2) 1]);
    Uncorr_erp(n,:,:)=(nanmean(tempERP,3));
    
    for ntask=1:2
        allRTs=test_res(:,10)-test_res(:,8);
        corrTr=find(((~isnan(test_res(:,11)) & test_res(:,11)==1) | (~isnan(test_res(:,12)) & test_res(:,12)==1)) & test_res(:,2)==ntask);
        corrRTs=allRTs(corrTr);
        corrTr(isnan(corrRTs) | corrRTs>0.9 | corrRTs<0.2)=[];
        corrRTs(isnan(corrRTs) | corrRTs>0.9 | corrRTs<0.2)=[];
        timeW=-0.4*D.fsample:0.3*D.fsample-1;
        tempERP=nan(63,length(timeW),length(corrRTs));
        for nRT=1:length(corrRTs)
            tempERP(:,:,nRT)=temp_data(:,timeW+round(corrRTs(nRT)*D.fsample)+D.indsample(0),corrTr(nRT));
        end
        tempERP=tempERP-repmat(nanmean(tempERP(:,1:0.2*D.fsample,:),2),[1 size(tempERP,2) 1]);
        if ntask==1
            Corr_erp_F(n,:,:)=(nanmean(tempERP,3));
        elseif ntask==2
            Corr_erp_D(n,:,:)=(nanmean(tempERP,3));
        end
            
        uncorrTr=find(((~isnan(test_res(:,11)) & test_res(:,11)==0) | (~isnan(test_res(:,12)) & test_res(:,12)==0)) & test_res(:,2)==ntask);
        uncorrRTs=allRTs(uncorrTr);
        corrRTs=allRTs(uncorrTr);
        corrTr(isnan(corrRTs) | corrRTs>0.9 | corrRTs<0.2)=[];
        corrRTs(isnan(corrRTs) | corrRTs>0.9 | corrRTs<0.2)=[];
        timeW=-0.4*D.fsample:0.3*D.fsample-1;
        tempERP=nan(63,length(timeW),length(corrRTs));
        for nRT=1:length(corrRTs)
            tempERP(:,:,nRT)=temp_data(:,timeW+round(corrRTs(nRT)*D.fsample)+D.indsample(0),corrTr(nRT));
        end
        tempERP=tempERP-repmat(nanmean(tempERP(:,1:0.2*D.fsample,:),2),[1 size(tempERP,2) 1]);
        if ntask==1
            Uncorr_erp_F(n,:,:)=(nanmean(tempERP,3));
        elseif ntask==2
            Uncorr_erp_D(n,:,:)=(nanmean(tempERP,3));
        end        
    end
end

%%
xTime=D.time(1):1/D.fsample:D.time(end);

this_ch=match_str(D.chanlabels,'Cz');
figure; set(gcf,'position',[440   247   958   551])
subplot(2,2,1); format_fig;
my_ERP=squeeze(STD_erp(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,[1 1 1]*0.5,0,'-',0.5,1,0,1,1);
my_ERP=squeeze(DEV_erp(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,[1 1 1]*0,0,'-',0.5,1,0,1,1);
xlabel('Time from stim (s)')
ylabel('Amplitude')
xlim([-0.1 1])

subplot(2,2,2); format_fig;
xTime2=-0.4:1/D.fsample:0.3-1/D.fsample;
my_ERP=squeeze(Corr_erp(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime2,my_ERP,0,'b',0,'-',0.5,1,0,1,1);
my_ERP=squeeze(Uncorr_erp(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime2,my_ERP,0,'r',0,'-',0.5,1,0,1,1);
xlabel('Time from resp. (s)')
ylabel('Amplitude')
xlim([-0.4 0.4])

subplot(2,2,3); format_fig;
my_ERP=squeeze(DEV_erp(:,this_ch,:))-squeeze(STD_erp(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,[1 1 1]*0.5,0,'-',0.5,1,0,1,1);
xlabel('Time from stim (s)')
ylabel('Amplitude')
xlim([-0.1 1])
line(xlim,[0 0],'LineStyle','--','Color',[1 1 1]*0.7)

subplot(2,2,4); format_fig;
my_ERP=squeeze(Uncorr_erp(:,this_ch,:))-squeeze(Corr_erp(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime2,my_ERP,0,'k',0,'-',0.5,1,0,1,1);
xlabel('Time from resp. (s)')
ylabel('Amplitude')
xlim([-0.4 0.4])
line(xlim,[0 0],'LineStyle','--','Color',[1 1 1]*0.7)
%%
Colors_Cond={[1 0.5 0.1],[0.4 0.65 0.8]};

figure; set(gcf,'Position',[ 2     1   462   804])
this_ch=match_str(D.chanlabels,'Pz');
subplot(3,1,1); format_fig; hp=[];
my_ERP=squeeze(STD_erp_F(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
[~,hp(1)]=simpleTplot(xTime,my_ERP,0,Colors_Cond{2},0,'-',0.5,1,10,1,1);
my_ERP=squeeze(DEV_erp_F(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
[~,hp(2)]=simpleTplot(xTime,my_ERP,0,Colors_Cond{1},0,'-',0.5,1,10,1,1);
xlabel('Time (s)')
ylabel('ERP amp.')
xlim([-0.2 1])
legend(hp,{'STD','DEV'})
title('Face Task')

subplot(3,1,2); format_fig;
my_ERP=squeeze(STD_erp_D(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,Colors_Cond{2},0,'-',0.5,1,10,1,1);
my_ERP=squeeze(DEV_erp_D(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,Colors_Cond{1},0,'-',0.5,1,10,1,1);
xlabel('Time (s)')
ylabel('ERP amp.')
xlim([-0.2 1])
title('Digit Task')

subplot(3,1,3); format_fig;  hp2=[];
my_ERP=squeeze(DEV_erp_F(:,this_ch,:))-squeeze(STD_erp_F(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
[~,hp2(1)]=simpleTplot(xTime,my_ERP,0,[0.8 0.4 0],[2 0.05 0.05 1000],'-',0.5,1,10,1,1);
my_ERP=squeeze(DEV_erp_D(:,this_ch,:))-squeeze(STD_erp_D(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
[~,hp2(2)]=simpleTplot(xTime,my_ERP,0,[0.8 0.1 0.2],[2 0.05 0.05 1000],':',0.5,1,10,1,1);
xlabel('Time (s)')
ylabel('DEV-STD')
xlim([-0.2 1])
legend(hp2,{'Face','Digit'})

%%
Colors_Cond={[1 0.5 0.1],[0.4 0.65 0.8]};
figure; set(gcf,'Position',[ 2     1   462   804])
this_ch=match_str(D.chanlabels,'Cz');
subplot(3,1,1); format_fig; hp=[];
my_ERP=squeeze(Corr_erp_F(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
[~,hp(1)]=simpleTplot(xTime2,my_ERP,0,Colors_Cond{2},0,'-',0.5,1,10,1,1);
my_ERP=squeeze(Uncorr_erp_F(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
[~,hp(2)]=simpleTplot(xTime2,my_ERP,0,Colors_Cond{1},0,'-',0.5,1,10,1,1);
xlabel('Time to resp (s)')
ylabel('ERP amp.')
xlim([-0.4 0.4])
legend(hp,{'Corr','Inc'})
title('Face Task')

subplot(3,1,2); format_fig;
my_ERP=squeeze(Corr_erp_D(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime2,my_ERP,0,Colors_Cond{2},0,'-',0.5,1,10,1,1);
my_ERP=squeeze(Uncorr_erp_D(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime2,my_ERP,0,Colors_Cond{1},0,'-',0.5,1,10,1,1);
xlabel('Time to resp (s)')
ylabel('ERP amp.')
xlim([-0.4 0.4])
title('Digit Task')

subplot(3,1,3); format_fig;  hp2=[];
my_ERP=squeeze(Uncorr_erp_F(:,this_ch,:))-squeeze(Corr_erp_F(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
[~,hp2(1)]=simpleTplot(xTime2,my_ERP,0,[0.8 0.4 0],[2 0.05 0.05 1000],'-',0.5,1,10,1,1);
my_ERP=squeeze(Uncorr_erp_D(:,this_ch,:))-squeeze(Corr_erp_D(:,this_ch,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
[~,hp2(2)]=simpleTplot(xTime2,my_ERP,0,[0.8 0.1 0.2],[2 0.05 0.05 1000],':',0.5,1,10,1,1);
xlabel('Time to resp (s)')
ylabel('Inc-Cor')
xlim([-0.4 0.4])
legend(hp2,{'Face','Digit'})

%% Cluster permutation
time_W=[0.4 0.6];
% retrieve electrode locations and layout
% % load('EasyCap64_PsychBuilding.mat')
% % elec = [];
% % elec.elecpos = [chaninfo.X ; chaninfo.Y ; chaninfo.Z]';
% % elec.label = {chaninfo(1:63).labels}';
% % elec.chanpos = [chaninfo.X ; chaninfo.Y ; chaninfo.Z]';
% % cfg=[];
% % cfg.elec=elec;
% % [layout] = ft_prepare_layout(cfg);
% % layout.chaninfo=chaninfo;
% % save(['EasyCap64_layout'],'layout');
load(['EasyCap64_layout'])
labels={layout.chaninfo(1:63).labels};

amp_DEV=[];amp_STD=[];
amp_DEV(:,:,1)=squeeze(nanmean(DEV_erp(:,:,xTime>time_W(1) & xTime<time_W(2)),3));
amp_STD(:,:,1)=squeeze(nanmean(STD_erp(:,:,xTime>time_W(1) & xTime<time_W(2)),3));

mypaths.spm=path_spm8;
mypaths.eeglab=path_eeglab;
[stat] = simple_topoplot_clusperm(amp_DEV,amp_STD,0.3,labels,D.fsample,0.1,0.05,1000,layout,'EasyCap64_layout',mypaths);

if ~isempty(stat.posclusters)
cluster_pos_ID=find([stat.posclusters.prob]<0.05);
elec_cluster_pos=(find(ismember(stat.posclusterslabelmat,cluster_pos_ID)));
else
cluster_pos_ID=[];
elec_cluster_pos=[];
end
if ~isempty(stat.negclusters)
cluster_neg_ID=find([stat.negclusters.prob]<0.05);
elec_cluster_neg=(find(ismember(stat.negclusterslabelmat,cluster_neg_ID)));
else
    cluster_neg_ID=[];
    elec_cluster_neg=[];
end
save('clusters_P3','time_W','elec_cluster_pos','elec_cluster_neg')

%%
time_W=[0.4 0.6];
load(['EasyCap64_layout'])
labels={layout.chaninfo(1:63).labels};

amp_DEV=[];amp_STD=[];
amp_DEV(:,:,1)=squeeze(nanmean(DEV_erp_F(:,:,xTime>time_W(1) & xTime<time_W(2)),3));
amp_STD(:,:,1)=squeeze(nanmean(STD_erp_F(:,:,xTime>time_W(1) & xTime<time_W(2)),3));

[stat] = simple_topoplot_clusperm(amp_DEV,amp_STD,0.3,labels,D.fsample,0.05,0.01,1000,layout,'EasyCap64_layout',mypaths);
caxis([-8 8])

time_W=[0.4 0.6];
amp_DEV=[];amp_STD=[];
amp_DEV(:,:,1)=squeeze(nanmean(DEV_erp_D(:,:,xTime>time_W(1) & xTime<time_W(2)),3));
amp_STD(:,:,1)=squeeze(nanmean(STD_erp_D(:,:,xTime>time_W(1) & xTime<time_W(2)),3));
[stat] = simple_topoplot_clusperm(amp_DEV,amp_STD,0.3,labels,D.fsample,0.05,0.01,1000,layout,'EasyCap64_layout',mypaths);
caxis([-8 8])

%%

time_W=[-0.1 0];
load(['EasyCap64_layout'])
labels={layout.chaninfo(1:63).labels};
mypaths.spm=path_spm8;
mypaths.eeglab=path_eeglab;

amp_CORR=[];amp_INC=[];
amp_CORR(:,:,1)=squeeze(nanmean(Corr_erp_F(:,:,xTime2>time_W(1) & xTime2<time_W(2)),3));
amp_INC(:,:,1)=squeeze(nanmean(Uncorr_erp_F(:,:,xTime2>time_W(1) & xTime2<time_W(2)),3));

[stat] = simple_topoplot_clusperm(amp_INC,amp_CORR,-0.1,labels,D.fsample,0.05,0.05,1000,layout,'EasyCap64_layout',mypaths);
caxis([-5 5])

amp_CORR=[];amp_INC=[];
amp_CORR(:,:,1)=squeeze(nanmean(Corr_erp_D(:,:,xTime2>time_W(1) & xTime2<time_W(2)),3));
amp_INC(:,:,1)=squeeze(nanmean(Uncorr_erp_D(:,:,xTime2>time_W(1) & xTime2<time_W(2)),3));
[stat] = simple_topoplot_clusperm(amp_INC,amp_CORR,-0.1,labels,D.fsample,0.05,0.05,1000,layout,'EasyCap64_layout',mypaths);
caxis([-5 5])
