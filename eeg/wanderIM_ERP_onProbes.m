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
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);

%% loop across trials for baseline blocks
onprobe_erp_MB=nan(length(bsl_files),63,601);
onprobe_logSNR_MB=nan(length(bsl_files),63,5001);
onprobe_logPow_MB=nan(length(bsl_files),63,5001);

onprobe_erp_MW=nan(length(bsl_files),63,601);
onprobe_logSNR_MW=nan(length(bsl_files),63,5001);
onprobe_logPow_MW=nan(length(bsl_files),63,5001);

onprobe_erp_ON=nan(length(bsl_files),63,601);
onprobe_logSNR_ON=nan(length(bsl_files),63,5001);
onprobe_logPow_ON=nan(length(bsl_files),63,5001);

onprobe_logSNR=nan(length(bsl_files),63,5001,60);
onprobe_logPow=nan(length(bsl_files),63,5001,60);
onprobe_erp=nan(length(bsl_files),63,601);

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
    
    param=[];
    param.method='fft'; % fast fourier transform
    param.mindist=0.5; % we want to be able to separate peaks separated by at least 1 Hz
    these_times=D.indsample(-20):D.indsample(0);
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    temp_data=temp_data-repmat(mean(temp_data,1),[size(temp_data,1) 1 1]); % re-reference the data to the average
    
    [logSNR, faxis, logpow]=get_logSNR(temp_data,D.fsample,param);
    
    onprobe_logSNR(n,:,:,:)=logSNR;
    onprobe_logPow(n,:,:,:)=logpow;
    
    
    % sanity check: computer ERP on block onset
    these_times2=D.indsample(-0.2):D.indsample(1);
    these_timesbs=D.indsample(-0.2):D.indsample(0);
    temp_data=D(1:63,these_times2,:)-repmat(mean(D(1:63,these_timesbs,:),2),[1 length(these_times2) 1]);
    temp_data=temp_data-repmat(mean(temp_data,1),[size(temp_data,1) 1 1]); % re-reference the data to the average
    onprobe_erp(n,:,:)=squeeze(mean(temp_data,3));
    
    these_trials=find(~cellfun(@isempty,regexp(D.conditions,'MW')));
    if length(these_trials)>=6
        onprobe_erp_MW(n,:,:)=squeeze(mean(temp_data(:,:,these_trials),3));
        onprobe_logSNR_MW(n,:,:)=squeeze(mean(logSNR(:,:,these_trials),3));
        onprobe_logPow_MW(n,:,:)=squeeze(mean(logpow(:,:,these_trials),3));
%     else
%         onprobe_erp_MW(n,:,:)=nan(1,63,length(these_times2));
%         onprobe_logSNR_MW(n,:,:)=nan(1,63,length(faxis));
%         onprobe_logPow_MW(n,:,:)=nan(1,63,length(faxis));
    end
    
    these_trials=find(~cellfun(@isempty,regexp(D.conditions,'ON')));
    if length(these_trials)>=6
        onprobe_erp_ON(n,:,:)=squeeze(mean(temp_data(:,:,these_trials),3));
        onprobe_logSNR_ON(n,:,:)=squeeze(mean(logSNR(:,:,these_trials),3));
        onprobe_logPow_ON(n,:,:)=squeeze(mean(logpow(:,:,these_trials),3));
%     else
%         onprobe_erp_ON(n,:,:)=nan(1,63,length(these_times2));
%         onprobe_logSNR_ON(n,:,:)=nan(1,63,length(faxis));
%         onprobe_logSNR_ON(n,:,:)=nan(1,63,length(faxis));
    end
    these_trials=find(~cellfun(@isempty,regexp(D.conditions,'MB')));
    if length(these_trials)>=6
        onprobe_erp_MB(n,:,:)=squeeze(mean(temp_data(:,:,these_trials),3));
        onprobe_logSNR_MB(n,:,:)=squeeze(mean(logSNR(:,:,these_trials),3));
        onprobe_logPow_MB(n,:,:)=squeeze(mean(logpow(:,:,these_trials),3));
%     else
%         onprobe_erp_MB(n,:,:)=nan(1,63,length(these_times2));
%         onprobe_logSNR_MB(n,:,:)=nan(1,63,length(faxis));
%         onprobe_logPow_MB(n,:,:)=nan(1,63,length(faxis));
    end
end

%%
figure; format_fig;
xTime=-0.2:1/D.fsample:1;
% average ERP across conditions and use only Oz
this_ch=match_str(D.chanlabels,'Cz');
my_ERP=squeeze(mean(onprobe_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'b',0,'-',0.5,1,0,1,1);


%%
figure; format_fig;
xTime=-0.2:1/D.fsample:1;
% average ERP across conditions and use only Oz
this_ch=match_str(D.chanlabels,'Cz');
my_ERP=squeeze(mean(onprobe_erp_MW(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'b',0,'-',0.5,1,0,1,1);
my_ERP=squeeze(mean(onprobe_erp_ON(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'k',0,'-',0.5,1,0,1,1);
my_ERP=squeeze(mean(onprobe_erp_MB(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'r',0,'-',0.5,1,0,1,1);

%%
myChanLabel='Cz';
figure;
subplot(2,2,1)
plot(faxis,squeeze(nanmean(onprobe_logSNR_MW(:,match_str(D.chanlabels,myChanLabel),:),1)),'b')
xlim([1 30])
hold on; plot(faxis,squeeze(nanmean(onprobe_logSNR_MB(:,match_str(D.chanlabels,myChanLabel),:),1)),'r')
hold on; plot(faxis,squeeze(nanmean(onprobe_logSNR_ON(:,match_str(D.chanlabels,myChanLabel),:),1)),'k')


subplot(2,2,2)
plot(faxis,squeeze(nanmean(onprobe_logPow_MW(:,match_str(D.chanlabels,myChanLabel),:),1)),'b')
xlim([1 30])
hold on; plot(faxis,squeeze(nanmean(onprobe_logPow_MB(:,match_str(D.chanlabels,myChanLabel),:),1)),'r')
hold on; plot(faxis,squeeze(nanmean(onprobe_logPow_ON(:,match_str(D.chanlabels,myChanLabel),:),1)),'k')

subplot(2,2,3)
plot(faxis,squeeze(nanmean(onprobe_logSNR_ON(:,match_str(D.chanlabels,myChanLabel),:)-onprobe_logSNR_MW(:,match_str(D.chanlabels,myChanLabel),:),1)),'b')
hold on
plot(faxis,squeeze(nanmean(onprobe_logSNR_ON(:,match_str(D.chanlabels,myChanLabel),:)-onprobe_logSNR_MB(:,match_str(D.chanlabels,myChanLabel),:),1)),'g')
xlim([1 30])

subplot(2,2,4)
plot(faxis,squeeze(nanmean(onprobe_logPow_ON(:,match_str(D.chanlabels,myChanLabel),:)-onprobe_logPow_MW(:,match_str(D.chanlabels,myChanLabel),:),1)),'b')
hold on
plot(faxis,squeeze(nanmean(onprobe_logPow_ON(:,match_str(D.chanlabels,myChanLabel),:)-onprobe_logPow_MB(:,match_str(D.chanlabels,myChanLabel),:),1)),'g')
xlim([1 30])