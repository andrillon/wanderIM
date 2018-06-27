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
    temp_data=D(match_str(D.chanlabels,'Oz'),these_times2,:);
    %     temp_data=temp_data-repmat(mean(temp_data,1),[size(temp_data,1) 1 1]); % re-reference the data to the average
    ontrial_erp(n,:,:,:)=(mean(temp_data,3));
    
    standardTr=find(test_res(:,5)~=test_res(:,6));
    onSTD_erp(n,:,:,:)=(mean(temp_data(:,:,standardTr),3));
    
    deviantTr=find(test_res(:,5)==test_res(:,6));
    onDEV_erp(n,:,:,:)=(mean(temp_data(:,:,deviantTr),3));
    fprintf('... ... %g standard and %g deviant trials\n',length(standardTr),length(deviantTr))
    
    temp_data=D(match_str(D.chanlabels,'FCz'),these_times2,:)-mean(D(match_str(D.chanlabels,{'TP7','TP8'}),these_times2,:),1);
    
    these_trials=find(test_res(:,5)==test_res(:,6) & test_res(:,11)==1);
    corrDEV_erp(n,:,:,:)=(mean(temp_data(:,:,these_trials),3));
    
    these_trials=find(test_res(:,5)==test_res(:,6) & test_res(:,11)==0);
    uncorrDEV_erp(n,:,:,:)=(mean(temp_data(:,:,these_trials),3));
    
    these_trials=find(test_res(:,5)~=test_res(:,6) & test_res(:,12)==1);
    corrSTD_erp(n,:,:,:)=(mean(temp_data(:,:,these_trials),3));
    
    these_trials=find(test_res(:,5)~=test_res(:,6) & test_res(:,12)==0);
    uncorrSTD_erp(n,:,:,:)=(mean(temp_data(:,:,these_trials),3));
end

%%
figure; format_fig;
xTime=-0.2:1/D.fsample:1;
% average ERP across conditions and use only Oz
this_ch=1; %match_str(D.chanlabels,'FCz');
my_ERP=ontrial_erp; %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'b',0,'-',0.5,1,0,1,1);


figure; subplot(1,2,1); format_fig;
xTime=-0.2:1/D.fsample:1;
my_ERP=onSTD_erp; %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'b',0,'-',0.5,1,0,1,1);
my_ERP=onDEV_erp; %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'r',0,'-',0.5,1,0,1,1);
xlabel('Time (s)')
ylabel('Amplitude')

subplot(1,2,2); format_fig;
my_ERP=onDEV_erp-onSTD_erp; %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'k',0,'-',0.5,1,0,1,1);
xlabel('Time (s)')
ylabel('Amplitude')

%%
figure; subplot(1,2,1); format_fig;
xTime=-0.2:1/D.fsample:1;
my_ERP=corrDEV_erp; %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'b',0,'-',0.5,1,0,1,1);
my_ERP=uncorrDEV_erp; %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'r',0,'-',0.5,1,0,1,1);
xlabel('Time (s)')
ylabel('Amplitude')

subplot(1,2,2); format_fig;
my_ERP=corrSTD_erp; %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'b',0,'-',0.5,1,0,1,1);
my_ERP=uncorrDEV_erp; %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'r',0,'-',0.5,1,0,1,1);
xlabel('Time (s)')
ylabel('Amplitude')

figure; format_fig;
xTime=-0.2:1/D.fsample:1;
simpleTplot(xTime,corrDEV_erp,0,'r',0,'-',0.5,0,0,0,3);
simpleTplot(xTime,uncorrDEV_erp,0,'r',0,'--',0.5,0,0,0,3);
simpleTplot(xTime,corrSTD_erp,0,'b',0,'-',0.5,0,0,0,3);
simpleTplot(xTime,uncorrSTD_erp,0,'b',0,'--',0.5,0,0,0,3);
xlabel('Time (s)')
ylabel('Amplitude')