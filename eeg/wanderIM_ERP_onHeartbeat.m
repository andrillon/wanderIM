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
bsl_files=dir([eeg_path filesep 'hb_nfEEG_S3*.mat']);

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
    these_times2=D.indsample(-0.1):D.indsample(0.4);
    these_timesbs=D.indsample(-0.1):D.indsample(0);
    %     temp_data=D(1:63,these_times2,:)-repmat(mean(D(1:63,these_timesbs,:),2),[1 length(these_times2) 1]);
    temp_data=D(match_str(D.chanlabels,{'Cz','ECG'}),these_times2,:);
    goodtr=squeeze(max(abs(temp_data(1,:,:)),[],2))<80 & squeeze(max(abs(temp_data(2,:,:)),[],2))<150;
    %     temp_data=temp_data-repmat(mean(temp_data,1),[size(temp_data,1) 1 1]); % re-reference the data to the average
    onHB_erp(n,:,:)=(mean(temp_data(:,:,goodtr),3));
   
    figure;
    plot(-0.1:1/D.fsample:0.4,squeeze(temp_data(2,:,goodtr))');
end

%%
figure; format_fig;
xTime=-0.1:1/D.fsample:0.4;
% average ERP across conditions and use only Oz
this_ch=1; %match_str(D.chanlabels,'FCz');
my_ERP=squeeze(onHB_erp(:,2,:)); %squeeze(mean(ontrial_erp(:,this_ch,:),2));
simpleTplot(xTime,my_ERP,0,'b',0,'-',0.5,1,0,1,1);

