%% ERP
clear all;
% close all;
run ../localdef_wanderIM
figure

% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
% addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);

for nE =1:63;
    
all_amp_Waves=[];
erp_all=[];

for n=1:length(bsl_files) % for each participant (generate a topography of correlation coefficients)
    
    filename=bsl_files(n).name;
    D=spm_eeg_load([eeg_path filesep filename]);
    fprintf('... processing subject %s\n',D.fname)
    
    
    % load behavioural results
    SubID=D.fname;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    
    load([behav_path filesep behav_file.name]);
    
%     left_freq(n)=SubjectInfo.FlickerR; % CAREFUL this is inverted
%     right_freq(n)=SubjectInfo.FlickerL;
    
    these_times=D.indsample(-20):D.indsample(0)-1;
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    temp_data=temp_data-repmat(mean(temp_data([10 21],:,:),1),[size(temp_data,1) 1 1]); % D contains the data with channels * time * trials
    load([eeg_path filesep 'wanderIM_twa2_' SubID])
    
    fz_waves = all_Waves(all_Waves(:,3)==nE,:);
%     thr_wave=prctile(all_Waves(:,4),80);
    thr_wave=prctile(all_Waves(all_Waves(:,3)==nE,4),80); 

    
 
    negzx_all = fz_waves(fz_waves(:,4)<=thr_wave,5);
    npr_all = fz_waves(fz_waves(:,4)<=thr_wave,2);
% negzx_all = fz_waves(:,5);
% npr_all = fz_waves(:,2);

    npr_all(negzx_all>=9250 | negzx_all<=750)=[];
    negzx_all(negzx_all>=9250)=[];
    negzx_all(negzx_all<=750)=[];
    
    
    
    
    erp=[];
    for npr = 1:length(npr_all)
           temp =  temp_data(2,negzx_all(npr)-750:negzx_all(npr)+750,npr_all(npr));
           erp = [erp; temp];
    end
%     
%     subplot(5,4,n);
%     plot(-1.5:0.002:1.5, mean(erp));
%     xlabel(SubID);
%     erp_all=[erp_all; erp];
    
end

subplot(7,9,nE);
plot(-1.5:0.002:1.5,mean(erp_all))
xlab = num2str(nE)
xlabel(xlab)

end

%
% suplabel('Samples (Fs=500)', 'x'); 
% suplabel('Amplitude', 'y');
% title('ERP for Non-Theta(80) Waves for Fz(Ch2) for each participant');
% 
% figure; simpleTplot(-1.5:0.002:1.5, mean(erp_all))
% xlabel('Time (s)')
% ylabel('Amplitude')
% title('Mean ERP for Non-Theta(80) Waves for Fz(Ch2) Across Participants');