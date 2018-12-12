
clear all;
% close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
% addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);

all_amp_Waves=[];
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
    a_t_all = zeros(63,60);

        for nE = 1:63 % for each channel (generate a correlation coefficients of alpha/theta ratio and fatigue)
    
            for npr=1:size(temp_data,3) % for each probe (generate an alpha/theta ratio)
      
                datainput=squeeze(temp_data(nE,:,npr));
                SamplingRate=500;
                DecibelsFlag=0;
                plotFlag=0;
                signal = squeeze(temp_data(nE, : , npr));
                [faxis pow] = get_PowerSpec(signal, SamplingRate, DecibelsFlag , plotFlag);
                pow_alpha=log(mean(pow(faxis>8 & faxis<13)));
                pow_theta=log(mean(pow(faxis>4 & faxis<7)));
                alpha_theta= pow_alpha/pow_theta;
                a_t_all(nE,npr) = alpha_theta; 
            end
        end
        
        save([eeg_path filesep 'at_ratio_wanderIM_twa2_' SubID],'a_t_all')
    
    end
    
    
    %% Original Data
   close all
   thr_Wave = thr_Wave1
   fz_waves = all_Waves(all_Waves(:,3)==2,:);
   a = fz_waves(fz_waves(:,4)>=thr_Wave,5);
   b = fz_waves(fz_waves(:,4)>=thr_Wave,6);
   c = fz_waves(fz_waves(:,4)>=thr_Wave,7);
   d = fz_waves(fz_waves(:,4)>=thr_Wave,2);
   
   a=a(d==1);% configure these; d==x where x is probe, n=x where x is wave#
   b=b(d==1);
   n = 1;
   nn=1;
   
figure
for n = 1:5
   
%    figure;
   subplot(5,2,nn)
   plot(temp_data(2,a(n):c(n),1));
   hold;
   e = zeros(length(1:c(n)-a(n)+1),1);
   plot(e,'-o','MarkerIndices',[1 b(n)-a(n)+1])

aa = round((a./500).*128);
bb = round((b./500).*128);
cc = round((c./500).*128);   

% figure;
subplot(5,2,nn+1)
plot(filtered_data(2, aa(n):cc(n),1));
hold;
e = zeros(length(1:cc(n)-aa(n)+1),1);
plot(e,'-o','MarkerIndices',[1 bb(n)-aa(n)+1]);

nn=nn+2;
end
%% alphatheta next step


clear all;
% close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
% addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);

all_corrs = zeros(63,60);

for n=1:length(bsl_files)
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);    
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'at_ratio_wanderIM_twa2_' SubID])
    for nE = 1:63
        a = a_t_all(nE,:)';
        b = probe_res(:,38);
        rho = corr(a,b);
        all_corrs(nE,n) = rho;
    end
    
end


%% Generate Topographies

figure;
format_fig;
addpath(genpath(path_eeglab));
temp_topo=mean(all_corrs.*100,2); % vector of 63 values
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
% caxis([0 2.5])
rmpath(genpath(path_eeglab));
colorbar;
caxis([0 max(temp_topo)])
 
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

all_amp_Waves=[];
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
    
    fz_waves = all_Waves(all_Waves(:,3)==2,:);
    thr_wave=prctile(all_Waves(all_Waves(:,3)==2,4),90);
    negzx_all = fz_waves(fz_waves(:,4)>=thr_wave,5);
    npr_all = fz_waves(fz_waves(:,4)>=thr_wave,2);
    npr_all(negzx_all>=9250 | negzx_all<=750)=[];
    negzx_all(negzx_all>=9250)=[];
    negzx_all(negzx_all<=750)=[];
    
    
    
    
    erp=[];
    for npr = 1:length(npr_all)
           temp =  temp_data(2,negzx_all(npr)-750:negzx_all(npr)+750,npr_all(npr));
           erp = [erp; temp];
    end
    
    subplot(5,4,n);
    plot(mean(erp));
    xlabel(SubID);
    
end

suplabel('Samples (Fs=500)', 'x'); 
suplabel('Amplitude', 'y');
title('ERP for Theta(90) Events for Fz(Ch2) for each participant');
    
    
