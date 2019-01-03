%% Generate Topographies for Alpha-Theta Ratio

figure;
format_fig;
addpath(genpath(path_eeglab));
temp_topo=mean(all_corrs,1); % vector of 63 values
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
% caxis([0 2.5])
rmpath(genpath(path_eeglab));
colorbar;
caxis([min(temp_topo) max(temp_topo)])
%% maxnegpkamp
clear all;
% close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
% addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
% addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);

% loop across trials for baseline blocks
all_amp_Waves=[];
histc_record = zeros(20,50); % Documenting purposes; rows = participants, columns = bin/counts
collect_all=[];

all_corrs=[];
for n=1:length(bsl_files)
    
    if n==5
        continue
    end
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    
    
    all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)
        
    for nE = 1:63
    nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
    thr_wave = prctile(nE_Waves(:,9),80);
    all_corrs = [all_corrs mean(nE_Waves(nE_Waves(:,9)>=thr_wave,9))];
       
    
    
    end
end

%% p2p amplitude
clear all;
% close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
% addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
% addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);

% loop across trials for baseline blocks
all_amp_Waves=[];
histc_record = zeros(20,50); % Documenting purposes; rows = participants, columns = bin/counts
collect_all=[];

all_corrs=[];
for n=1:length(bsl_files)
    
    if n==5
        continue
    end
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    
    
    all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)
        
    for nE = 1:63
    nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
    thr_wave = prctile(nE_Waves(:,9),80);
    all_corrs = [all_corrs mean(nE_Waves(nE_Waves(:,9)>=thr_wave,4))];
       
    
    
    end
end

all_corrs = vec2mat(all_corrs,63);

%% nTheta
clear all;
% close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
% addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
% addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_nfEEG_S3*.mat']);

% loop across trials for baseline blocks
all_amp_Waves=[];
histc_record = zeros(20,50); % Documenting purposes; rows = participants, columns = bin/counts
collect_all=[];

all_corrs=[];
for n=1:length(bsl_files)
    
    if n==5
        continue
    end
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    
    
    all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)
        
    for nE = 1:63
    nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
    thr_wave = prctile(nE_Waves(:,9),80);
    all_corrs = [all_corrs length(nE_Waves(nE_Waves(:,9)>=thr_wave,4))];
       
    
    
    end
end

all_corrs = vec2mat(all_corrs,63);
        