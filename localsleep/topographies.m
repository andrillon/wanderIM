%% This files contains the code for generating the topography plots
% 1. Master Topography Code for Plotting
% 2. Plot Maxnegpkamp

%% Plot Topography: 'all_corrs'
figure;
format_fig;
addpath(genpath(path_eeglab));
temp_topo=mean(all_corrs,1); % vector of 63 values, the 'mean' depends on original matrix
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
% caxis([0 2.5])
rmpath(genpath(path_eeglab));
colorbar;
caxis([min(temp_topo) max(temp_topo)])
% caxis([-0.1 0.1])

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
    collect_all=[];    
    for nE = 1:63
    nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
    thr_wave = prctile(nE_Waves(:,9),80);
    all_corrs = [all_corrs mean(nE_Waves(nE_Waves(:,9)>=thr_wave,9))]; %topography for mean across participant
    collect_all=[collect_all mean(nE_Waves(nE_Waves(:,9)>=thr_wave,9))]; %topography for each participant           
    end
    
    subplot(5,4,n);
format_fig;
addpath(genpath(path_eeglab));
temp_topo=collect_all; % vector of 63 values
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
% caxis([0 2.5])
rmpath(genpath(path_eeglab));
colorbar;
caxis([min(temp_topo) max(temp_topo)])
xlabel(SubID);
    
    
end

all_corrs = vec2mat(all_corrs,63);

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
        collect_all=[];
    for nE = 1:63
    nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
    thr_wave = prctile(nE_Waves(:,9),80);
    all_corrs = [all_corrs mean(nE_Waves(nE_Waves(:,9)>=thr_wave,4))];
       
     collect_all=[collect_all mean(nE_Waves(nE_Waves(:,9)>=thr_wave,4))]; %topography for each participant           
    end
    
    subplot(5,4,n);
format_fig;
addpath(genpath(path_eeglab));
temp_topo=collect_all; % vector of 63 values
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
% caxis([0 2.5])
rmpath(genpath(path_eeglab));
colorbar;
caxis([4 7])
xlabel(SubID);
    
end

all_corrs = vec2mat(all_corrs,63);

%% nTheta
clear all;
close all;
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


all_corrs=[];
for n=1:length(bsl_files)
    
    if n==5
        continue
    end
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    
    
%     all_Waves_old=all_Waves;
%     thrP2P=60;
%     all_Waves(all_Waves(:,4)>thrP2P,:)=[];
%         fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P) %
    collect_all=[];
    for nE = 1:63
    nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
    thr_wave = prctile(nE_Waves(:,9),80);
%     all_corrs = [all_corrs length(nE_Waves(nE_Waves(:,9)>=thr_wave,4))];
%        
    
    collect_all=[collect_all length(nE_Waves(nE_Waves(:,9)>=thr_wave,4))]; %topography for each participant           
    
    end

    collect_all = zscore(collect_all); %Transform into zscore to account for inter-participant variability (Hung et al.)
    all_corrs = [all_corrs collect_all];
    subplot(5,4,n);
format_fig;
addpath(genpath(path_eeglab));
temp_topo=collect_all; % vector of 63 values
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
% caxis([0 2.5])
rmpath(genpath(path_eeglab));
colorbar;
caxis([min(temp_topo) max(temp_topo)])
xlabel(SubID);
end

all_corrs = vec2mat(all_corrs,63);

%% corr p2p/maxnegpkamp and sleepiness
clear all;
 close all;
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



% loop across trials for baseline blocks
% Vectors for x and y axis


    
maxnegpkamp = [];
sleepiness = [];

all_corrs = [];
for n=1:length(bsl_files)
    
    if n==5
        continue
    end
    
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    filename=bsl_files(n).name;
%     D=spm_eeg_load([eeg_path filesep filename]);
    
    corr_per=[];

    all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)

    

   
    sleepiness_n = (probe_res(:,38)); %sleepiness rating indexed by probe for the participant in the loop
    sleepiness_n = zscore(sleepiness_n);
    sleepiness = [sleepiness ; sleepiness_n];
    
    for nE = 1:63
           nE_waves = all_Waves(all_Waves(:,3)==nE,:);
           thr_wave=prctile(abs(all_Waves(all_Waves(:,3)==nE,4)),80);  
           maxnegpkcorr=[];
            for npr = 1:60
                maxnegpkcorr = [maxnegpkcorr mean(nE_waves(nE_waves(:,2)==npr,4)    )];
            end
        [r,pV]=corr(maxnegpkcorr', probe_res(:,38),'type','spearman');
        all_corrs = [all_corrs r];
    end

end
all_corrs = vec2mat(all_corrs,63);
    
%% Corr: nTheta and sleepiness (slide:47)


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

all_corrs = zeros(63,20);
all_nTheta=[];
all_sleepiness=[];
for n=1:length(bsl_files)
    if n==5;
        continue
    end
    
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);    
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    
    all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)   
    
    
    
    for nE = 1:63
        nE_wave = all_Waves(all_Waves(:,3)==nE,:);
        thr_wave = prctile(all_Waves(all_Waves(:,3)==nE,4),80);
        nE_wave = nE_wave(nE_wave(:,4)>=thr_wave,:);
        a = zeros(60,1);
        for nPr = 1:60
%             a(nPr) =  mean(nE_wave(nE_wave(:,2)==nPr,4));
              a(nPr) = length(nE_wave(nE_wave(:,2)==nPr,1));  
%             if isnan(a(nPr))==true
%                 a(nPr)=0;
%             end
        end
        
        
        
            b = probe_res(:,38); % sleepiness rating
            all_nTheta = [all_nTheta ; a];
            all_sleepiness =[all_sleepiness ; b];
%         b(a==0)=[];
%         a(a==0)=[];
        [rho] = corr(a,b,'type', 'spearman'); % correlation of theta number and sleepiness rating
        all_corrs(nE,n)=rho;
    end
end

%% Corr: nTheta and sleepiness (slide:47) with z score


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

all_corrs = zeros(63,20);
all_nTheta=[];
all_sleepiness=[];
for n=1:length(bsl_files)
    if n==5;
        continue
    end
    
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);    
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    
    all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)   
    collect_all=[];
    all_nTheta = [];
    
    all_sleepiness = []; 
    
    for nE = 1:63
        nE_wave = all_Waves(all_Waves(:,3)==nE,:);
        thr_wave = prctile(all_Waves(all_Waves(:,3)==nE,9),80);
        nE_wave = nE_wave(nE_wave(:,9)>=thr_wave,:);
        a = zeros(60,1);
        for nPr = 1:60
%             a(nPr) =  mean(nE_wave(nE_wave(:,2)==nPr,4));
              a(nPr) = length(nE_wave(nE_wave(:,2)==nPr,1));  
%             if isnan(a(nPr))==true
%                 a(nPr)=0;
%             end
        end
        
        
        
            b = probe_res(:,38); % sleepiness rating
            
            all_nTheta = [all_nTheta ; a];
            all_sleepiness =[all_sleepiness ; b];
%         b(a==0)=[];
%         a(a==0)=[];
        
        
    end
    
    all_nTheta = zscore(all_nTheta);
    all_nTheta = vec2mat(all_nTheta,60);
    all_sleepiness = zscore(all_sleepiness);
    all_sleepiness = vec2mat(all_sleepiness,60);
    
    for nE = 1:63
        a = all_nTheta(nE,:);
        b = all_sleepiness(nE,:);
        [rho] = corr(a(:),b(:),'type','spearman');
        all_corrs(nE,n)=rho;
        collect_all(nE)=rho;
    end
    
 
    subplot(5,4,n);
format_fig;
addpath(genpath(path_eeglab));
temp_topo=collect_all; % vector of 63 values
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
% caxis([0 2.5])
rmpath(genpath(path_eeglab));
colorbar;
caxis([min(temp_topo) max(temp_topo)])
xlabel(SubID);
    
end 

all_corrs(:,5)=[];


%% nTheta trend over time (probes)




clear all;
close all;
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

all_corrs = zeros(63,20);
all_nTheta=[];
alln_nTheta=[];
for n=1:length(bsl_files)
    if n==5;
        continue
    end
    
    
    
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);    
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    load([eeg_path filesep 'New folder' filesep 'wanderIM_twa2_' SubID])
    
    all_Waves_old=all_Waves;
    thrP2P=60;
    all_Waves(all_Waves(:,4)>thrP2P,:)=[];
        fprintf('... ... discarding %3.0f %% of waves (P2P thr: %guV)\n',mean(all_Waves_old(:,4)>thrP2P)*100,thrP2P)   
    
       all_nTheta=[];   
    for nE = 1:63    
        
        nTheta=[];
           nE_Waves = all_Waves(all_Waves(:,3)==nE,:);
           thr_wave = prctile(nE_Waves(:,4),80);
            for npr = 1:60
                npr_wave = nE_Waves(nE_Waves(:,2)==npr,:);
                nTheta = [nTheta length(npr_wave(npr_wave(:,4)>=thr_wave,2))];
                all_nTheta=[all_nTheta length(npr_wave(npr_wave(:,4)>=thr_wave,2))];    
            end
            xnpr = 1:60;
            
%             plot(xnpr,nTheta);
       
       end
       all_nTheta = zscore(all_nTheta);     
%        
%       subplot(5,4,n);
%       plot(nTheta)
%       axis([0 60 -3 3])
%       xlabel(SubID)
       
       
    
    alln_nTheta = [alln_nTheta; all_nTheta];
    
end

%% Topographies; at differing frequencies
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
figure;

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
        
    
    Fz_Waves = all_Waves(all_Waves(:,3)==2,:);
    thr_wave = prctile(Fz_Waves(:,4),80);
    halfw_samples = Fz_Waves(Fz_Waves(:,4)>=thr_wave,6) - Fz_Waves(Fz_Waves(:,4)>=thr_wave,5);
%     halfw_samples = Fz_Waves(:,6)-Fz_Waves(:,5);
    
    period = 1./(halfw_samples/500); % ^ Half-period in samples --> convert to Hz
   
    [N edges] = histcounts(period,0:0.2:10); % output 'edges' includes rightmost bin
    histc_record(n,:)= N;
    N=(N./(sum(N))*100); % convert to percentage
    % 20 participants, 5 by 4 subplotting
    collect_all = [collect_all ; N];
    subplot(5,4,n)
    
    xlabel(SubID);
    plot(edges(1:length(edges)-1),N);
    axis([0 10 0 10]);
    xlabel(SubID);
end    