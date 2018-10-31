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

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'lprobe_nfEEG_S3*.mat']);

%% loop across trials for baseline blocks
all_Waves_byProbes=[];
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
    
    left_freq(n)=SubjectInfo.FlickerR; % CAREFUL this is inverted
    right_freq(n)=SubjectInfo.FlickerL;
    
    param=[];
    param.method='fft'; % fast fourier transform
    param.mindist=1; % we want to be able to separate peaks separated by at least 1 Hz
    these_times=D.indsample(-20):D.indsample(0)-1;
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    
    % all_Waves
    % nSub nProbe nE P2Pamp posX negX WaveEnd MaxNegPeak MaxPosPeak
    % PaxsPosPeakAmp MaxDownSlope MaxUpSlope
    
    load([eeg_path filesep 'wanderIM_twa_' SubID])
    for nE=1:63
        thr_Wave1(n,nE)=prctile(all_Waves(all_Waves(:,3)==nE,4),80);
        thr_Wave2(n,nE)=prctile(all_Waves(all_Waves(:,3)==nE,4),90);
        
        num_Waves(n,nE)=sum(all_Waves(all_Waves(:,3)==nE,4)>thr_Wave2(n,nE));
        upSlope_Waves(n,nE)=mean(all_Waves(all_Waves(:,3)==nE & all_Waves(:,4)>thr_Wave2(n,nE),12));
        downSlope_Waves(n,nE)=mean(all_Waves(all_Waves(:,3)==nE & all_Waves(:,4)>thr_Wave2(n,nE),11));
        
        for npr=1:60
            tp_num_Waves=sum(all_Waves(all_Waves(:,2)==npr & all_Waves(:,3)==nE,4)>thr_Wave2(n,nE));
            tp_upSlope_Waves=mean(all_Waves(all_Waves(:,2)==npr & all_Waves(:,3)==nE & all_Waves(:,4)>thr_Wave2(n,nE),12));
            tp_downSlope_Waves=mean(all_Waves(all_Waves(:,2)==npr & all_Waves(:,3)==nE & all_Waves(:,4)>thr_Wave2(n,nE),11));
            
            all_Waves_byProbes=[all_Waves_byProbes ; [n npr probe_res(npr,[4 5 6 32 38]) nE probe_res(npr,32) tp_num_Waves tp_upSlope_Waves tp_downSlope_Waves]];
        end
    end
    
    %                 num_Waves(nE,npr)=length(cell2mat(twa_results.channels(nE).maxnegpkamp));
    %             amp_Waves(nE,npr)=mean(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))+abs(cell2mat(twa_results.channels(nE).maxpospkamp)));
end

%%
load(['../eeg/EasyCap64_layout']);

figure;
for np=1:2
    subplot(1,2,np); format_fig;
    if np==1
        temp_topo=mean(upSlope_Waves);
    elseif np==2
        temp_topo=mean(downSlope_Waves);
    end
    addpath(genpath(path_eeglab));
    topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off');
    rmpath(genpath(path_eeglab));
end

%%
figure;
for nt=1:2
    for nstate=1:3
        subplot(2,3,3*(nt-1)+nstate); format_fig;
        temp_topo=nan(1,63);
        for nE=1:63
        temp_topo(nE)=nanmean(all_Waves_byProbes(all_Waves_byProbes(:,[8])==nE & all_Waves_byProbes(:,[4])==nt & all_Waves_byProbes(:,[6])==nstate,9+np));
        end
        addpath(genpath(path_eeglab));
        topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off');
        cmap=colormap('hot'); cmap=flipud(cmap); colormap(cmap);
            caxis([0 5])
        rmpath(genpath(path_eeglab));
    end
end

%%
figure;
for np=1:3
    for nstate=1:3
        subplot(3,3,3*(np-1)+nstate); format_fig;
        temp_topo=nan(1,63);
        for nE=1:63
        temp_topo(nE)=nanmean(all_Waves_byProbes(all_Waves_byProbes(:,[8])==nE & all_Waves_byProbes(:,[6])==nstate,9+np));
        end
        addpath(genpath(path_eeglab));
        topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
        cmap=colormap('hot'); cmap=flipud(cmap); colormap(cmap);
        if np==1
            caxis([1 4])
        else
            caxis([0 10000])
        end
        rmpath(genpath(path_eeglab));
    end
end

%%
%[n npr probe_res(npr,[4 5 6 32 38]) nE probe_res(npr,32) tp_num_Waves tp_upSlope_Waves tp_downSlope_Waves]
tbl=array2table(all_Waves_byProbes,'VariableNames',{'SubID','nProbe','nBlock','Task','nTrial','State','Vig','Chan','State2','nWave','UpSlow','DownSlope'});
tbl(tbl.State==4,:)=[];
tbl.SubID=categorical(tbl.SubID);
tbl.nBlock=categorical(tbl.nBlock);
tbl.Task=categorical(tbl.Task);
tbl.State=categorical(tbl.State);
tbl.Chan=categorical(tbl.Chan);

mdl1= fitlme(tbl,'Vig~nBlock+(1|SubID)');
mdl2= fitlme(tbl,'nWave~nBlock+Task*Chan*State+(1|SubID)');
