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

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'lprobe_nfEEG_S3*.mat']);

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
    
    left_freq(n)=SubjectInfo.FlickerR; % CAREFUL this is inverted
    right_freq(n)=SubjectInfo.FlickerL;
    
    param=[];
    param.method='fft'; % fast fourier transform
    param.mindist=1; % we want to be able to separate peaks separated by at least 1 Hz
    these_times=D.indsample(-20):D.indsample(0)-1;
%     these_times=D.indsample(0)+1:D.indsample(20);
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    
    [logSNR, faxis, logpow]=get_logSNR(temp_data,D.fsample,param);
    
    onprobe_logSNR(n,:,:)=mean(logSNR,3);
    onprobe_logPow(n,:,:)=mean(logpow,3);
    
    % Split between Face and Digit
    idx_trials=find_trials(D.conditions,'FA');
    onprobe_logSNR_FA(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_FA(n,:,:)=mean(logpow(:,:,idx_trials),3);
    
    idx_trials=find_trials(D.conditions,'DT');
    onprobe_logSNR_DG(n,:,:)=mean(logSNR(:,:,idx_trials),3);
    onprobe_logPow_DG(n,:,:)=mean(logpow(:,:,idx_trials),3);
    
    for nT=1:2
        if nT==1
            idx_trials=find_trials(D.conditions,'FA');
        elseif nT==2
            idx_trials=find_trials(D.conditions,'DT');
        end
        temp_data=D(1:63,these_times,idx_trials); % D contains the data with channels * time * trials
        % RESS logSNR
        myFreqs=[6 7.5 12 15 13.5];
        for nf=1:length(myFreqs)
            paramRESS=[];
            paramRESS.peakfreq1=myFreqs(nf);
            paramRESS.peakwidt=0.2;
            paramRESS.neighfreq=1;
            paramRESS.neighwidt=1;
            paramRESS.fft_res=0.1; %0.1 default
            paramRESS.mindist=0.3; %0.5Hz default
            paramRESS.lowerfreq=1; %2Hz default
            [snrR, snrE, faxis2, maps, components]=get_logSNR_RESS(temp_data,D.fsample,paramRESS);
            
            param=[];
            param.method='fft'; % fast fourier transform
            param.mindist=1;
            RESS_comp(1,:,:)=components;
            [logSNR_comp, faxis_comp, logpow_comp]=get_logSNR(RESS_comp,D.fsample,param);
            
            onprobe_logSNR_RESS(n,nT,nf,:,:)=logSNR_comp;
            onprobe_Comp_RESS(n,nT,nf,:,:)=components;
            onprobe_Maps_RESS(n,nT,nf,:)=maps;
        end
        [logSNR_dat, faxis_dat, logpow_dat]=get_logSNR(temp_data,D.fsample,param);
        onprobe_logSNR_chan(n,nT,:,:,:)=logSNR_dat;
    end
end

%%
figure;
subplot(1,2,1)
plot(faxis,squeeze(mean(onprobe_logSNR(:,match_str(D.chanlabels,'Oz'),:),1)),'b')
xlim([1 30])

subplot(1,2,2)
plot(faxis,squeeze(nanmean(onprobe_logPow(:,match_str(D.chanlabels,'Oz'),:),1)),'r')
xlim([1 30])

%%
figure;
subplot(1,2,1)
plot(faxis,squeeze(mean(onprobe_logSNR_FA(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','b')
hold on;
plot(faxis,squeeze(mean(onprobe_logSNR_DG(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','r')
xlim([1 30])
legend({'Face','Digit'})

subplot(1,2,2)
plot(faxis,squeeze(mean(onprobe_logPow_FA(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','b')
hold on;
plot(faxis,squeeze(mean(onprobe_logPow_DG(:,match_str(D.chanlabels,'Oz'),:),1)),'Color','r')
xlim([1 30])
legend({'Face','Digit'})

%%
load('../BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
myFreq=6;

figure;
subplot(2,2,1); format_fig;
[closestf, idxclosest]=findclosest(faxis,myFreq);
temp_topo=squeeze(mean(onprobe_logSNR(:,:,idxclosest),1));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-2 2]); colorbar;

myFreq=7.5;

subplot(2,2,2); format_fig;
[closestf, idxclosest]=findclosest(faxis,myFreq);
temp_topo=squeeze(mean(onprobe_logSNR(:,:,idxclosest),1));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-2 2]); colorbar;

myFreq=7.5-6;

subplot(2,2,3); format_fig;
[closestf, idxclosest]=findclosest(faxis,myFreq);
temp_topo=squeeze(mean(onprobe_logSNR(:,:,idxclosest),1));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-2 2]); colorbar;

myFreq=7.5+6;

subplot(2,2,4); format_fig;
[closestf, idxclosest]=findclosest(faxis,myFreq);
temp_topo=squeeze(mean(onprobe_logSNR(:,:,idxclosest),1));
simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
caxis([-2 2]); colorbar;

%%
figure;
subplot(2,2,1); format_fig;

temp_topo=[];
for n=1:length(left_freq)
    [~,idxclosest]=findclosest(faxis,left_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(onprobe_logSNR(n,:,idxclosest),1));
end
simpleTopoPlot2(mean(temp_topo,1), pos', labels,0,[],0,lay,[]);
caxis([-2 2])
title('Left F')


subplot(2,2,2); format_fig;

temp_topo=[];
for n=1:length(right_freq)
    [~,idxclosest]=findclosest(faxis,right_freq(n)/2);
    temp_topo(n,:)=squeeze(mean(onprobe_logSNR(n,:,idxclosest),1));
end
simpleTopoPlot2(mean(temp_topo,1), pos', labels,0,[],0,lay,[]);
caxis([-2 2])
title('Right F')

subplot(2,2,3); format_fig;

temp_topo=[];
for n=1:length(left_freq)
    [~,idxclosest]=findclosest(faxis,left_freq(n));
    temp_topo(n,:)=squeeze(mean(onprobe_logSNR(n,:,idxclosest),1));
end
simpleTopoPlot2(mean(temp_topo,1), pos', labels,0,[],0,lay,[]);
caxis([-2 2])
title('Left 2F')


subplot(2,2,4); format_fig;

temp_topo=[];
for n=1:length(right_freq)
    [~,idxclosest]=findclosest(faxis,right_freq(n));
    temp_topo(n,:)=squeeze(mean(onprobe_logSNR(n,:,idxclosest),1));
end
simpleTopoPlot2(mean(temp_topo,1), pos', labels,0,[],0,lay,[]);
caxis([-2 2])
title('Right 2F')

%%
figure;
format_fig;

temp_topo=zscore(squeeze(mean(mean(onprobe_logPow(:,:,faxis>8 & faxis<11 & (faxis~=9 | faxis~=10.5)),3),1)));

simpleTopoPlot2(mean(temp_topo,1), pos', labels,0,[],0,lay,[]);
% caxis([-2 2])
title('Alpha')


%% Where are the tags on the scalp ?
% retrieve the channels position
load('../BrainVision_63ChLayout.mat') % the position are not ideal here - to be modified
% pos=D.coor2D';
% labels=D.chanlabels(1:63);
figure;
for nfre=1:5
    subplot(2,5,nfre); format_fig; % left fondamental
    temp_topo=(squeeze(mean(zscore(onprobe_Maps_RESS(:,1,nfre,:),[],4),1)));
    simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
    caxis([-1 1]*max(max((mean(zscore(onprobe_Maps_RESS(:,1,:,:),[],4),1)))))
    % title('Left F')
    
    subplot(2,5,5+nfre); format_fig; % left fondamental
    temp_topo=(squeeze(mean(zscore(onprobe_Maps_RESS(:,2,nfre,:),[],4),1)));
    simpleTopoPlot2(temp_topo, pos', labels,0,[],0,lay,[]);
    caxis([-1 1]*max(max((mean(zscore(onprobe_Maps_RESS(:,2,:,:),[],4),1)))))
end
