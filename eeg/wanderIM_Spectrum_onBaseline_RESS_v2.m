%% Init
clear all;
close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
addpath(genpath(RESS_path))

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'lbasel_nfEEG_S3*.mat']);

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
    
    for nT=1:2
        these_times=D.indsample(0):D.indsample(30)-1;
        if nT==1
            temp_data=D(1:63,these_times,[1 2 5 6]); % D contains the data with channels * time * trials
        elseif nT==2
            temp_data=D(1:63,these_times,[3 4 7 8]); % D contains the data with channels * time * trials
        end
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
            
            baseline_logSNR_RESS(n,nT,nf,:,:)=logSNR_comp;
            baseline_Comp_RESS(n,nT,nf,:,:)=components;
            baseline_Maps_RESS(n,nT,nf,:)=maps;
        end
        [logSNR_dat, faxis_dat, logpow_dat]=get_logSNR(temp_data,D.fsample,param);
        baseline_logSNR_chan(n,nT,:,:,:)=logSNR_dat;
    end
end

%% Power spectrum and log SNR - do we have signal?
this_ch=match_str(D.chanlabels,'Oz');
[~,f1_idx]=findclosest(faxis_comp,6);
[~,f2_idx]=findclosest(faxis_comp,7.5);
[~,harm1_idx]=findclosest(faxis_comp,6*2);
[~,harm2_idx]=findclosest(faxis_comp,7.5*2);
[~,IM_idx]=findclosest(faxis_comp,6+7.5);

%%%%
names_blocks={'Face','FaceGap','Square','InvSquare'};
figure;
format_fig;
for nCond=1:4
    SNR_Method1(:,1,nCond)=squeeze(mean(baseline_logSNR_chan(:,:,this_ch,f1_idx,nCond),2));
    SNR_Method1(:,2,nCond)=squeeze(mean(baseline_logSNR_chan(:,:,this_ch,f2_idx,nCond),2));
    SNR_Method1(:,3,nCond)=squeeze(mean(baseline_logSNR_chan(:,:,this_ch,harm1_idx,nCond),2));
    SNR_Method1(:,4,nCond)=squeeze(mean(baseline_logSNR_chan(:,:,this_ch,harm2_idx,nCond),2));
    SNR_Method1(:,5,nCond)=squeeze(mean(baseline_logSNR_chan(:,:,this_ch,IM_idx,nCond),2));
    
    SNR_Method2(:,1,nCond)=squeeze(mean(baseline_logSNR_RESS(:,:,1,f1_idx,nCond),2));
    SNR_Method2(:,2,nCond)=squeeze(mean(baseline_logSNR_RESS(:,:,2,f2_idx,nCond),2));
    SNR_Method2(:,3,nCond)=squeeze(mean(baseline_logSNR_RESS(:,:,3,harm1_idx,nCond),2));
    SNR_Method2(:,4,nCond)=squeeze(mean(baseline_logSNR_RESS(:,:,4,harm2_idx,nCond),2));
    SNR_Method2(:,5,nCond)=squeeze(mean(baseline_logSNR_RESS(:,:,5,IM_idx,nCond),2));
end
for nCond=1:4
    subplot(2,2,nCond); format_fig;
    for nfreq=1:5
        if nCond==1 || nCond==3
            hb(1)=simpleBarPlot(nfreq-0.2,SNR_Method1(:,nfreq,nCond),'k',0.38,'r',[],2);
            hb(2)=simpleBarPlot(nfreq+0.2,SNR_Method2(:,nfreq,nCond),[0.5 0.5 0.5],0.38,'r',[],2);
        else
            hb(1)=simpleBarPlot(nfreq-0.2,SNR_Method1(:,nfreq,nCond),'k',0.38,'r',{0 SNR_Method1(:,nfreq,nCond-1) 0.05},2);
            hb(2)=simpleBarPlot(nfreq+0.2,SNR_Method2(:,nfreq,nCond),[0.5 0.5 0.5],0.38,'r',{0 SNR_Method2(:,nfreq,nCond-1) 0.05},2);
        end
    end
    title(names_blocks{nCond})
    if nCond==1, legend(hb,{'best chan','RESS'}); end
    set(gca,'XTick',1:5,'XTickLabel',{'f1','f2','2f1','2f2','IM'})
end

%% RESS component
% F1, F2, IM for left
% F1, F2, IM for right
% myFreqs=[6 7.5 12 15 13.5];
L_F_FMaps=[];
R_F_FMaps=[];
IM_FMaps=[];
L_F_DMaps=[];
R_F_DMaps=[];
IM_DMaps=[];

for n=1:length(bsl_files)
    if left_freq(n)==15
        L_F_FMaps=[L_F_FMaps ; real(squeeze(mean(baseline_Maps_RESS(n,1,2,:),2)))'];
        R_F_FMaps=[R_F_FMaps ; real(squeeze(mean(baseline_Maps_RESS(n,1,1,:),2)))'];
    else
        L_F_FMaps=[L_F_FMaps ; real(squeeze(mean(baseline_Maps_RESS(n,1,1,:),2)))'];
        R_F_FMaps=[R_F_FMaps ; real(squeeze(mean(baseline_Maps_RESS(n,1,2,:),2)))'];
    end
    IM_FMaps=[IM_FMaps ; real(squeeze(mean(baseline_Maps_RESS(n,1,5,:),2)))'];
    
        if left_freq(n)==15
        L_F_DMaps=[L_F_DMaps ; real(squeeze(mean(baseline_Maps_RESS(n,2,2,:),2)))'];
        R_F_DMaps=[R_F_DMaps ; real(squeeze(mean(baseline_Maps_RESS(n,2,1,:),2)))'];
    else
        L_F_DMaps=[L_F_DMaps ; real(squeeze(mean(baseline_Maps_RESS(n,2,1,:),2)))'];
        R_F_DMaps=[R_F_DMaps ; real(squeeze(mean(baseline_Maps_RESS(n,2,2,:),2)))'];
    end
    IM_DMaps=[IM_DMaps ; real(squeeze(mean(baseline_Maps_RESS(n,2,5,:),2)))'];
end

figure;
subplot(2,3,1); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(L_F_FMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('FACE - Right Fund')

subplot(2,3,2); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(R_F_FMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('FACE - Left Fund')

subplot(2,3,3); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(IM_FMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('FACE - IM: F2-F1')

subplot(2,3,4); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(L_F_DMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('DIGIT - Right Fund')

subplot(2,3,5); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(R_F_DMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('DIGIT - Left Fund')

subplot(2,3,6); format_fig;
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
topoplot(mean(IM_DMaps), layout.chaninfo,'style','both','whitebk','on','electrodes','off');
rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
title('DIGIT - IM: F2-F1')

%% Where are the tags on the scalp ?
% retrieve the channels position
% pos=D.coor2D';
% labels=D.chanlabels(1:63);
figure;
for nfre=1:5
    subplot(2,5,nfre); format_fig; % left fondamental
    temp_topo=real(squeeze(mean(baseline_Maps_RESS(:,1,nfre,:),1)));
    
    addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
    topoplot(temp_topo, layout.chaninfo,'style','both','whitebk','on','electrodes','off');
    rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
    
    caxis([-1 1]*max(max((mean(baseline_Maps_RESS(:,1,:,:),1))))/2)
    % title('Left F')
    
    subplot(2,5,5+nfre); format_fig; % left fondamental
    temp_topo=real(squeeze(mean(baseline_Maps_RESS(:,2,nfre,:),1)));
addpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
    topoplot(temp_topo, layout.chaninfo,'style','both','whitebk','on','electrodes','off');
    rmpath(genpath('/Users/tand0009/Work/local/eeglab14_1_2b/'));
    
    caxis([-1 1]*max(max((mean(baseline_Maps_RESS(:,2,:,:),1))))/2)
end

