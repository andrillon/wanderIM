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
eeg_path2=[root_path filesep 'preproc_ica'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path2 filesep 'probe_infEEG_S3*.mat']);

prticle_Thr=90; % 80 or 90 or 95
LimFrqW=[1 4]; % [1 4] or [4 10]
AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
thisChannel=17;

%% loop across trials for baseline blocks
all_Waves_byProbes=[];
all_Waves_byProbes2=[];
all_len=[];
all_Waves_byProbes3=[];
nc=0;
for n=1:length(bsl_files)
    if n==5
        continue;
    end
    nc=nc+1;
    % load file with spm
    filename=bsl_files(n).name;
    D=spm_eeg_load([bsl_files(n).folder filesep filename]);
    these_times=D.indsample(-20):D.indsample(0)-1;
    temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    temp_data=temp_data-repmat(mean(temp_data([10 21],:,:),1),[size(temp_data,1) 1 1]);
    
    % load behavioural results
    SubID=bsl_files(n).name;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    fprintf('... processing subject %s\n',SubID)
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    
    left_freq(n)=SubjectInfo.FlickerR; % CAREFUL this is inverted
    right_freq(n)=SubjectInfo.FlickerL;
    
    %     param=[];
    %     param.method='fft'; % fast fourier transform
    %     param.mindist=1; % we want to be able to separate peaks separated by at least 1 Hz
    %     these_times=D.indsample(-20):D.indsample(0)-1;
    %     temp_data=D(1:63,these_times,:); % D contains the data with channels * time * trials
    
    % all_Waves
    % nSub nProbe nE P2Pamp negX posX WaveEnd MaxNegPeak MaxPosPeak
    % PaxsPosPeakAmp MaxDownSlope MaxUpSlope
    
%     load([eeg_path filesep 'wanderIM_twa3_' SubID])
    load([eeg_path2 filesep 'wanderIM_twa4_' SubID])
    
    for nE=1:63
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./500;
        temp_freq=1./temp_len;
        temp_abs=1./temp_len;
        temp_p2p=thisE_Waves(:,AmpCriterionIdx);
        
        %         all_len=[all_len ; [repmat([n nE],length(temp_len),1) temp_len temp_abs]];
        %         thr_Wave1(n,nE)=prctile(all_Waves(all_Waves(:,3)==nE,4),80);
        %         thr_Wave2(n,nE)=prctile(temp_p2p(temp_freq>5),80);
        thr_Wave2(nc,nE)=prctile(all_Waves((temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),AmpCriterionIdx),prticle_Thr);
        num_Waves(nc,nE)=sum(temp_p2p>thr_Wave2(nc,nE));
        %           len_Waves(n,nE)=sum(all_Waves(all_Waves(:,3)==nE,4)>thr_Wave2(n,nE));
        upSlope_Waves(nc,nE)=mean(thisE_Waves(temp_p2p>thr_Wave2(nc,nE),12));
        downSlope_Waves(nc,nE)=mean(thisE_Waves(temp_p2p>thr_Wave2(nc,nE),11));
        
        probe_res=[probe_res (probe_res(:,4)-1)*10+probe_res(:,1)];
        for npr=1:60
            thisE_Waves=all_Waves(all_Waves(:,3)==nE & all_Waves(:,2)==npr,:);
            temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./500;
            temp_freq=1./temp_len;
            temp_abs=1./temp_len;
            temp_p2p=thisE_Waves(:,4);
            
            tp_num_Start=thisE_Waves((temp_p2p>thr_Wave2(nc,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2))),5);
            tp_num_Waves=sum(temp_p2p>thr_Wave2(nc,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)));
            tp_P2P_Waves=nanmean(temp_p2p((temp_p2p>thr_Wave2(nc,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)))));
            tp_posAmp_Waves=nanmean(thisE_Waves((temp_p2p>thr_Wave2(nc,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2))),11));
            tp_negAmp_Waves=nanmean(thisE_Waves((temp_p2p>thr_Wave2(nc,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2))),9));
            tp_numE_Waves=length(unique(thisE_Waves((temp_p2p>thr_Wave2(nc,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2))),3)));
            %             tp_negAmp_Waves=nanmean(temp_p2p((temp_p2p>thr_Wave2(n,nE) & temp_freq>4)));
            
            tp_upSlope_Waves=nanmean(thisE_Waves(thisE_Waves(:,2)==npr & temp_p2p>thr_Wave2(nc,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),12));
            tp_downSlope_Waves=nanmean(thisE_Waves(thisE_Waves(:,2)==npr & temp_p2p>thr_Wave2(nc,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),11));
            
            these_probes=probe_res(probe_res(:,end)==npr,:);
            these_trials=test_res(test_res(:,1)==these_probes(4),:);
            this_pr_tridx=these_probes(1,6);
            last_pr_tridx=this_pr_tridx-21;
            temp_testres=these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,:);
            %             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            tcorr_go=nanmean(temp_testres(:,12));%/corr_go(n,these_probes(npr,5));
            tcorr_nogo=nanmean(temp_testres(:,11));%/corr_nogo(n,these_probes(npr,5));
            num_go=sum(~isnan(temp_testres(:,12)));
            num_nogo=sum(~isnan(temp_testres(:,11)));
            rt_go=nanmean(temp_testres(~isnan(temp_testres(:,12)),10)-temp_testres(~isnan(temp_testres(:,12)),8));
            rt_nogo=nanmean(temp_testres(~isnan(temp_testres(:,11)),10)-temp_testres(~isnan(temp_testres(:,11)),8));
            
            all_Waves_byProbes=[all_Waves_byProbes ; [n npr probe_res(npr,[4 5 6 32 38]) nE probe_res(npr,32) tp_num_Waves tp_upSlope_Waves tp_downSlope_Waves tp_P2P_Waves tp_posAmp_Waves tp_negAmp_Waves tp_numE_Waves rt_go tcorr_go tcorr_nogo]];
            all_Waves_byProbes2=[all_Waves_byProbes2 ; [n npr probe_res(npr,[4 5 6 32 33 34 35 36 37 38]) nE probe_res(npr,32) tp_num_Waves tp_upSlope_Waves tp_downSlope_Waves tp_P2P_Waves tp_posAmp_Waves]];
            all_Waves_byProbes3=[all_Waves_byProbes3 ; [repmat([n npr probe_res(npr,[4 5 6 32 33 34 35 36 37 38]) nE probe_res(npr,32)],length(tp_num_Start),1) tp_num_Start]];
        end
        
        if nE==24 %Cz=24
            
            thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
            temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./500;
            temp_freq=1./temp_len;
            temp_abs=1./temp_len;
            temp_p2p=thisE_Waves(:,4);
            onset_Waves=thisE_Waves(temp_p2p>thr_Wave2(nc,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),5);
            probes_Waves=thisE_Waves(temp_p2p>thr_Wave2(nc,nE) & (temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),2);
            
            temp_ERP_Waves=[];
            for nW=1:length(onset_Waves)
                if onset_Waves(nW)-1*D.fsample<1 || onset_Waves(nW)+2*D.fsample>size(temp_data,2)-1
                    temp_ERP_Waves(nW,:)=nan(1,length((-1*D.fsample:2*D.fsample)));
                    continue;
                end
                temp=temp_data(nE,onset_Waves(nW)+(-1*D.fsample:2*D.fsample),probes_Waves(nW));
                if max(abs(temp))>200
                    temp_ERP_Waves(nW,:)=nan(1,length((-1*D.fsample:2*D.fsample)));
                    continue;
                end
                temp_ERP_Waves(nW,:)=temp;
                temp_ERP_Waves(nW,:)=temp_ERP_Waves(nW,:)-mean(temp_ERP_Waves(nW,:));
            end
            numERP_Waves(nc)=size(temp_ERP_Waves,1);
            if numERP_Waves(nc)<50
                ERP_Waves(nc,:)=nan;
            else
                ERP_Waves(nc,:)=nanmean(temp_ERP_Waves);
            end
        end
    end
    
    %                 num_Waves(nE,npr)=length(cell2mat(twa_results.channels(nE).maxnegpkamp));
    %             amp_Waves(nE,npr)=mean(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))+abs(cell2mat(twa_results.channels(nE).maxpospkamp)));
end

%%
figure;
format_fig;
temp_topo=nanmean(num_Waves)/20/60;
addpath(genpath(path_eeglab));
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','on');
cmap=colormap('hot'); cmap=flipud(cmap); colormap(cmap);
caxis([0 0.6])
rmpath(genpath(path_eeglab));
cb=colorbar; set(cb,'FontSize',24)

figure; format_fig
simpleTplot(-1:1/D.fsample:2,squeeze((ERP_Waves(:,:))),0,'k',[0],'-',0.5,1,2,[],2);
xlim([-0.3 0.6]);
xlabel('Time from wave onset (s)')
ylabel('Amplitude (\muV)')
title('20% higher amplitude')

%%
%[n npr probe_res(npr,[4 5 6 32 38]) nE probe_res(npr,32) tp_num_Waves tp_upSlope_Waves tp_downSlope_Waves]
tbl=array2table(all_Waves_byProbes,'VariableNames',{'SubID','nProbe','nBlock','Task','nTrial','State','Vig','Chan','State2','nWave','UpSlope','DownSlope','P2P','PosAmp','NegAmp','numbE','RT','GO','NOGO'});
tbl(tbl.State==4,:)=[];
tbl.SubID=categorical(tbl.SubID);
% tbl.nBlock=categorical(tbl.nBlock);
tbl.Task=categorical(tbl.Task);
tbl.State=categorical(tbl.State);
% tbl.Chan=categorical(tbl.Chan);
tbl.pWave=tbl.nWave;
tbl.pWave(tbl.nWave>0)=1;
tbl.pWave(tbl.nWave<=0)=0;

for nE=1:63
    [r_NOGO(nE) pV_NOGO(nE)]=corr(tbl.nWave(tbl.Chan==nE),1-tbl.NOGO(tbl.Chan==nE),'type','pearson');
    [r_MISS(nE) pV_MISS(nE)]=corr(tbl.nWave(tbl.Chan==nE),1-tbl.GO(tbl.Chan==nE),'type','pearson');
    temp=tbl.RT(tbl.Chan==nE);
    temp(temp<0.2)=nan;
    [r_RT(nE) pV_RT(nE)]=corr(tbl.nWave(tbl.Chan==nE),temp,'type','pearson','rows','pairwise');
end
% filename = '/Users/tand0009/Data/WanderIM/preproc_eeg/wanderIM_table_localsleep_results3.txt';
% res=import_localSleepRes_fromR(filename);

figure;
subplot(1,3,1)
addpath(genpath(path_eeglab));
temp_all=max([r_MISS; r_NOGO]);
temp_topo=r_MISS;
tempplotpV=pV_MISS;
stat_thr=fdr(tempplotpV,0.01);
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{find(tempplotpV<=stat_thr),'.','w',24,2});
rmpath(genpath(path_eeglab));
cb=colorbar; colormap('parula'); set(cb,'FontSize',20); caxis([-1 1]*0.2)
format_fig;
title('Misses')

subplot(1,3,2)
addpath(genpath(path_eeglab));
temp_topo=r_NOGO;
tempplotpV=pV_NOGO;
stat_thr=fdr(tempplotpV,0.01);
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{find(tempplotpV<=stat_thr),'.','w',24,2});
rmpath(genpath(path_eeglab));
cb=colorbar; colormap('parula'); set(cb,'FontSize',20); caxis([-1 1]*0.2)
format_fig;
title('False Alarms')

subplot(1,3,3)
addpath(genpath(path_eeglab));
temp_topo=r_RT;
tempplotpV=pV_RT;
stat_thr=fdr(tempplotpV,0.01);
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{find(tempplotpV<=stat_thr),'.','w',24,2});
rmpath(genpath(path_eeglab));
cb=colorbar; colormap('parula'); set(cb,'FontSize',20); caxis([-0.2 0.2])
format_fig;title('RT')

% %%
% thisChannel=17;
% figure; format_fig;
% simpleCorPlotsetbin(tbl.nWave(tbl.Chan==thisChannel),tbl.RT(tbl.Chan==thisChannel),0:1:8,{'o','k',[1 1 1]*0.7,72,3});
% xlabel('Number of Waves')
% ylabel('RT')
% set(gca,'XTick',0:1:8,'XTickLabel',0:8)
% xlim([-0.5 8.5])

%%
tbl=array2table(all_Waves_byProbes,'VariableNames',{'SubID','nProbe','nBlock','Task','nTrial','State','Vig','Chan','State2','nWave','UpSlope','DownSlope','P2P','PosAmp','NegAmp','numbE','RT','GO','NOGO'});

%%
thisCh=2;
allwaves_perProbe=[];
allvig_perProbe=[];
myS=unique(tbl.SubID);
for nS=1:length(myS)
    for nB=1:6
        waves_perProbe(nS,nB)=nansum(tbl.nWave(tbl.SubID==myS(nS) & tbl.nBlock==nB & tbl.Chan==thisCh));
        if isempty(unique(tbl.Vig(tbl.SubID==myS(nS) & tbl.nBlock==nB & tbl.Chan==thisCh)))
            vig_perProbe(nS,nB)=nan;
        else
            vig_perProbe(nS,nB)=nanmean(tbl.Vig(tbl.SubID==myS(nS) & tbl.nBlock==nB & tbl.Chan==thisCh));
        end
    end
    [r_waveVig(nS), pV_waveVig(nS)]=corr(waves_perProbe(nS,:)',vig_perProbe(nS,:)');
    allwaves_perProbe=[allwaves_perProbe waves_perProbe(nS,:)];
    allvig_perProbe=[allvig_perProbe (vig_perProbe(nS,:))];
%     waves_perProbe(nS,:)=minmax(waves_perProbe(nS,:));
%     vig_perProbe(nS,:)=minmax(vig_perProbe(nS,:));
end

figure; format_fig; hold on;
tempA=2*waves_perProbe./repmat(waves_perProbe(:,1),1,size(waves_perProbe,2))-1; 
plot(mean(tempA),'LineWidth',3,'Color','b')
% plot(mean(tempA)+sem(tempA),'LineWidth',1,'Color','r','LineStyle','--')
% plot(mean(tempA)-sem(tempA),'LineWidth',1,'Color','r','LineStyle','--')
hold on
tempB=vig_perProbe./repmat(vig_perProbe(:,1),1,size(vig_perProbe,2));
plot(mean(tempB),'LineWidth',3,'Color','k')
% plot(mean(tempB)+sem(tempB),'LineWidth',1,'Color','k','LineStyle','--')
% plot(mean(tempB)-sem(tempB),'LineWidth',1,'Color','k','LineStyle','--')
xlabel('Block'); format_fig;
% plot(mean(waves_perProbe))
% hold on
% plot(mean(vig_perProbe))
ylim([1 1.6])

%%
figure; format_fig; hold on;
tempA=waves_perProbe./repmat(waves_perProbe(:,1),1,size(waves_perProbe,2)); 
tempB=vig_perProbe./repmat(vig_perProbe(:,1),1,size(vig_perProbe,2));
[axout,h1,h2]=plotyy(1:6,mean(tempA),1:6,mean(tempB));
set(h1,'Color','k','LineWidth',3)
set(h2,'Color','r','LineWidth',3)
set(axout(1),'YLim',[1 1.3])
set(axout(2),'YLim',[1 1.6])
[axout,h1,h2]=plotyy(1:6,mean(tempA)+sem(tempA),1:6,mean(tempB)+sem(tempB));
set(h1,'Color','k','LineWidth',1,'LineStyle','--')
set(h2,'Color','r','LineWidth',1,'LineStyle','--')
set(axout(1),'YLim',[1 1.3])
set(axout(2),'YLim',[1 1.6])
[axout,h1,h2]=plotyy(1:6,mean(tempA)-sem(tempA),1:6,mean(tempB)-sem(tempB));
set(h1,'Color','k','LineWidth',1,'LineStyle','--')
set(h2,'Color','r','LineWidth',1,'LineStyle','--')
set(axout(1),'YLim',[1 1.3])
set(axout(2),'YLim',[1 1.6])
xlabel('Block'); format_fig;
%%
figure;
for nState=1:3
    subplot(1,3,nState); format_fig;
    temp_topo=nan(1,63);
    for nE=1:63
        temp_topo(nE)=nanmean(tbl.nWave(tbl.State==nState & tbl.Chan==nE));
    end
    addpath(genpath(path_eeglab));
    topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off');
    cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
    caxis([0 5])
    rmpath(genpath(path_eeglab));
%     cb=colorbar; 
    set(cb,'FontSize',24)
    
%     if nState~=1
%           subplot(2,3,nState+3); format_fig;
%   temp_topo=nan(1,63);
%         for nE=1:63
%             temp_topo(nE)=nanmean(tbl.nWave(tbl.State==nState & tbl.Chan==nE))-nanmean(tbl.nWave(tbl.State==1 & tbl.Chan==nE));
%         end
%         addpath(genpath(path_eeglab));
%         topoplot((temp_topo), layout.chaninfo,'style','map','whitebk','on','electrodes','on');
%         cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
% %         caxis([-1 1]*3.5)
%         rmpath(genpath(path_eeglab));
%         cb=colorbar; set(cb,'FontSize',24)
%     end
end

%%
figure;
for nState=2:3
    subplot(1,2,nState-1); format_fig;
    temp_topo=nan(1,63);
    for nE=1:63
        temp_topo(nE)=nanmean(tbl.P2P(tbl.State==nState & tbl.Chan==nE))-nanmean(tbl.P2P(tbl.State==1 & tbl.Chan==nE));
    end
    addpath(genpath(path_eeglab));
    topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off');
    cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
    caxis([-1 1]*10)
    rmpath(genpath(path_eeglab));
%     cb=colorbar; 
    set(cb,'FontSize',24)
    
%     if nState~=1
%           subplot(2,3,nState+3); format_fig;
%   temp_topo=nan(1,63);
%         for nE=1:63
%             temp_topo(nE)=nanmean(tbl.nWave(tbl.State==nState & tbl.Chan==nE))-nanmean(tbl.nWave(tbl.State==1 & tbl.Chan==nE));
%         end
%         addpath(genpath(path_eeglab));
%         topoplot((temp_topo), layout.chaninfo,'style','map','whitebk','on','electrodes','on');
%         cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
% %         caxis([-1 1]*3.5)
%         rmpath(genpath(path_eeglab));
%         cb=colorbar; set(cb,'FontSize',24)
%     end
end

%%
myS=unique(tbl.SubID);
figure;
for nState=1:3
    subplot(1,3,nState); format_fig;
    temp_topo=nan(length(myS),63);
    for nS=1:length(myS)
        for nE=1:63
            temp_topo(nS,nE)=nanmean(tbl.nWave(tbl.SubID==myS(nS) & tbl.State==nState & tbl.Chan==nE));
        end
        temp_ntopo(nS)=nansum(tbl.SubID==myS(nS) & tbl.State==nState & tbl.Chan==1);
    end
    addpath(genpath(path_eeglab));
    av_topo=nansum(temp_topo.*temp_ntopo',1)./sum(temp_ntopo);
    topoplot(av_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off');
    cmap=colormap('parula'); %cmap=flipud(cmap); colormap(cmap);
    caxis([0 5])
    rmpath(genpath(path_eeglab));
%     cb=colorbar; 
%     set(cb,'FontSize',24)
    
end