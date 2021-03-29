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
eeg_path=[root_path filesep 'preproc_ica'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_infEEG_S3*.mat']);


%% loop across trials for baseline blocks
prticle_Thr=90; % 80 or 90 or 95
LimFrqW=[1 4]; % [1 4] or [4 10]
AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
% fixThr=75;
fixThr=[];
art_ampl=150;
max_posampl=75;
max_Freq=7;
frontalElecs=[1 32 33 60];
WavesAndVig=[];
nSc=0;

elecs_ERP=[2, 24, 13, 17];
for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    SubID=filename;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    if ~ismember(SubID,GoodSubID)
        continue;
    end
    
    D=spm_eeg_load([bsl_files(n).folder filesep filename]);
    these_times=D.indsample(D.time);
    temp_data=D(1:63,:,:); % D contains the data with channels * time * trials
    temp_data=temp_data-repmat(mean(temp_data([10 21],:,:),1),[size(temp_data,1) 1 1]);
    
    % load behavioural results
    if exist([eeg_path filesep 'wanderIM_twa5_noica_bigwnd_' SubID '.mat'])~=0
        fprintf('... load local sleep detection for subject %s\n',SubID)
        
        load([eeg_path filesep 'wanderIM_twa5_noica_bigwnd_' SubID]);
        %         load([eeg_path2 filesep 'wanderIM_twa4_' SubID]);
    else
        fprintf('... load local sleep detection subject %s DOES NOT EXIST\n',SubID)
        continue;
    end
    
    nSc=nSc+1;
    GoodSudID{nSc}=SubID;
    if AmpCriterionIdx==9
        all_Waves(:,AmpCriterionIdx)=-all_Waves(:,AmpCriterionIdx);
    end
    %     hdr=ft_read_header([eeg_path filesep filename]);
    
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    
    all_Waves(all_Waves(:,7)>0,:)=[];
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./D.fsample);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,AmpCriterionIdx)>art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>max_posampl | all_Waves(:,14)>art_ampl| abs(all_Waves(:,15))>art_ampl)*100)
    all_Waves(all_freq>max_Freq | all_Waves(:,AmpCriterionIdx)>art_ampl | all_Waves(:,11)>max_posampl| all_Waves(:,14)>art_ampl| abs(all_Waves(:,15))>art_ampl,:)=[];
    
    all_Wpos=[];
    slow_Waves=[];
    vicinity_matrix2(nSc,:)=zeros(1,63);
    for nE=1:63
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./D.fsample;
        temp_freq=1./temp_len;
        temp_abs=1./temp_len;
        temp_p2p=thisE_Waves(:,AmpCriterionIdx);
        temp_p2p2=thisE_Waves(:,4);
        
        
        if ~isempty(fixThr)
            thr_Wave(nSc,nE)=fixThr;
        else
            %             thr_Wave2(n,nE)=prctile(all_Waves((temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),AmpCriterionIdx),prticle_Thr);
            thr_Wave(nSc,nE)=prctile(thisE_Waves(:,AmpCriterionIdx),prticle_Thr);
        end
        dens_Waves(nSc,nE)=sum(temp_p2p>thr_Wave(nSc,nE));
        %         dens_Waves(n,nE)=sum(temp_p2p>thr_Wave(n,nE) & temp_p2p2>75)/(hdr.nSamples/D.fsample);
        %           len_Waves(n,nE)=sum(all_Waves(all_Waves(:,3)==nE,4)>thr_Wave2(n,nE));
        upSlope_Waves(nSc,nE)=mean(thisE_Waves(temp_p2p>thr_Wave(nSc,nE),12));
        downSlope_Waves(nSc,nE)=mean(thisE_Waves(temp_p2p>thr_Wave(nSc,nE),11));
        
        if ismember(nE,elecs_ERP) %Cz=24
            sel_Waves=thisE_Waves(temp_p2p>thr_Wave(nSc,nE),:);
            onset_Waves=sel_Waves(:,5)/D.fsample;
            temp_ERP_Waves=[];
            for nW=1:length(onset_Waves)
                if onset_Waves(nW)-1<-20 || onset_Waves(nW)+2>0
                    temp_ERP_Waves(nW,:)=nan(1,length((-1*D.fsample:2*D.fsample)));
                    continue;
                end
                temp=temp_data(nE,D.indsample(onset_Waves(nW))+(-1*D.fsample:2*D.fsample),sel_Waves(nW,2));
                if max(abs(temp))>200
                    temp_ERP_Waves(nW,:)=nan(1,length((-1*D.fsample:2*D.fsample)));
                    continue;
                end
                temp_ERP_Waves(nW,:)=temp;
                temp_ERP_Waves(nW,:)=temp_ERP_Waves(nW,:)-mean(temp_ERP_Waves(nW,:));
            end
            numERP_Waves{nE}(nSc)=size(temp_ERP_Waves,1);
            if numERP_Waves{nE}(nSc)<50
                ERP_Waves{nE}(nSc,:)=nan;
            else
                ERP_Waves{nE}(nSc,:)=nanmean(temp_ERP_Waves);
            end
        end
        
        
        % concomittant
        chanA_waves=thisE_Waves(temp_p2p>thr_Wave(nSc,nE),:);
        for nE2=1:63
            if nE2==nE || isempty(chanA_waves)
                vicinity_matrix(nSc,nE,nE2)=nan;
                continue;
            end
            % chanel B
            thisE2_Waves=all_Waves(all_Waves(:,3)==nE2,:);
            temp2_len=abs((thisE2_Waves(:,5)-thisE2_Waves(:,7)))./D.fsample;
            temp2_freq=1./temp2_len;
            temp2_abs=1./temp2_len;
            temp2_p2p=thisE2_Waves(:,AmpCriterionIdx);
            chanB_waves=thisE2_Waves(temp2_p2p>thr_Wave(nSc,nE),:);
            all_posW_B=chanB_waves(:,5);
            neigh=nan(1,size(chanA_waves,1));
            neigh2=nan(1,size(chanA_waves,1));
            for nW=1:size(chanA_waves,1)
                posW=chanA_waves(nW,5);
                neigh(nW)=~isempty(find(all_posW_B>posW-0.5*D.fsample & all_posW_B<posW+0.5*D.fsample & chanB_waves(:,2)==chanA_waves(nW,2)));
            end
            vicinity_matrix(nSc,nE,nE2)=nanmean(neigh);
        end
        % chanel ALL
        thisE2_Waves=all_Waves(all_Waves(:,3)~=nE,:);
        temp2_len=abs((thisE2_Waves(:,5)-thisE2_Waves(:,7)))./D.fsample;
        temp2_freq=1./temp2_len;
        temp2_abs=1./temp2_len;
        temp2_p2p=thisE2_Waves(:,AmpCriterionIdx);
        chanB_waves=thisE2_Waves(temp2_p2p>thr_Wave(nSc,nE),:);
        all_posW_B=chanB_waves(:,5);
        neigh2=nan(1,size(chanA_waves,1));
        for nW=1:size(chanA_waves,1)
            posW=chanA_waves(nW,5);
            neigh2(nW)=length(unique(chanB_waves(find(all_posW_B>posW-0.5*D.fsample & all_posW_B<posW+0.5*D.fsample),3)));
        end
        vicinity_matrix2(n,nE)=nanmean(neigh2);
        
        temp=chanA_waves(:,5);
        temp(find(diff(temp)/D.fsample<0.5)+1)=[];
        %         all_Wpos=[all_Wpos ; [temp nE*ones(length(temp),1)]];
        VigWaves=nan(size(chanA_waves,1),1);
        for nW=1:size(chanA_waves,1)
            VigWaves(nW)=probe_res(chanA_waves(nW,2),38);
        end
        slow_Waves=[slow_Waves ; [chanA_waves neigh2' VigWaves nanzscore(VigWaves)]];
        
        
    end
    
    Vig=probe_res(:,38);
    zVig=nanzscore(probe_res(:,38));
    for nPr=1:60
        for nE=1:63
            
            WavesAndVig=[WavesAndVig ; [n nPr probe_res(nPr,1) probe_res(nPr,4) probe_res(nPr,5) nE sum(slow_Waves(:,2)==nPr & slow_Waves(:,3)==nE)/(20) Vig(nPr) zVig(nPr)]];
        end
    end
    %     save([eeg_path filesep 'wanderIM_cont_twa2_' SubID '_selected'],'slow_Waves');
    
    %     [~,idx] = sort(all_Wpos(:,1));
    %     all_Wpos=all_Wpos(idx,:);
    %     diffPos=[NaN ; diff(all_Wpos(:,1)/D.fsample)];
    %     breakSeq=find(diffPos>0.5);
    %     begS=1;
    %     endS=breakSeq(1)-1;
    %     begE=all_Wpos(1,2);
    %     endE=all_Wpos(endS,2);
    %     begP=all_Wpos(1,1);
    %     endP=all_Wpos(endS,1);
    %     temp_traj=[begS endS begE endE begP endP (endP-begP)/D.fsample length(unique(all_Wpos(1:endS,2)))];
    %     for nSeq=2:length(breakSeq)
    %         begS=breakSeq(nSeq-1);
    %         endS=breakSeq(nSeq)-1;
    %         begE=all_Wpos(begS,2);
    %         endE=all_Wpos(endS,2);
    %         begP=all_Wpos(begS,1);
    %         endP=all_Wpos(endS,1);
    %         temp_traj=[temp_traj ; [begS endS begE endE begP endP (endP-begP)/D.fsample length(unique(all_Wpos(begS:endS,2)))]];
    %     end
    %     temp_traj(temp_traj(:,end-1)>1,:)=[];
end

%%
WavesAndVig(ismember(WavesAndVig(:,3),[10 21]),:)=[];
table=array2table(WavesAndVig,'VariableNames',{'SubID','nPr','nBlock','nPr2','Task','Elec','SWdens','Vig','zVig'});
table.Elec=categorical(table.Elec);
table.Task=categorical(table.Task);

%
% dens_Waves(:,[10 21])=NaN; %Mastoids
% thr_Wave(:,[10 21])=NaN; %Mastoids
% vicinity_matrix(:,[10 21],:)=NaN;
% vicinity_matrix(:,:,[10 21])=NaN;
% vicinity_matrix2(:,[10 21])=NaN;
%%
figure;
% subplot(1,3,1)
addpath((path_fieldtrip)); ft_defaults;
temp_topo=mean(dens_Waves,1); %in seconds but for all probes (60) so equivalent of in minutes
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'on',[],0,1);
colorbar; caxis([0 9])
title('SW density /min')
format_fig;
%export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo.fig'])
%export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo.eps'],'-r 300')
% 
% %%
% figure;
% % subplot(1,3,1)
% addpath((path_fieldtrip)); ft_defaults;
% temp_topo=squeeze(nanmean(vicinity_matrix2,1));
% simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'no',[],0,0);
% 

%%
allLabels=layout.label;
load(path_PsychFTlayout);
figure; set(gcf,'Position',[185         206        1394         756]);
myLabels={'Fz','Cz','Pz','Oz'};
for nE=1:4
    myE=find(ismember(layout.label,myLabels{nE}));
    subplot(1,4,nE)
    temp_topo=squeeze(nanmean(vicinity_matrix(:,myE,:),1));
    temp_topo(myE)=NaN;
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    hold on;
    scatter(layout.pos(match_str(layout.label,myLabels{nE}),1),...
        layout.pos(match_str(layout.label,myLabels{nE}),2),'Marker','o','MarkerFaceColor','r',...
        'MarkerEdgeColor','r','SizeData',72);
    colorbar('Location','EastOutside');
    caxis([0 0.55]);
    title(myLabels{nE}); format_fig
end
%export_fig([path_fig filesep 'LocalSleep_SWppties_NeighWaves.fig'])
%export_fig([path_fig filesep 'LocalSleep_SWppties_NeighWaves.eps'],'-r 300')

%%
load([root_path filesep 'SWsleepComp' filesep 'ERP_SW_Sleep.mat']);

figure; hold on;
subplot(1,2,1); hold on;
timeP1=-1:1/500:2;
plot(timeP1,nanmean(ERP_Waves{elecs_ERP(2)}),'LineWidth',4,'Color',[1 0 0]*0.8);
xlim([-.5 1])
xlabel('Time fron onset (s)')
ylabel('Voltage (\muV)')
format_fig;

subplot(1,2,2); hold on;
plot(timeP,squeeze(nanmean(allSO(countSO>16,:,65),1)),'LineWidth',4,'Color',[0 0 1]*0.8);
xlim([-.5 1])
xlabel('Time fron onset (s)')
ylabel('Voltage (\muV)')
format_fig;

%export_fig([path_fig filesep 'LocalSleep_SWppties_Waveform.fig'])
%export_fig([path_fig filesep 'LocalSleep_SWppties_Waveform.eps'],'-r 300')


 

%%
temp_topo=squeeze(nanmean(countSObyElec(countSO>16,:),1));
cfg=[];
cfg.layout = 'egi64_GSN_HydroCel_vTA.sfp';
layout=ft_prepare_layout(cfg);

    figure; set(gcf,'Position',[ 562   669   325   316]);
% subplot(1,3,1)
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('SW density')
format_fig;

  hb=colorbar; 
    set(hb,'YTick',[]);
    
% export_fig([path_fig filesep 'LocalSleep_SWppties_SleepTopo.fig'])
% export_fig([path_fig filesep 'LocalSleep_SWppties_SleepTopo.eps'],'-r 300')

%%
motorERP=load('/Users/tand0009/Data/WanderIM/ERPComp/wanderIM_motorERPs');
visualERP=load('/Users/tand0009/Data/WanderIM/ERPComp/wanderIM_visualERPs');
load('/Users/tand0009/Data/WanderIM/ERPComp/wanderIM_swERPs');

figure; set(gcf,'Position',[231         473        1239         483]);
elecs_ERP_labels={'Fz','Cz','Pz','Oz'};
cmap=cbrewer('qual','Set1',3);
for nEl=2 %1:length(elecs_ERP)
    subplot(1,3,1);
    hold on;
    xTime1=-1:1/500:2;
    my_ERP=(ERP_Waves{elecs_ERP(nEl)});
    [~, hplot1]=simpleTplot(xTime1,my_ERP,0,cmap(1,:),0,'-',0.5,1,10,1,1);

    xlim([-.2 .6])
    ylim([-14 6])
  xlabel('Time (s)')
    ylabel('Voltage (\muV)')
    format_fig;
    line(xlim,[0 0],'Color',[1 1 1]*0.7,'LineStyle','--');
         line([0 0],ylim,'Color',[1 1 1]*0.7,'LineStyle','--');
  legend(hplot1,'SW')
    
    subplot(1,3,2);
    hold on;
    xTime2=visualERP.xTime;
    this_ch=match_str(visualERP.layout.label,elecs_ERP_labels{nEl});
    my_ERP=squeeze(visualERP.STD_erp_D(:,this_ch,:));
%     plot(xTime2,(my_ERP),'LineWidth',4,'Color',[0 1 0]*0.9);
        [~, hplot2(1)]=simpleTplot(xTime2,my_ERP,0,cmap(2,:),0,'-',0.5,1,10,1,1);

    this_ch=match_str(visualERP.layout.label,elecs_ERP_labels{nEl});
    my_ERP=squeeze(visualERP.STD_erp_F(:,this_ch,:));
%     plot(xTime2,(my_ERP),'LineWidth',4,'Color',[0 1 0]*0.7);
        [~, hplot2(2)]=simpleTplot(xTime2,my_ERP,0,cmap(2,:),0,':',0.5,1,10,1,1);

    xlim([-.2 .6])
      ylim([-14 6])
  xlabel('Time (s)')
    ylabel('Voltage (\muV)')
    format_fig;
        line(xlim,[0 0],'Color',[1 1 1]*0.7,'LineStyle','--');
       line([0 0],ylim,'Color',[1 1 1]*0.7,'LineStyle','--');
    legend(hplot2,{'Digits','Faces'})

    subplot(1,3,3);
    hold on;
    xTime3=motorERP.xTime;
    this_ch=match_str(motorERP.layout.label,elecs_ERP_labels{nEl});
    my_ERP=squeeze(motorERP.Corr_erp(:,this_ch,:));
%     plot(xTime3,(my_ERP),'LineWidth',4,'Color',[0 0 1]*0.8);
        [~, hplot3]=simpleTplot(xTime3,my_ERP,0,cmap(3,:),0,'-',0.5,1,10,1,1);

        xlim([-.4 .4])
      ylim([-14 6])
  xlabel('Time (s)')
    ylabel('Voltage (\muV)')
    format_fig;
      line(xlim,[0 0],'Color',[1 1 1]*0.7,'LineStyle','--');
       line([0 0],ylim,'Color',[1 1 1]*0.7,'LineStyle','--');
     legend(hplot3,{'Go Resp.'})

end

export_fig([path_fig filesep 'LocalSleep_compERPs_VisualMotor.fig'])
export_fig([path_fig filesep 'LocalSleep_compERPs_VisualMotor.eps'],'-r 300')

%% Topographies
% figure;
% [~,index1]=findclosest(xTime1,0.1);
% temp_topo=squeeze(mean(visualERP.DGT_erp(:,:,index1),1));
% simpleTopoPlot_ft(temp_topo', visualERP.layout,'on',[],0,1);
% title('SW')
% format_fig;
% hb=colorbar;
% 
figure;
temp_topo=squeeze(mean(mean((visualERP.STD_erp_D(:,:,xTime2==0.2)),3),1));
simpleTopoPlot_ft(temp_topo', visualERP.layout,'on',[],0,1);
title('digits')
format_fig;
hb=colorbar;
caxis([-1 1]*2)
set(hb,'Position',[0.85 0.6 0.04 0.3])
export_fig([path_fig filesep 'LocalSleep_compERPs_DigitTopo.eps'],'-r 300')

% set(hb,'YTick',[]);
    
figure;
temp_topo=squeeze(mean(mean((visualERP.STD_erp_F(:,:,xTime2==0.2)),3),1));
simpleTopoPlot_ft(temp_topo', visualERP.layout,'on',[],0,1);
title('faces')
format_fig;
hb=colorbar;
caxis([-1 1]*2)
set(hb,'Position',[0.85 0.6 0.04 0.3])
export_fig([path_fig filesep 'LocalSleep_compERPs_FaceTopo.eps'],'-r 300')

figure;
temp_topo=squeeze(mean(mean((motorERP.Corr_erp(:,:,xTime3==-0.05)),3),1));
simpleTopoPlot_ft(temp_topo', motorERP.layout,'on',[],0,1);
title('motor')
format_fig;
hb=colorbar;
caxis([-1 1]*3)
set(hb,'Position',[0.85 0.6 0.04 0.3])
export_fig([path_fig filesep 'LocalSleep_compERPs_MotorTopo.eps'],'-r 300')
