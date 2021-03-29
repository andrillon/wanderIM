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

elecs_ERP=1:63;
for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    SubID=filename;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    if ~strcmp(SubID,'305')
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
%                 if max(abs(temp))>200
%                     temp_ERP_Waves(nW,:)=nan(1,length((-1*D.fsample:2*D.fsample)));
%                     continue;
%                 end
                temp_ERP_Waves(nW,:)=temp;
                temp_ERP_Waves(nW,:)=temp_ERP_Waves(nW,:)-mean(temp_ERP_Waves(nW,:));
            end
            numERP_Waves{nE}(nSc)=size(temp_ERP_Waves,1);
            if numERP_Waves{nE}(nSc)<50
                ERP_Waves{nE}(nSc,:)=nan;
            else
                ERP_Waves{nE}=(temp_ERP_Waves);
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
end


%%
visualERP2=load('/Users/tand0009/Data/WanderIM/ERPComp/wanderIM_visualERPs_305');
motorERP2=load('/Users/tand0009/Data/WanderIM/ERPComp/wanderIM_motorERPs_305');


% %%
% bins=-200:1:200;
% cmap=cbrewer('qual','Set1',3);
% figure;
% this_ch=match_str(visualERP2.layout.label,'Cz');
% [nout,~]=histc(visualERP2.DGT_erp_avwin(:,this_ch),bins);
% plot(bins,nout/max(nout),'Color',cmap(2,:),'LineWidth',2);
% 
% hold on;
% this_ch=match_str(visualERP2.layout.label,'Cz');
% [nout,~]=histc(visualERP2.FAC_erp_avwin(:,this_ch),bins);
% plot(bins,nout/max(nout),'Color',cmap(2,:),'LineWidth',2,'LineStyle','--');
% 
% hold on;
% this_ch=match_str(visualERP2.layout.label,'Cz');
% [nout,~]=histc(motorERP2.Corr_erp_avwin(:,this_ch),bins);
% plot(bins,nout/max(nout),'Color',cmap(3,:),'LineWidth',2);
% 
% hold on;
%     xTime1=-1:1/500:2;
% this_ch=match_str(layout.label,'Cz');
% temp=mean(ERP_Waves(:,xTime1>0 & xTime1<0.1),2);
% [nout,~]=histc(temp,bins);
% plot(bins,nout/max(nout),'Color',cmap(1,:),'LineWidth',2);


%%

figure; set(gcf,'Position',[231         473        1239         483]);
tight_subplot(1,4);
elecs_ERP_labels={'Fz','Cz','Pz','Oz'};
cLIM=50;
for nEl=2 %1:length(elecs_ERP)
    subplot(1,3,1);
    hold on;
    xTime1=-1:1/500:2;
    my_ERP=(ERP_Waves{match_str(layout.label,elecs_ERP_labels{nEl})});
    my_ERP(sum(isnan(my_ERP),2)>0,:)=[];
    Amp=mean(my_ERP(:,xTime1>0.0 & xTime1<0.1),2);
    [~,orderIdx]=sort(Amp);
    imagesc(xTime1,1:size(my_ERP,1),my_ERP);
    
    xlim([-.2 .6])
    ylim([1 120])
    xlabel('Time (s)')
%     ylabel('Trials')
    format_fig;
    line([0 0],ylim,'Color','k','LineStyle','--');
    caxis([-1 1]*cLIM)
    %
    subplot(1,3,2);
    hold on;
    xTime2=visualERP2.xTime;
    this_ch=match_str(visualERP2.layout.label,elecs_ERP_labels{nEl});
    my_ERP=[squeeze(visualERP2.STD_erp_D(:,this_ch,:)) ; squeeze(visualERP2.STD_erp_F(:,this_ch,:))];
    Amp=mean(my_ERP(:,xTime2>0.15 & xTime2<0.25),2);
    [~,orderIdx]=sort(Amp);
    imagesc(xTime2,1:size(my_ERP,1),my_ERP);
    
    xlim([-.2 .6])
    ylim([1 size(my_ERP,1)])
    xlabel('Time (s)')
%     ylabel('Trials')
    format_fig;
    line([0 0],ylim,'Color','k','LineStyle','--');
    caxis([-1 1]*cLIM)
    
%     subplot(1,4,3);
%     hold on;
%     xTime2=visualERP2.xTime;
%     this_ch=match_str(visualERP2.layout.label,elecs_ERP_labels{nEl});
%     my_ERP=[squeeze(visualERP2.STD_erp_D(:,this_ch,:))];
%     Amp=mean(my_ERP(:,xTime2>0.15 & xTime2<0.25),2);
%     [~,orderIdx]=sort(Amp);
%     imagesc(xTime2,1:size(my_ERP,1),my_ERP(orderIdx,:));
%     
%     xlim([-.2 .6])
%     ylim([1 size(my_ERP,1)])
%     xlabel('Time (s)')
% %     ylabel('Trials')
%     format_fig;
%     line([0 0],ylim,'Color','k','LineStyle','--');
%     caxis([-1 1]*75)
    ylim([1 120])

    
    subplot(1,3,3);
    hold on;
    xTime3=motorERP2.xTime;
    this_ch=match_str(motorERP2.layout.label,elecs_ERP_labels{nEl});
    my_ERP=squeeze(motorERP2.Corr_erp(this_ch,:,:))';
    Amp=mean(my_ERP(:,xTime3>-0.1 & xTime3<0.0),2);
    [~,orderIdx]=sort(Amp);
    imagesc(xTime3,1:size(my_ERP,1),my_ERP);
    
    
    xlim([-.4 .4])
    ylim([1 size(my_ERP,1)])
    xlabel('Time (s)')
%     ylabel('Trials')
    format_fig;
    line([0 0],ylim,'Color','k','LineStyle','--');
    caxis([-1 1]*cLIM)
    ylim([1 120])
end
colorbar('Position',[0.92 0.7 0.02 0.2])
export_fig([path_fig filesep 'LocalSleep_ERPcomp_sub305.png'],'-r 300')
