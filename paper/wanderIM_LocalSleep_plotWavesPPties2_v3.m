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
Ppties_Waves=[];
nSc=0;

for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    SubID=filename;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    if ~ismember(SubID,GoodSubID)
        continue;
    end
    
    
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
    if AmpCriterionIdx==9
        all_Waves(:,AmpCriterionIdx)=-all_Waves(:,AmpCriterionIdx);
    end
    %     hdr=ft_read_header([eeg_path filesep filename]);
    
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' SubID '_*.mat']);
    load([behav_path filesep behav_file.name]);
    all_Waves=double(all_Waves);
    all_Waves(all_Waves(:,5)<-20*500 | all_Waves(:,7)>0,:)=[];
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./500);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,AmpCriterionIdx)>art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>max_posampl | all_Waves(:,14)>art_ampl| abs(all_Waves(:,15))>art_ampl)*100)
    all_Waves(all_freq>max_Freq | all_Waves(:,AmpCriterionIdx)>art_ampl | all_Waves(:,11)>max_posampl| all_Waves(:,14)>art_ampl| abs(all_Waves(:,15))>art_ampl,:)=[];
    
    
    for nE=1:63
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./500;
        temp_freq=1./temp_len;
        temp_abs=1./temp_len;
        temp_p2p=thisE_Waves(:,AmpCriterionIdx);
        temp_p2p2=thisE_Waves(:,4);
        
        
        if ~isempty(fixThr)
            thr_Wave(nSc,nE)=fixThr;
        else
            thr_Wave(nSc,nE)=prctile(thisE_Waves(:,AmpCriterionIdx),prticle_Thr);
        end
        for nPr=1:60
            sel_Waves=thisE_Waves(thisE_Waves(:,2)==nPr & temp_p2p>thr_Wave(nSc,nE),:);
            if size(sel_Waves,1)~=0
            temp_len2=abs((sel_Waves(:,5)-sel_Waves(:,7)))./500;
        temp_freq2=1./temp_len2;
            dens_Waves(nSc,nE)=size(sel_Waves,1);
            Ppties_Waves=[Ppties_Waves ; [n nPr nE probe_res(nPr,1) probe_res(nPr,4) probe_res(nPr,5) probe_res(nPr,32) probe_res(nPr,38) dens_Waves(nSc,nE) mean(temp_len2) mean(temp_freq2) mean(sel_Waves(:,[4:15]),1)]];            
            else
                Ppties_Waves=[Ppties_Waves ; [n nPr nE probe_res(nPr,1) probe_res(nPr,4) probe_res(nPr,5) probe_res(nPr,32) probe_res(nPr,38) 0 nan(1,14)]];
            end
            end
    end
   
end

%%
addpath((path_fieldtrip)); ft_defaults;

%
table_all=array2table(Ppties_Waves,'VariableNames',{'SubID','nPr','Elec','ProbeN','BlockN','Task','State','Vig','nWave','Dur','Freq','P2P','xNeg','xPos','endPos','xNegPeak','NegPeak','xPosPeak','PosPeak','DownSlope','UpSlope','MaxAmpW','MinAmpW'});
table_all.State(table_all.State==4)=3;
% table_av.Vig=abs(table_av.Vig-5);
table_all.SubID=categorical(table_all.SubID);
table_all.Task=categorical(table_all.Task);
table_all.State=categorical(table_all.State);

table_all.MWvsON=nan(size(table_all,1),1);
table_all.MWvsON(table_all.State=='2')=1;
table_all.MWvsON(table_all.State=='1')=0;
table_all.MBvsON=nan(size(table_all,1),1);
table_all.MBvsON(table_all.State=='3')=1;
table_all.MBvsON(table_all.State=='1')=0;

table_all.MWvsMB=nan(size(table_all,1),1);
table_all.MWvsMB(table_all.State=='2')=1;
table_all.MWvsMB(table_all.State=='3')=0;

fprintf('E:%2.0f/63\n',0)
totperm=1000;
MWvsON_P2P_est=cell(1,2);
MWvsON_DwS_est=cell(1,2);
MBvsON_P2P_est=cell(1,2);
MBvsON_DwS_est=cell(1,2);
for nE=1:63
    tic;
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    [real_out, perm_out]=lme_perm(table_all(table_all.Elec==nE,:),'P2P','MWvsON~1+pred+Task+(1|SubID)',totperm);
    MWvsON_P2P_est{1}=[MWvsON_P2P_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MWvsON_P2P_est{2}=[MWvsON_P2P_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    
    [real_out, perm_out]=lme_perm(table_all(table_all.Elec==nE,:),'DownSlope','MWvsON~1+pred+Task+(1|SubID)',totperm);
    MWvsON_DwS_est{1}=[MWvsON_DwS_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MWvsON_DwS_est{2}=[MWvsON_DwS_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    
    [real_out, perm_out]=lme_perm(table_all(table_all.Elec==nE,:),'P2P','MBvsON~1+pred+Task+(1|SubID)',totperm);
    MBvsON_P2P_est{1}=[MBvsON_P2P_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MBvsON_P2P_est{2}=[MBvsON_P2P_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    
    [real_out, perm_out]=lme_perm(table_all(table_all.Elec==nE,:),'DownSlope','MBvsON~1+pred+Task+(1|SubID)',totperm);
    MBvsON_DwS_est{1}=[MBvsON_DwS_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MBvsON_DwS_est{2}=[MBvsON_DwS_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
toc;
end
save('/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm.mat','MWvsON_P2P_est','MWvsON_DwS_est','MBvsON_P2P_est','MBvsON_DwS_est')

%% Compute monte cartlo p
clus_alpha=0.01;
montecarlo_alpha=0.05/8;

addpath((path_fieldtrip)); ft_defaults;
cfg_neighb=[];
cfg_neighb.method = 'tri';
cfg_neighb.layout=path_PsychFTlayout;
neighbours = ft_prepare_neighbours(cfg_neighb);

[MWvsON_P2P_clus]=get_clusterperm_lme(MWvsON_P2P_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[MWvsON_DwS_clus]=get_clusterperm_lme(MWvsON_DwS_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[MBvsON_P2P_clus]=get_clusterperm_lme(MBvsON_P2P_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[MBvsON_DwS_clus]=get_clusterperm_lme(MBvsON_DwS_est,clus_alpha,montecarlo_alpha,totperm,neighbours);

%%
% for nTask=1:2
Titles={'P2P','DownSlope'};
AxesLim={[-5 5],[-5 5]};
for nP=1:2
    figure;
    if nP==1
        temp_topo=MWvsON_P2P_est{1}(:,3);
        temp_clus=MWvsON_P2P_clus;
    else
        temp_topo=MWvsON_DwS_est{1}(:,3);
        temp_clus=MWvsON_DwS_clus;
    end

    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    colorbar; caxis(AxesLim{nP}); %caxis([-1 1]*max(abs(temp_topo)))
    
    load(path_PsychFTlayout);
   if ~isempty(temp_clus)
        for nclus=1:length(temp_clus)
        ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','no')
        fprintf('... ... found %s cluster (%g) of %g electrodes (tval cluster=%g, Pmc=%g)\n',temp_clus{nclus}{1},nclus,length(temp_clus{nclus}{2}),temp_clus{nclus}{3},temp_clus{nclus}{4})
        end
    end
    title([Titles{nP}])
    export_fig([path_fig filesep 'LocalSleep_MWvsON_' Titles{nP} '.fig'])
    export_fig([path_fig filesep 'LocalSleep_MWvsON_' Titles{nP} '.eps'],'-r 300')

    figure;
   if nP==1
        temp_topo=MBvsON_P2P_est{1}(:,3);
        temp_clus=MBvsON_P2P_clus;
    else
        temp_topo=MBvsON_DwS_est{1}(:,3);
        temp_clus=MBvsON_DwS_clus;
    end
    
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    colorbar; caxis(AxesLim{nP}); % caxis([-1 1]*max(abs(temp_topo)))
    
    load(path_PsychFTlayout);
    if ~isempty(temp_clus)
        for nclus=1:length(temp_clus)
        ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','no')
        fprintf('... ... found %s cluster (%g) of %g electrodes (tval cluster=%g, Pmc=%g)\n',temp_clus{nclus}{1},nclus,length(temp_clus{nclus}{2}),temp_clus{nclus}{3},temp_clus{nclus}{4})
        end
    end
    title([Titles{nP}])
    export_fig([path_fig filesep 'LocalSleep_MBvsON_' Titles{nP} '.fig'])
    export_fig([path_fig filesep 'LocalSleep_MBvsON_' Titles{nP} '.eps'],'-r 300')

%      subplot(1,3,3)
%     temp_topo=Model_State_res2{1,1}(:,nP+8);
%     temp_pV=Model_State_res2{2,1}(:,nP+8);
%     
%     simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
%     colorbar; caxis([-1 1]*max(abs(temp_topo)))
%     
%     load(path_PsychFTlayout);
%     if ~isempty(find(temp_pV<0.05))
%         ft_plot_lay_me(layout, 'chanindx',find(temp_pV<0.05),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','no')
%     end
%     title([Titles{nP} ' - FDR: ' num2str(fdr(temp_pV,0.05))])
end
% end
