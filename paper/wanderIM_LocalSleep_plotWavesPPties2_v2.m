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
Model_State_res=[];
Model_State_res2=[];
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
%     for nTask=1:2
        nWav_mdl1= fitlme(table_all(table_all.Elec==nE,:),sprintf('nWave~1+Task+State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        nP2P_mdl1= fitlme(table_all(table_all.Elec==nE,:),sprintf('P2P~1+Task+State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        nDSl_mdl1= fitlme(table_all(table_all.Elec==nE,:),sprintf('DownSlope~1+Task+State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        nUSl_mdl1= fitlme(table_all(table_all.Elec==nE,:),sprintf('UpSlope~1+Task+State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        
        nWav_mdl2= fitlme(table_all(table_all.Elec==nE,:),sprintf('MWvsON~1+Task+nWave+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        nP2P_mdl2= fitlme(table_all(table_all.Elec==nE,:),sprintf('MWvsON~1+Task+P2P+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        nDSl_mdl2= fitlme(table_all(table_all.Elec==nE,:),sprintf('MWvsON~1+Task+DownSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        nUSl_mdl2= fitlme(table_all(table_all.Elec==nE,:),sprintf('MWvsON~1+Task+UpSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        
        nWav_mdl3= fitlme(table_all(table_all.Elec==nE,:),sprintf('MBvsON~1+Task+nWave+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        nP2P_mdl3= fitlme(table_all(table_all.Elec==nE,:),sprintf('MBvsON~1+Task+P2P+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        nDSl_mdl3= fitlme(table_all(table_all.Elec==nE,:),sprintf('MBvsON~1+Task+DownSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        nUSl_mdl3= fitlme(table_all(table_all.Elec==nE,:),sprintf('MBvsON~1+Task+UpSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        
        nWav_mdl4= fitlme(table_all(table_all.Elec==nE,:),sprintf('MWvsMB~1+Task+nWave+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        nP2P_mdl4= fitlme(table_all(table_all.Elec==nE,:),sprintf('MWvsMB~1+Task+P2P+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        nDSl_mdl4= fitlme(table_all(table_all.Elec==nE,:),sprintf('MWvsMB~1+Task+DownSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        nUSl_mdl4= fitlme(table_all(table_all.Elec==nE,:),sprintf('MWvsMB~1+Task+UpSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
        
        Model_State_res{1,1}(nE,:)=[double(nWav_mdl1.Coefficients(3,4)) double(nP2P_mdl1.Coefficients(3,4)) double(nDSl_mdl1.Coefficients(3,4)) double(nUSl_mdl1.Coefficients(3,4)) ...
            double(nWav_mdl1.Coefficients(4,4)) double(nP2P_mdl1.Coefficients(4,4)) double(nDSl_mdl1.Coefficients(4,4)) double(nUSl_mdl1.Coefficients(4,4))];
        Model_State_res{2,1}(nE,:)=[double(nWav_mdl1.Coefficients(3,6)) double(nP2P_mdl1.Coefficients(3,6)) double(nDSl_mdl1.Coefficients(3,6)) double(nUSl_mdl1.Coefficients(3,6)) ...
            double(nWav_mdl1.Coefficients(4,6)) double(nP2P_mdl1.Coefficients(4,6)) double(nDSl_mdl1.Coefficients(4,6)) double(nUSl_mdl1.Coefficients(4,6))];
        
        Model_State_res2{1,1}(nE,:)=[double(nWav_mdl2.Coefficients(3,4)) double(nP2P_mdl2.Coefficients(3,4)) double(nDSl_mdl2.Coefficients(3,4)) double(nUSl_mdl2.Coefficients(3,4)) ...
            double(nWav_mdl3.Coefficients(3,4)) double(nP2P_mdl3.Coefficients(3,4)) double(nDSl_mdl3.Coefficients(3,4)) double(nUSl_mdl3.Coefficients(3,4)) ...
            double(nWav_mdl4.Coefficients(3,4)) double(nP2P_mdl4.Coefficients(3,4)) double(nDSl_mdl4.Coefficients(3,4)) double(nUSl_mdl4.Coefficients(3,4))];
        Model_State_res2{2,1}(nE,:)=[double(nWav_mdl2.Coefficients(3,6)) double(nP2P_mdl2.Coefficients(3,6)) double(nDSl_mdl2.Coefficients(3,6)) double(nUSl_mdl2.Coefficients(3,6)) ...
            double(nWav_mdl3.Coefficients(3,6)) double(nP2P_mdl3.Coefficients(3,6)) double(nDSl_mdl3.Coefficients(3,6)) double(nUSl_mdl3.Coefficients(3,6)) ...
            double(nWav_mdl4.Coefficients(3,6)) double(nP2P_mdl4.Coefficients(3,6)) double(nDSl_mdl4.Coefficients(3,6)) double(nUSl_mdl4.Coefficients(3,6))];
%     end
end

%%
% for nTask=1:2
Titles={'nWave','P2P','DownSlope','UpSlope'};
AxesLim={[],[-5 5],[-5 5],[-5 5]};
for nP=[3]
    figure;
%     subplot(1,3,1)
    temp_topo=Model_State_res2{1,1}(:,nP);
    temp_pV=Model_State_res2{2,1}(:,nP);
    
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    colorbar; caxis(AxesLim{nP}); %caxis([-1 1]*max(abs(temp_topo)))
    
    load(path_PsychFTlayout);
    if ~isempty(find(temp_pV<fdr(temp_pV,0.05)))
        ft_plot_lay_me(layout, 'chanindx',find(temp_pV<fdr(temp_pV,0.05)),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','no')
    end
    title([Titles{nP} ' - FDR: ' num2str(fdr(temp_pV,0.05))])
    export_fig([path_fig filesep 'LocalSleep_MWvsON_' Titles{nP} '.fig'])
    export_fig([path_fig filesep 'LocalSleep_MWvsON_' Titles{nP} '.eps'],'-r 300')

    figure;
    temp_topo=Model_State_res2{1,1}(:,nP+4);
    temp_pV=Model_State_res2{2,1}(:,nP+4);
    
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    colorbar; caxis(AxesLim{nP}); % caxis([-1 1]*max(abs(temp_topo)))
    
    load(path_PsychFTlayout);
    if ~isempty(find(temp_pV<fdr(temp_pV,0.05)))
        ft_plot_lay_me(layout, 'chanindx',find(temp_pV<fdr(temp_pV,0.05)),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','no')
    end
    title([Titles{nP} ' - FDR: ' num2str(fdr(temp_pV,0.05))])
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
%%
figure; set(gcf,'Position',[ 562   669   325*3   316]);
for nState=1:3
    subplot(2,3,nState)
    temp_topo=grpstats(table_all.nWave(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
    temp_topo([10 21])=NaN;
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    colorbar; caxis([0 2.5])
    % title('Down Slope \muV.s^-^1')
    format_fig;
    
    if nState>1
        subplot(2,3,nState+3)
        temp_topo=grpstats(table_all.nWave(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState)))-...
            grpstats(table_all.nWave(table_all.State==num2str(1)),table_all.Elec(table_all.State==num2str(1))); %in seconds but for all probes (60) so equivalent of in minutes
        temp_topo([10 21])=NaN;
        simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
        colorbar; caxis([-1 1]*0.4)
        % title('Down Slope \muV.s^-^1')
        format_fig;
    end
end

figure; set(gcf,'Position',[ 562   669   325*3   316]);
for nState=1:3
    subplot(2,3,nState)
    temp_topo=grpstats(table_all.P2P(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
    temp_topo([10 21])=NaN;
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    colorbar; caxis([0 60])
    % title('Down Slope \muV.s^-^1')
    format_fig;
    
    if nState>1
        subplot(2,3,nState+3)
        temp_topo=grpstats(table_all.P2P(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState)))-...
            grpstats(table_all.P2P(table_all.State==num2str(1)),table_all.Elec(table_all.State==num2str(1))); %in seconds but for all probes (60) so equivalent of in minutes
        temp_topo([10 21])=NaN;
        simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
        colorbar; caxis([-1 1]*12)
        % title('Down Slope \muV.s^-^1')
        format_fig;
    end
end

figure; set(gcf,'Position',[ 562   669   325*3   316]);
for nState=1:3
    subplot(2,3,nState)
    temp_topo=grpstats(table_all.DownSlope(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
    temp_topo([10 21])=NaN;
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    colorbar; caxis([200 550])
    % title('Down Slope \muV.s^-^1')
    format_fig;
    
    if nState>1
        subplot(2,3,nState+3)
        temp_topo=grpstats(table_all.DownSlope(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState)))-...
            grpstats(table_all.DownSlope(table_all.State==num2str(1)),table_all.Elec(table_all.State==num2str(1))); %in seconds but for all probes (60) so equivalent of in minutes
        temp_topo([10 21])=NaN;
        simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
        colorbar; caxis([-1 1]*250)
        % title('Down Slope \muV.s^-^1')
        format_fig;
    end
end
%%
figure; set(gcf,'Position',[ 562   669   325   316]);
% subplot(1,3,1)
temp_topo=grpstats(table_all.nWave,table_all.Elec); %in seconds but for all probes (60) so equivalent of in minutes
temp_topo([10 21])=NaN;
simpleTopoPlot_ft(temp_topo/20*60, path_PsychFTlayout,'on',[],0,1);
colorbar; %caxis([0 9])
title('SW density /min')
format_fig;
export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo.fig'])
export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo.eps'],'-r 300')

figure; set(gcf,'Position',[ 562   669   325   316]);
% subplot(1,3,1)
temp_topo=grpstats(table_all.P2P,table_all.Elec); %in seconds but for all probes (60) so equivalent of in minutes
temp_topo([10 21])=NaN;
simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
colorbar; %caxis([0 9])
title('P2P amp \muV')
format_fig;
export_fig([path_fig filesep 'LocalSleep_SWppties_P2PTopo.fig'])
export_fig([path_fig filesep 'LocalSleep_SWppties_P2PTopo.eps'],'-r 300')

figure; set(gcf,'Position',[ 562   669   325   316]);
% subplot(1,3,1)
temp_topo=grpstats(table_all.DownSlope,table_all.Elec); %in seconds but for all probes (60) so equivalent of in minutes
temp_topo([10 21])=NaN;
simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
colorbar; %caxis([0 9])
title('Down Slope \muV.s^-^1')
format_fig;
export_fig([path_fig filesep 'LocalSleep_SWppties_DownSlopeTopo.fig'])
export_fig([path_fig filesep 'LocalSleep_SWppties_DownSlopeTopo.eps'],'-r 300')

figure; set(gcf,'Position',[ 562   669   325   316]);
% subplot(1,3,1)
temp_topo=grpstats(table_all.UpSlope,table_all.Elec); %in seconds but for all probes (60) so equivalent of in minutes
temp_topo([10 21])=NaN;
simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
colorbar; %caxis([0 9])
title('Up Slope \muV.s^-^1')
format_fig;
export_fig([path_fig filesep 'LocalSleep_SWppties_UpSlopeTopo.fig'])
export_fig([path_fig filesep 'LocalSleep_SWppties_UpSlopeTopo.eps'],'-r 300')

%%
fprintf('E:%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    nWav_mdl1= fitlme(table_all(table_all.Elec==nE,:),sprintf('Vig~1+nWave+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
    nP2P_mdl1= fitlme(table_all(table_all.Elec==nE,:),sprintf('Vig~1+P2P+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
    nDSl_mdl1= fitlme(table_all(table_all.Elec==nE,:),sprintf('Vig~1+DownSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
    nUSl_mdl1= fitlme(table_all(table_all.Elec==nE,:),sprintf('Vig~1+UpSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
    
    Model_res{1}(nE,:)=[double(nWav_mdl1.Coefficients(2,4)) double(nP2P_mdl1.Coefficients(2,4)) double(nDSl_mdl1.Coefficients(2,4)) double(nUSl_mdl1.Coefficients(2,4))];
    Model_res{2}(nE,:)=[double(nWav_mdl1.Coefficients(2,6)) double(nP2P_mdl1.Coefficients(2,6)) double(nDSl_mdl1.Coefficients(2,6)) double(nUSl_mdl1.Coefficients(2,6))];
end

%%
% Titles={'mWave','P2P','DownSlope','UpSlope'};
% for nP=1:length(Titles)
%     figure;
%     temp_topo=Model_res{1}(:,nP);
%     temp_pV=Model_res{2}(:,nP);
%     
%     simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
%     colorbar; caxis([-1 1]*max(abs(temp_topo)))
%     
%     load(path_PsychFTlayout);
%     if ~isempty(find(temp_pV<fdr(temp_pV,0.05)))
%         ft_plot_lay_me(layout, 'chanindx',find(temp_pV<fdr(temp_pV,0.05)),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','no')
%     end
%     title([Titles{nP} ' - FDR: ' num2str(fdr(temp_pV,0.05))])
% %     export_fig([path_fig filesep 'LocalSleep_CorrVig_' Titles{nP} '.fig'])
% %     export_fig([path_fig filesep 'LocalSleep_CorrVig_' Titles{nP} '.eps'],'-r 300')
% end
%%
uniqueSubs=unique(Ppties_Waves(:,1));
Ppties_Waves_av=    [];
for nS=1:length(uniqueSubs)
    for nPr=1:60
    Ppties_Waves_av=[Ppties_Waves_av ; nanmean(Ppties_Waves(Ppties_Waves(:,1)==uniqueSubs(nS) & Ppties_Waves(:,2)==nPr,:))];
    end
end

%%
table_av=array2table(Ppties_Waves_av,'VariableNames',{'SubID','nPr','Elec','ProbeN','BlockN','Task','State','Vig','nWave','Dur','Freq','P2P','xNeg','xPos','endPos','xNegPeak','NegPeak','xPosPeak','PosPeak','DownSlope','UpSlope','MaxAmpW','MinAmpW'});
table_av.State(table_av.State==4)=3;
% table_av.Vig=abs(table_av.Vig-5);
table_av.SubID=categorical(table_av.SubID);
table_av.Task=categorical(table_av.Task);
table_av.State=categorical(table_av.State);

%%
myS=unique(table_av.SubID);
all_zval=[]; all_val=[]; all_val2=[];
for nS=1:length(myS)
    for nTask=1:2
        all_zval=[all_zval ; [zscore(table_av.Vig(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            zscore(table_av.nWave(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            zscore(table_av.P2P(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            zscore(table_av.DownSlope(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            zscore(table_av.UpSlope(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            double(table_av.State(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask)))]];
        
        temp=[(table_av.Vig(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            (table_av.nWave(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            (table_av.P2P(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            (table_av.DownSlope(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            (table_av.UpSlope(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            double(table_av.SubID(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            double(table_av.State(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask)))];
        
%         temp(:,2:end-2)=temp(:,2:end-2)-repmat(mean(temp(:,2:end-2),1),size(temp,1),1);
all_val=[all_val ; temp];
        
%         all_val2=[all_val2 ; [(table_all.Vig(table_all.SubID==myS(nS) & table_all.Task==num2str(nTask))) ...
%             (table_all.nWave(table_all.SubID==myS(nS) & table_all.Task==num2str(nTask))) ...
%             (table_all.P2P(table_all.SubID==myS(nS) & table_all.Task==num2str(nTask))) ...
%             (table_all.DownSlope(table_all.SubID==myS(nS) & table_all.Task==num2str(nTask))) ...
%             (table_all.UpSlope(table_all.SubID==myS(nS) & table_all.Task==num2str(nTask))) ...
%             double(table_all.SubID(table_all.SubID==myS(nS) & table_all.Task==num2str(nTask))) ...
%             double(table_all.State(table_all.SubID==myS(nS) & table_all.Task==num2str(nTask)))]];
    end
end
 
%%
Titles={'nWave','P2P','DownSlope','UpSlope'};
for nT=1:length(Titles)
    figure; set(gcf,'Position',[ 562   669   325   316]);
    simpleCorPlotbin(-1*all_zval(:,1),all_zval(:,nT+1),{'o','k',[1 1 1]*0.5,288,10,4},'Pearson',0)
    format_fig;
    xlabel('Vig'); ylabel(Titles{nT});
    export_fig([path_fig filesep 'LocalSleep_SWppties_CorrZvalAllE_' Titles{nT} '.fig'])
export_fig([path_fig filesep 'LocalSleep_SWppties_CorrZvalAllE_' Titles{nT} '.eps'],'-r 300')
end

%%
mdl1_State_PPwave=[];
mdl1_State_PPwave{1,4}=fitlme(table_av,'Vig~1+Task+(1|SubID)');
mdl1_State_PPwave{2,4}=fitlme(table_av,'nWave~1+Task+(1|SubID)');
mdl1_State_PPwave{3,4}=fitlme(table_av,'P2P~1+Task+(1|SubID)');
mdl1_State_PPwave{4,4}=fitlme(table_av,'DownSlope~1+Task+(1|SubID)');

mdl1_State_PPwave{1,3}=fitlme(table_av,'Vig~1+Task+State+(1|SubID)');
mdl1_State_PPwave{2,3}=fitlme(table_av,'nWave~1+Task+State+(1|SubID)');
mdl1_State_PPwave{3,3}=fitlme(table_av,'P2P~1+Task+State+(1|SubID)');
mdl1_State_PPwave{4,3}=fitlme(table_av,'DownSlope~1+Task+State+(1|SubID)');

for nP=1:4
   ComModel=compare(mdl1_State_PPwave{nP,4},mdl1_State_PPwave{nP,3}); 
   CompModel_pval(nP)=(ComModel.pValue(2));
end
for nTask=1:2
%     mdl0_State_nWave2{nTask}=fitlme(table_av,'nWave~1+Task+(1|SubID)');
%     mdl0_State_P2P2{nTask}=fitlme(table_av,'P2P~1+Task+(1|SubID)');
%     mdl0_State_DownSlope2{nTask}=fitlme(table_av,'DownSlope~1+Task+(1|SubID)');
    
    mdl1_State_PPwave{1,nTask}=fitlme(table_av(table_av.Task==num2str(nTask),:),'Vig~1+State+(1|SubID)');
    mdl1_State_PPwave{2,nTask}=fitlme(table_av(table_av.Task==num2str(nTask),:),'nWave~1+State+(1|SubID)');
    mdl1_State_PPwave{3,nTask}=fitlme(table_av(table_av.Task==num2str(nTask),:),'P2P~1+State+(1|SubID)');
    mdl1_State_PPwave{4,nTask}=fitlme(table_av(table_av.Task==num2str(nTask),:),'DownSlope~1+State+(1|SubID)');
    
end
% mdl2_State_nWave=fitlme(table_av,'nWave~1+Task*State+(1|SubID)');
% mdl2_State_P2P=fitlme(table_av,'P2P~1+Task*State+(1|SubID)');
% mdl2_State_DownSlope=fitlme(table_av,'DownSlope~1+Task*State+(1|SubID)');


%%
Titles={'Vig','nWave','P2P','DownSlope'};
figure; set(gcf,'Position',[680   246   492   434])
MarkerT={'d','o'};
YLims=[[-0.4 1]*1; [-0.4 1]*0.2; [-0.4 1]*1.2; [-0.4 1]*20];
for nP=1:4
    subplot(2,2,nP); hold on;
    for nState=2:3
        line((nState-1).*[1 1],double(mdl1_State_PPwave{nP,3}.Coefficients(nState+1,2))+double(mdl1_State_PPwave{nP,3}.Coefficients(nState+1,3)).*[-1 1],'Color',Colors(nState,:),'LineWidth',3);
        if double(mdl1_State_PPwave{nP,3}.Coefficients(nState+1,6))<0.05
            MarkerCol=Colors(nState,:);
        else
            MarkerCol='w';
        end
        fprintf('... ... %s state %g pV=%1.6f\n',Titles{nP},nState,double(mdl1_State_PPwave{nP,3}.Coefficients(nState+1,6)))
        scatter((nState-1),double(mdl1_State_PPwave{nP,3}.Coefficients(nState+1,2)),'Marker','o','SizeData',200,'LineWidth',3,'MarkerFaceColor',MarkerCol,'MarkerEdgeColor',Colors(nState,:),'MarkerFaceAlpha',0.7)
    end
    format_fig;
    xlim([0.3 2.7])
    line(xlim,[0 0],'Color',[1 1 1]*0.7,'LineStyle','--','LineWidth',2)
    set(gca,'XTick',1:2,'XTickLabel',{'MW','MB'})
    ylim(YLims(nP,:))
    ylabel(Titles{nP})
end
export_fig([path_fig filesep 'LocalSleep_SWppties_LMEbyProbes_ModelRes.fig'])
export_fig([path_fig filesep 'LocalSleep_SWppties_LMEbyProbes_ModelRes.eps'],'-r 300')
%%
figure; set(gcf,'Position',[680   380   253   598])
MarkerT={'d','o'};
for nP=1:3
subplot(3,1,nP); hold on;
for nTask=1:2
for nState=2:3
line((nState-1).*[1 1]+(2*nTask-3)*0.2,double(mdl1_State_PPwave{nP,nTask}.Coefficients(nState,2))+double(mdl1_State_PPwave{nP,nTask}.Coefficients(nState,3)).*[-1 1],'Color',Colors(nState,:),'LineWidth',3);
if double(mdl1_State_PPwave{nP,nTask}.Coefficients(nState,6))<0.05
    MarkerCol=Colors(nState,:);
else
    MarkerCol='w';
end
scatter((nState-1)+(2*nTask-3)*0.2,double(mdl1_State_PPwave{nP,nTask}.Coefficients(nState,2)),'Marker',MarkerT{nTask},'SizeData',200,'LineWidth',3,'MarkerFaceColor',MarkerCol,'MarkerEdgeColor',Colors(nState,:))
end
end
format_fig;
xlim([0.3 2.7])
line(xlim,[0 0],'Color',[1 1 1]*0.7,'LineStyle','--','LineWidth',2)
set(gca,'XTick',1:2,'XTickLabel',{'MW','MB'})
end

%%
Titles={'nWave','P2P','DownSlope'};
TitlesCol=[2 3 4];
figure; set(gcf,'Position',[680    10   407   795]);
YLims={[0 6],[20 80],[100 1300]};
barplot=0;
violinplot=1;
boxxplot=0;
for nT=1:length(Titles)
    subplot(3,1,nT); format_fig; hold on;
    for nState=1:3
        if barplot
            simpleBarPlot(nState,all_val(all_zval(:,size(all_zval,2))==nState,TitlesCol(nT)),Colors(nState,:),0.9,'k',[],3);
            ylim(YLims{nT})
        elseif violinplot
            temp=all_val(all_val(:,size(all_val,2))==nState,TitlesCol(nT));
            simple_violinplot(temp,[],nState,'ViolinColor',Colors(nState,:),'Width',.38,'BoxWidth',0.1,'ShowMean',true);
        elseif boxxplot
            temp=all_val(all_val(:,size(all_val,2))==nState,TitlesCol(nT));
            boxplot(temp,'positions',nState);
        else
        h1{nState} = raincloud_plot(all_val(all_val(:,size(all_val,2))==nState,TitlesCol(nT)), 'box_on', 1, 'color', Colors(nState,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .15+(nState-1)*0.2, 'dot_dodge_amount', .15+(nState-1)*0.2,...
            'box_col_match', 0); %h1{nState}.SizeData=72; h1{nState}.MarkerFaceAlpha=0.5;
%         set(gca,'Xlim',YLims{nT});
        end
        %         h2 = raincloud_plot(data{2,ntask}, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
        %             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
        %         h2{2}.SizeData=72; h2{2}.MarkerFaceAlpha=0.5;
        %         h3 = raincloud_plot(data{3,ntask}, 'box_on', 1, 'color', Colors(3,:), 'alpha', 0.5,...
        %             'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
        %         h3{2}.SizeData=72; h3{2}.MarkerFaceAlpha=0.5;
    end
end
% %%
% data=[];
% for i = 1:4
%     for j = 1:2
%         temp=table_av.nWave(table_av.Task==num2str(j) & table_av.Vig==i);
%         tempS=table_av.SubID(table_av.Task==num2str(j) & table_av.Vig==i);
%         tempbyS=[]; myS=unique(tempS);
%         for nS=1:length(myS)
%             tempbyS(nS)=nanmean(temp(tempS==myS(nS)));
%         end
%         data{i, j} = tempbyS;
%     end
% end
% 
% % make figure
% figure;
% set(gcf,'position',[11   513   375   445])
% h   = rm_raincloud(data, Colors(1:2,:));
% set(gca, 'XLim', [1 4],'YTick','');
% 
% % change one subset to new colour and alter dot size
% for i=1:4
%     % scatter
%     h.s{i, 2}.MarkerFaceColor   = Colors(i,:);
%     h.s{i, 2}.MarkerEdgeColor   = [1 1 1]*0.5;
%     h.s{i, 2}.SizeData          = 100;
%     h.s{i, 2}.LineWidth          = 2;
%     h.s{i,2}.YData=h.s{i,2}.YData-0.1;
%     
%     % patch
%     h.p{i, 2}.FaceColor     = Colors(i,:);
%     h.p{i, 2}.LineWidth          = 2;
%     
%     % scatter mean
%   h.m(i,2).MarkerFaceAlpha=0.7;
%             h.m(i,2).MarkerFaceColor=[1 1 1]*0.5;
%             h.m(i,2).MarkerEdgeColor=[0 0 0];%[190 90 255]./255;
%             h.m(i,2).SizeData=300;
%             h.m(i,2).LineWidth=4;
%     % line mean
%     if i<3
%         h.l(i,2).Color=[1 1 1]*0;
%     end
% end
% 
% for i=1:4
%     h.s{i, 1}.Marker         = 'd';
%     h.s{i, 1}.MarkerEdgeColor   = [1 1 1]*0.5;
%     h.s{i, 1}.MarkerFaceColor   = Colors(i,:);
%     h.s{i, 1}.SizeData          = 100;
%     h.s{i, 1}.LineWidth          = 2;
%     h.s{i,1}.YData=h.s{i,1}.YData+0.1;
%     
%     h.m(i, 1).MarkerEdgeColor   = 'none';
%     h.m(i, 1).MarkerFaceColor   = Colors(i,:);
%     
%     h.p{i, 1}.FaceAlpha   = 0;
%     h.p{i, 1}.EdgeColor= Colors(i,:);
%     h.p{i, 1}.LineWidth          = 1;
%     
%     
%     h.m(i,1).Marker='d';
%     h.m(i,1).MarkerFaceAlpha=0.7;
%     h.m(i,1).MarkerFaceColor=[1 1 1]*0.5;
%     h.m(i,1).MarkerEdgeColor=[0 0 0];
%     h.m(i,1).SizeData=300;
%     h.m(i,1).LineWidth=4;
%     if i<3
%         h.l(i,1).Color=[1 1 1]*0;
%         h.l(i,1).LineStyle='-';
%     end
% end
% format_fig;
% % export_fig([path_fig filesep 'WanderIM_behav_fatigue_byStateAndTask.eps'],'-r 300')