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
Ppties_Waves_bin=[];
nSc=0;

Boundaries=[-20 -15;-15 -10; -10 -5; -5 0];
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
                
                for nWin=1:size(Boundaries,1)
                    sel_Waves2=sel_Waves(sel_Waves(:,15)>=Boundaries(nWin,1) & sel_Waves(:,15)<Boundaries(nWin,2),:);
                    if size(sel_Waves2,1)~=0
                        temp_len2b=abs((sel_Waves2(:,5)-sel_Waves2(:,7)))./500;
                        temp_freq2b=1./temp_len2b;
                        dens_Waves2(nSc,nE)=size(sel_Waves2,1);
                        Ppties_Waves_bin=[Ppties_Waves_bin ; [n nPr nE probe_res(nPr,1) probe_res(nPr,4) probe_res(nPr,5) probe_res(nPr,32) probe_res(nPr,38) dens_Waves2(nSc,nE) mean(temp_len2b) mean(temp_freq2b) mean(sel_Waves2(:,[4:15]),1) nWin]];
                    end
                end
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
AxesLim={[-5 5],[-5 5],[-5 5],[-5 5]};
for nP=[1 2 3 4]
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
    %export_fig([path_fig filesep 'LocalSleep_MWvsON_' Titles{nP} '.fig'])
    %export_fig([path_fig filesep 'LocalSleep_MWvsON_' Titles{nP} '.eps'],'-r 300')
    
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
    %export_fig([path_fig filesep 'LocalSleep_MBvsON_' Titles{nP} '.fig'])
    %export_fig([path_fig filesep 'LocalSleep_MBvsON_' Titles{nP} '.eps'],'-r 300')
    
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
hb=colorbar; %caxis([0 9])
caxis([3 8])
set(hb,'YTick',[]);
title('SW density /min')
format_fig;
export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo.fig'])
export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo.eps'],'-r 300')

figure; set(gcf,'Position',[ 562   669   325   316]);
% subplot(1,3,1)
temp_topo=grpstats(table_all.P2P,table_all.Elec); %in seconds but for all probes (60) so equivalent of in minutes
temp_topo([10 21])=NaN;
simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
hb=colorbar; caxis([25 60])
set(hb,'YTick',[]);

title('P2P amp \muV')
format_fig;
export_fig([path_fig filesep 'LocalSleep_SWppties_P2PTopo.fig'])
export_fig([path_fig filesep 'LocalSleep_SWppties_P2PTopo.eps'],'-r 300')

figure; set(gcf,'Position',[ 562   669   325   316]);
% subplot(1,3,1)
temp_topo=grpstats(table_all.DownSlope,table_all.Elec); %in seconds but for all probes (60) so equivalent of in minutes
temp_topo([10 21])=NaN;
simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
hb=colorbar; caxis([250 480])
set(hb,'YTick',[]);
title('Down Slope \muV.s^-^1')
format_fig;
export_fig([path_fig filesep 'LocalSleep_SWppties_DownSlopeTopo.fig'])
export_fig([path_fig filesep 'LocalSleep_SWppties_DownSlopeTopo.eps'],'-r 300')


figure; set(gcf,'Position',[ 562   669   325   316]);
% subplot(1,3,1)
temp_topo=grpstats(table_all.UpSlope,table_all.Elec); %in seconds but for all probes (60) so equivalent of in minutes
temp_topo([10 21])=NaN;
simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
hb=colorbar; caxis([250 530])
set(hb,'YTick',[]);
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
% %     %export_fig([path_fig filesep 'LocalSleep_CorrVig_' Titles{nP} '.fig'])
% %     %export_fig([path_fig filesep 'LocalSleep_CorrVig_' Titles{nP} '.eps'],'-r 300')
% end
%%
filename = '/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_localsleep_amp_pup_thrE90P2P_Dec21_v5.txt';
% filename = '/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_localsleep_amp_pup_thrE90neg_Dec21_v5.txt';

delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

fclose(fileID);

res_table = table;
res_table.SubID = dataArray{:, 1};
res_table.BlockN = dataArray{:, 2};
res_table.TrialN = dataArray{:, 3};
res_table.ProbeN = dataArray{:, 4};
res_table.DistProbe = dataArray{:, 5};
res_table.Task = dataArray{:, 6};
res_table.StimCat = dataArray{:, 7};
res_table.Perf = dataArray{:, 8};
res_table.RCode = dataArray{:, 9};
res_table.RT = dataArray{:, 10};
res_table.W_Fp1 = dataArray{:, 11};
res_table.W_Fz = dataArray{:, 12};
res_table.W_F3 = dataArray{:, 13};
res_table.W_F7 = dataArray{:, 14};
res_table.W_FT9 = dataArray{:, 15};
res_table.W_FC5 = dataArray{:, 16};
res_table.W_FC1 = dataArray{:, 17};
res_table.W_C3 = dataArray{:, 18};
res_table.W_T7 = dataArray{:, 19};
res_table.W_TP9 = dataArray{:, 20};
res_table.W_CP5 = dataArray{:, 21};
res_table.W_CP1 = dataArray{:, 22};
res_table.W_Pz = dataArray{:, 23};
res_table.W_P3 = dataArray{:, 24};
res_table.W_P7 = dataArray{:, 25};
res_table.W_O1 = dataArray{:, 26};
res_table.W_Oz = dataArray{:, 27};
res_table.W_O2 = dataArray{:, 28};
res_table.W_P4 = dataArray{:, 29};
res_table.W_P8 = dataArray{:, 30};
res_table.W_TP10 = dataArray{:, 31};
res_table.W_CP6 = dataArray{:, 32};
res_table.W_CP2 = dataArray{:, 33};
res_table.W_Cz = dataArray{:, 34};
res_table.W_C4 = dataArray{:, 35};
res_table.W_T8 = dataArray{:, 36};
res_table.W_FT10 = dataArray{:, 37};
res_table.W_FC6 = dataArray{:, 38};
res_table.W_FC2 = dataArray{:, 39};
res_table.W_F4 = dataArray{:, 40};
res_table.W_F8 = dataArray{:, 41};
res_table.W_Fp2 = dataArray{:, 42};
res_table.W_AF7 = dataArray{:, 43};
res_table.W_AF3 = dataArray{:, 44};
res_table.W_F1 = dataArray{:, 45};
res_table.W_F5 = dataArray{:, 46};
res_table.W_FT7 = dataArray{:, 47};
res_table.W_FC3 = dataArray{:, 48};
res_table.W_C1 = dataArray{:, 49};
res_table.W_C5 = dataArray{:, 50};
res_table.W_TP7 = dataArray{:, 51};
res_table.W_CP3 = dataArray{:, 52};
res_table.W_P1 = dataArray{:, 53};
res_table.W_P5 = dataArray{:, 54};
res_table.W_PO7 = dataArray{:, 55};
res_table.W_PO3 = dataArray{:, 56};
res_table.W_POz = dataArray{:, 57};
res_table.W_PO4 = dataArray{:, 58};
res_table.W_PO8 = dataArray{:, 59};
res_table.W_P6 = dataArray{:, 60};
res_table.W_P2 = dataArray{:, 61};
res_table.W_CPz = dataArray{:, 62};
res_table.W_CP4 = dataArray{:, 63};
res_table.W_TP8 = dataArray{:, 64};
res_table.W_C6 = dataArray{:, 65};
res_table.W_C2 = dataArray{:, 66};
res_table.W_FC4 = dataArray{:, 67};
res_table.W_FT8 = dataArray{:, 68};
res_table.W_F6 = dataArray{:, 69};
res_table.W_AF8 = dataArray{:, 70};
res_table.W_AF4 = dataArray{:, 71};
res_table.W_F2 = dataArray{:, 72};
res_table.W_FCz = dataArray{:, 73};
res_table.A_Fp1 = dataArray{:, 74};
res_table.A_Fz = dataArray{:, 75};
res_table.A_F3 = dataArray{:, 76};
res_table.A_F7 = dataArray{:, 77};
res_table.A_FT9 = dataArray{:, 78};
res_table.A_FC5 = dataArray{:, 79};
res_table.A_FC1 = dataArray{:, 80};
res_table.A_C3 = dataArray{:, 81};
res_table.A_T7 = dataArray{:, 82};
res_table.A_TP9 = dataArray{:, 83};
res_table.A_CP5 = dataArray{:, 84};
res_table.A_CP1 = dataArray{:, 85};
res_table.A_Pz = dataArray{:, 86};
res_table.A_P3 = dataArray{:, 87};
res_table.A_P7 = dataArray{:, 88};
res_table.A_O1 = dataArray{:, 89};
res_table.A_Oz = dataArray{:, 90};
res_table.A_O2 = dataArray{:, 91};
res_table.A_P4 = dataArray{:, 92};
res_table.A_P8 = dataArray{:, 93};
res_table.A_TP10 = dataArray{:, 94};
res_table.A_CP6 = dataArray{:, 95};
res_table.A_CP2 = dataArray{:, 96};
res_table.A_Cz = dataArray{:, 97};
res_table.A_C4 = dataArray{:, 98};
res_table.A_T8 = dataArray{:, 99};
res_table.A_FT10 = dataArray{:, 100};
res_table.A_FC6 = dataArray{:, 101};
res_table.A_FC2 = dataArray{:, 102};
res_table.A_F4 = dataArray{:, 103};
res_table.A_F8 = dataArray{:, 104};
res_table.A_Fp2 = dataArray{:, 105};
res_table.A_AF7 = dataArray{:, 106};
res_table.A_AF3 = dataArray{:, 107};
res_table.A_F1 = dataArray{:, 108};
res_table.A_F5 = dataArray{:, 109};
res_table.A_FT7 = dataArray{:, 110};
res_table.A_FC3 = dataArray{:, 111};
res_table.A_C1 = dataArray{:, 112};
res_table.A_C5 = dataArray{:, 113};
res_table.A_TP7 = dataArray{:, 114};
res_table.A_CP3 = dataArray{:, 115};
res_table.A_P1 = dataArray{:, 116};
res_table.A_P5 = dataArray{:, 117};
res_table.A_PO7 = dataArray{:, 118};
res_table.A_PO3 = dataArray{:, 119};
res_table.A_POz = dataArray{:, 120};
res_table.A_PO4 = dataArray{:, 121};
res_table.A_PO8 = dataArray{:, 122};
res_table.A_P6 = dataArray{:, 123};
res_table.A_P2 = dataArray{:, 124};
res_table.A_CPz = dataArray{:, 125};
res_table.A_CP4 = dataArray{:, 126};
res_table.A_TP8 = dataArray{:, 127};
res_table.A_C6 = dataArray{:, 128};
res_table.A_C2 = dataArray{:, 129};
res_table.A_FC4 = dataArray{:, 130};
res_table.A_FT8 = dataArray{:, 131};
res_table.A_F6 = dataArray{:, 132};
res_table.A_AF8 = dataArray{:, 133};
res_table.A_AF4 = dataArray{:, 134};
res_table.A_F2 = dataArray{:, 135};
res_table.A_FCz = dataArray{:, 136};
res_table.State = dataArray{:, 137};
res_table.Vig = dataArray{:, 138};
res_table.Pup = dataArray{:, 139};
res_table.pcPup = dataArray{:, 140};
res_table.pcRT = dataArray{:, 141};
res_table.stimulus = dataArray{:, 142};
res_table.response = dataArray{:, 143};

res_table.W_all=nanmean(table2array(res_table(:,[11:73])),2);

res_table(res_table.SubID==334,:)=[];

% clean RTs
res_table(res_table.RT<0.3,:)=[];
% res_table(res_table.SubID==29,:)=[];
res_mat=table2array(res_table);

res_table.State=categorical(res_table.State);
res_table.SubID=categorical(res_table.SubID);
res_table.StimCat=categorical(res_table.StimCat);
res_table.Task=categorical(res_table.Task);

uniqueSubs=unique(Ppties_Waves(:,1));
uniqueSubs2=unique(res_table.SubID);
Ppties_Waves_av=[];
for nS=1:length(uniqueSubs)
    for nPr=1:60
        temp=nanmean(res_table.pcPup(res_table.SubID==uniqueSubs2(nS) & res_table.ProbeN==nPr & res_table.DistProbe>-20));
        Ppties_Waves_av=[Ppties_Waves_av ; [nanmean(Ppties_Waves(Ppties_Waves(:,1)==uniqueSubs(nS) & Ppties_Waves(:,2)==nPr,:)) temp]];
    end
end

%%
table_av=array2table(Ppties_Waves_av,'VariableNames',{'SubID','nPr','Elec','ProbeN','BlockN','Task','State','Vig','nWave','Dur','Freq','P2P','xNeg','xPos','endPos','xNegPeak','NegPeak','xPosPeak','PosPeak','DownSlope','UpSlope','MaxAmpW','MinAmpW','pcPup'});
table_av.State(table_av.State==4)=3;
table_av.Vig=abs(table_av.Vig-5);
table_av.SubID=categorical(table_av.SubID);
table_av.Task=categorical(table_av.Task);
table_av.State=categorical(table_av.State);

%%
myS=unique(table_av.SubID);
all_zval=[];
for nS=1:length(myS)
    for nTask=1:2
        all_zval=[all_zval ; [zscore(table_av.Vig(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            zscore(table_av.nWave(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            zscore(table_av.P2P(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            zscore(table_av.DownSlope(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask))) ...
            zscore(table_av.UpSlope(table_av.SubID==myS(nS) & table_av.Task==num2str(nTask)))]];
    end
end

Titles={'nWave','P2P','DownSlope','UpSlope'};
for nT=1:length(Titles)
    figure; set(gcf,'Position',[ 562   669   325   316]);
    simpleCorPlotbin(all_zval(:,1),all_zval(:,nT+1),{'o','k',[1 1 1]*0.5,288,10,4},'Pearson',0)
    format_fig;
    xlabel('Vig'); ylabel(Titles{nT});
    %export_fig([path_fig filesep 'LocalSleep_SWppties_CorrZvalAllE_' Titles{nT} '.fig'])
    %export_fig([path_fig filesep 'LocalSleep_SWppties_CorrZvalAllE_' Titles{nT} '.eps'],'-r 300')
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
% % %export_fig([path_fig filesep 'WanderIM_behav_fatigue_byStateAndTask.eps'],'-r 300')
%%

mdlall0= fitlme(table_av,sprintf('Vig~1+Task+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));

nWav_mdlall1= fitlme(table_av,sprintf('Vig~1+Task+nWave+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
nP2P_mdlall1= fitlme(table_av,sprintf('Vig~1+Task+P2P+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
nDSl_mdlall1= fitlme(table_av,sprintf('Vig~1+Task+DownSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
nUSl_mdlall1= fitlme(table_av,sprintf('Vig~1+Task+UpSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));


mdlallb0= fitlme(table_av,sprintf('pcPup~1+Task+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));

nWav_mdlallb1= fitlme(table_av,sprintf('pcPup~1+Task+nWave+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
nP2P_mdlallb1= fitlme(table_av,sprintf('pcPup~1+Task+P2P+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
nDSl_mdlallb1= fitlme(table_av,sprintf('pcPup~1+Task+DownSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
nUSl_mdlallb1= fitlme(table_av,sprintf('pcPup~1+Task+UpSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));


nWav_mdlav0= fitlme(table_av,sprintf('nWave~1+Task+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
nP2P_mdlav0= fitlme(table_av,sprintf('P2P~1+Task+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
nDSl_mdlav0= fitlme(table_av,sprintf('DownSlope~1+Task+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
nUSl_mdlav0= fitlme(table_av,sprintf('UpSlope~1+Task+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));

nWav_mdlav1= fitlme(table_av,sprintf('nWave~1+Task+State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
nP2P_mdlav1= fitlme(table_av,sprintf('P2P~1+Task+State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
nDSl_mdlav1= fitlme(table_av,sprintf('DownSlope~1+Task+State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
nUSl_mdlav1= fitlme(table_av,sprintf('UpSlope~1+Task+State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));

% %%
% mdlall0= fitlme(table_av,sprintf('pcPup~1+Task+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
% 
% nWav_mdlall1= fitlme(table_av,sprintf('Vig~1+Task+nWave+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
% nP2P_mdlall1= fitlme(table_av,sprintf('Vig~1+Task+P2P+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
% nDSl_mdlall1= fitlme(table_av,sprintf('Vig~1+Task+DownSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
% nUSl_mdlall1= fitlme(table_av,sprintf('Vig~1+Task+UpSlope+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
% 

%%
table_all.ElecC=table_all.Elec;
table_all.ElecC=categorical(table_all.ElecC);

nW_state_mdlall0= fitlme(table_all,sprintf('nWave~1+Task+ElecC+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
P2P_state_mdlall0= fitlme(table_all,sprintf('P2P~1+Task+ElecC+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
DS_state_mdlall0= fitlme(table_all,sprintf('DownSlope~1+Task+ElecC+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
US_state_mdlall0= fitlme(table_all,sprintf('UpSlope~1+Task+ElecC+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));


nW_state_mdlall1= fitlme(table_all,sprintf('nWave~1+Task+ElecC*State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
P2P_state_mdlall1= fitlme(table_all,sprintf('P2P~1+Task+ElecC*State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
DS_state_mdlall1= fitlme(table_all,sprintf('DownSlope~1+Task+ElecC*State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
US_state_mdlall1= fitlme(table_all,sprintf('UpSlope~1+Task+ElecC*State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));

nW_modelcomp=compare(nW_state_mdlall0,nW_state_mdlall1);
P2P_modelcomp=compare(P2P_state_mdlall0,P2P_state_mdlall1);
DS_modelcomp=compare(DS_state_mdlall0,DS_state_mdlall1);
US_modelcomp=compare(US_state_mdlall0,US_state_mdlall1);

table_all.State=reordercats(table_all.State,[2 1 3]);
nW_state_mdlall1b= fitlme(table_all,sprintf('nWave~1+Task+ElecC*State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
P2P_state_mdlall1b= fitlme(table_all,sprintf('P2P~1+Task+ElecC*State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
DS_state_mdlall1b= fitlme(table_all,sprintf('DownSlope~1+Task+ElecC*State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
US_state_mdlall1b= fitlme(table_all,sprintf('UpSlope~1+Task+ElecC*State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
table_all.State=reordercats(table_all.State,[2 1 3]);

%%

nW_state_mdlall0b= fitlme(table_all,sprintf('nWave~1+Task+ElecC+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
P2P_state_mdlall0b= fitlme(table_all,sprintf('P2P~1+Task+ElecC+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
DS_state_mdlall0b= fitlme(table_all,sprintf('DownSlope~1+Task+ElecC+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
US_state_mdlall0b= fitlme(table_all,sprintf('UpSlope~1+Task+ElecC+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));


nW_state_mdlall1b= fitlme(table_all,sprintf('nWave~1+Task+ElecC+Vig+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
P2P_state_mdlall1b= fitlme(table_all,sprintf('P2P~1+Task+ElecC+Vig+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
DS_state_mdlall1b= fitlme(table_all,sprintf('DownSlope~1+Task+ElecC+Vig+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
US_state_mdlall1b= fitlme(table_all,sprintf('UpSlope~1+Task+ElecC+Vig+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));

%%
figure; set(gcf,'Position',[680   662   560*1.33   316]);
subplot(1,4,1); hold on; format_fig;
temp=nW_state_mdlall1;
for nState=2:3
    line([1 1]*nState,[-1 1]*double(temp.Coefficients(nState+1,3))+double(temp.Coefficients(nState+1,2)),'Color',Colors(nState,:),'LineWidth',3);
    scatter(nState,double(temp.Coefficients(nState+1,2)),'MarkerEdgeColor',Colors(nState,:),'MarkerFaceColor',Colors(nState,:),'MarkerFaceAlpha',0.7,'LineWidth',3,'SizeData',400);
end
xlim([1.5 3.5])
ylim([-0.1 1]*0.12)
line(xlim,[0 0],'Color',[1 1 1]*.25,'LineStyle','--')
set(gca,'XTick',2:3,'XTickLabel',{'MW','MB'});
title('Density')

subplot(1,4,2); hold on; format_fig;
temp=P2P_state_mdlall1;
for nState=2:3
    line([1 1]*nState,[-1 1]*double(temp.Coefficients(nState+1,3))+double(temp.Coefficients(nState+1,2)),'Color',Colors(nState,:),'LineWidth',3);
    scatter(nState,double(temp.Coefficients(nState+1,2)),'MarkerEdgeColor',Colors(nState,:),'MarkerFaceColor',Colors(nState,:),'MarkerFaceAlpha',0.7,'LineWidth',3,'SizeData',400);
end
xlim([1.5 3.5])
ylim([-0.1 1]*1.2)
line(xlim,[0 0],'Color',[1 1 1]*.25,'LineStyle','--')
set(gca,'XTick',2:3,'XTickLabel',{'MW','MB'});
title('P2P')

subplot(1,4,3); hold on; format_fig;
temp=DS_state_mdlall1;
for nState=2:3
    line([1 1]*nState,[-1 1]*double(temp.Coefficients(nState+1,3))+double(temp.Coefficients(nState+1,2)),'Color',Colors(nState,:),'LineWidth',3);
    scatter(nState,double(temp.Coefficients(nState+1,2)),'MarkerEdgeColor',Colors(nState,:),'MarkerFaceColor',Colors(nState,:),'MarkerFaceAlpha',0.7,'LineWidth',3,'SizeData',400);
end
xlim([1.5 3.5])
ylim([-0.1 1]*12)
line(xlim,[0 0],'Color',[1 1 1]*.25,'LineStyle','--')
set(gca,'XTick',2:3,'XTickLabel',{'MW','MB'});
title('Down Slope')

subplot(1,4,4); hold on; format_fig;
temp=US_state_mdlall1;
for nState=2:3
    line([1 1]*nState,[-1 1]*double(temp.Coefficients(nState+1,3))+double(temp.Coefficients(nState+1,2)),'Color',Colors(nState,:),'LineWidth',3);
    scatter(nState,double(temp.Coefficients(nState+1,2)),'MarkerEdgeColor',Colors(nState,:),'MarkerFaceColor',Colors(nState,:),'MarkerFaceAlpha',0.7,'LineWidth',3,'SizeData',400);
end
xlim([1.5 3.5])
ylim([-0.1 1]*16)
line(xlim,[0 0],'Color',[1 1 1]*.25,'LineStyle','--')
set(gca,'XTick',2:3,'XTickLabel',{'MW','MB'});
title('Up Slope')
%export_fig([path_fig filesep 'LocalSleep_SWppties_LMEbyProbes_ModelRes_v2.fig'])
%export_fig([path_fig filesep 'LocalSleep_SWppties_LMEbyProbes_ModelRes_v2.eps'],'-r 300')

%%

for nE=1:63
    mdl1= fitlme(table_all(table_all.Elec==nE,:),sprintf('nWave~1+Task+State+(1|SubID)'));%,GO_table.Properties.VariableNames{10+nE}));
    mdl1_MWvsON(nE,:)=double(mdl1.Coefficients(3,2:end));
    mdl1_MBvsON(nE,:)=double(mdl1.Coefficients(4,2:end));
end

%%
States={'ON','MW','MB'};
for nState=1:3
    figure; set(gcf,'Position',[ 562   669   325   316]);
    temp_topo=grpstats(table_all.nWave(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
    temp_topo([10 21])=NaN;
    simpleTopoPlot_ft(temp_topo/20*60, path_PsychFTlayout,'on',[],0,1);
    hb=colorbar; caxis([3 8])
    set(hb,'YTick',[]);
    title(sprintf('D - %s',States{nState}))
    format_fig; %set(gca,'FontSize',26,'FontWeight','normal');
    %     export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo_PerState.fig'])
    export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo_' States{nState} '.eps'],'-r 300')
end


for nState=1:3
    figure; set(gcf,'Position',[ 562   669   325   316]);
    temp_topo=grpstats(table_all.P2P(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
    temp_topo([10 21])=NaN;
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    hb=colorbar; caxis([25 60])
    set(hb,'YTick',[]);
    title(sprintf('A - %s',States{nState}))
    format_fig; %set(gca,'FontSize',22);
    
    
    %     export_fig([path_fig filesep 'LocalSleep_SWppties_P2PTopo_' '.fig'])
    export_fig([path_fig filesep 'LocalSleep_SWppties_P2PTopo_' States{nState} '.eps'],'-r 300')
end

for nState=1:3
    figure; set(gcf,'Position',[ 562   669   325   316]);
    temp_topo=grpstats(table_all.DownSlope(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
    temp_topo([10 21])=NaN;
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    hb=colorbar; caxis([250 480])
    set(hb,'YTick',[]);
    title(sprintf('S - %s',States{nState}))
    format_fig; %set(gca,'FontSize',22);
    
    %     export_fig([path_fig filesep 'LocalSleep_SWppties_UpSlopeTopo_PerState.fig'])
    export_fig([path_fig filesep 'LocalSleep_SWppties_DownSlopeTopo_' States{nState} '.eps'],'-r 300')
end

for nState=1:3
    figure; set(gcf,'Position',[ 562   669   325   316]);
    temp_topo=grpstats(table_all.UpSlope(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
    temp_topo([10 21])=NaN;
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    hb=colorbar; caxis([250 530])
    set(hb,'YTick',[]);
    title(sprintf('S - %s',States{nState}))
    format_fig; %set(gca,'FontSize',22);
    
    %     export_fig([path_fig filesep 'LocalSleep_SWppties_UpSlopeTopo_PerState.fig'])
    export_fig([path_fig filesep 'LocalSleep_SWppties_UpSlopeTopo_' States{nState} '.eps'],'-r 300')
end

%% Contrast on slow-waves properties (Rev1 - Comment R1)
States={'ON','MW','MB'};
figure; set(gcf,'Position',[ 562   669   325*3   316*2]);
count=0;
for nState=1:2
    for nState2=nState+1:3
        count=count+1;
        subplot(2,2,count)
        temp_topo1=grpstats(table_all.nWave(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
        temp_topo2=grpstats(table_all.nWave(table_all.State==num2str(nState2)),table_all.Elec(table_all.State==num2str(nState2))); %in seconds but for all probes (60) so equivalent of in minutes
        
        temp_topo=temp_topo2-temp_topo1;
        temp_topo([10 21])=NaN;
        simpleTopoPlot_ft(temp_topo/20*60, path_PsychFTlayout,'on',[],0,1);
        hb=colorbar; 
        caxis([-1 1])
%         set(hb,'YTick',[]);
        title(sprintf('D - %s vs %s',States{nState2},States{nState}))
        format_fig; %set(gca,'FontSize',26,'FontWeight','normal');
        %     export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo_PerState.fig'])
        %     export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo_' States{nState} '.eps'],'-r 300')
    end
end

figure; set(gcf,'Position',[ 562   669   325*3   316*2]);
count=0;
for nState=1:2
    for nState2=nState+1:3
        count=count+1;
        subplot(2,2,count)
        temp_topo1=grpstats(table_all.P2P(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
        temp_topo2=grpstats(table_all.P2P(table_all.State==num2str(nState2)),table_all.Elec(table_all.State==num2str(nState2))); %in seconds but for all probes (60) so equivalent of in minutes
        
        temp_topo=temp_topo2-temp_topo1;
        temp_topo([10 21])=NaN;
        simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
        hb=colorbar; 
        caxis([-1 1]*15)
%         set(hb,'YTick',[]);
        title(sprintf('A - %s vs %s',States{nState2},States{nState}))
        format_fig; %set(gca,'FontSize',26,'FontWeight','normal');
        %     export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo_PerState.fig'])
        %     export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo_' States{nState} '.eps'],'-r 300')
    end
end

figure; set(gcf,'Position',[ 562   669   325*3   316*2]);
count=0;
for nState=1:2
    for nState2=nState+1:3
        count=count+1;
        subplot(2,2,count)
        temp_topo1=grpstats(table_all.DownSlope(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
        temp_topo2=grpstats(table_all.DownSlope(table_all.State==num2str(nState2)),table_all.Elec(table_all.State==num2str(nState2))); %in seconds but for all probes (60) so equivalent of in minutes
        
        temp_topo=temp_topo2-temp_topo1;
        temp_topo([10 21])=NaN;
        simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
        hb=colorbar; 
        caxis([-1 1]*250)
%         set(hb,'YTick',[]);
        title(sprintf('DS - %s vs %s',States{nState2},States{nState}))
        format_fig; %set(gca,'FontSize',26,'FontWeight','normal');
        %     export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo_PerState.fig'])
        %     export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo_' States{nState} '.eps'],'-r 300')
    end
end

figure; set(gcf,'Position',[ 562   669   325*3   316*2]);
count=0;
for nState=1:2
    for nState2=nState+1:3
        count=count+1;
        subplot(2,2,count)
        temp_topo1=grpstats(table_all.UpSlope(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
        temp_topo2=grpstats(table_all.UpSlope(table_all.State==num2str(nState2)),table_all.Elec(table_all.State==num2str(nState2))); %in seconds but for all probes (60) so equivalent of in minutes
        
        temp_topo=temp_topo2-temp_topo1;
        temp_topo([10 21])=NaN;
        simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
        hb=colorbar; 
        caxis([-1 1]*200)
%         set(hb,'YTick',[]);
        title(sprintf('US - %s vs %s',States{nState2},States{nState}))
        format_fig; %set(gca,'FontSize',26,'FontWeight','normal');
        %     export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo_PerState.fig'])
        %     export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo_' States{nState} '.eps'],'-r 300')
    end
end
%%
MarkersT={'d','o'};
plotNames={'D','A','S'};
data=[];
XlimPlot=[1 3; 20 60; 200 450];
    figure; set(gcf,'Position',[ 440   315   800   360]);
for nvar=1:3
    for i = 1:3
        if nvar==1
        temp=table_all.nWave(table_all.State==num2str(i) & table_all.ElecC=='63');
        elseif nvar==2
        temp=table_all.P2P(table_all.State==num2str(i) & table_all.ElecC=='63');
        elseif nvar==3
        temp=table_all.DownSlope(table_all.State==num2str(i) & table_all.ElecC=='63');
        end
        tempS=double(table_all.SubID(table_all.State==num2str(i) & table_all.ElecC=='63'));
        tempbyS=[]; myS=unique(tempS);
        tempbyS_n=[];
        for nS=1:length(myS)
            tempbyS(nS)=nanmean(temp(tempS==myS(nS)));
            tempbyS_n(nS)=sum((tempS==myS(nS)));
        end
        data{i,nvar} = tempbyS;
        data_n{i,nvar} = tempbyS_n;
    end
    
    subplot(1,3,nvar)
    h1 = raincloud_plot(data{1,nvar}, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
        'box_col_match', 0,'size_data',data_n{1,nvar});
    h1{2}.MarkerFaceAlpha=0.5; h1{2}.Marker=MarkersT{1};
    %         h1{7}.MarkerFaceAlpha=0.5; h1{7}.Marker=MarkersT{ntask};  h1{7}.MarkerFaceColor=[1 1 1]*0.5;   h1{7}.MarkerEdgeColor=Colors(1,:);    h1{7}.SizeData=144;
    h2 = raincloud_plot(data{2,nvar}, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'size_data',data_n{2,nvar});
    h2{2}.MarkerFaceAlpha=0.5; h2{2}.Marker=MarkersT{1};
    h3 = raincloud_plot(data{3,nvar}, 'box_on', 1, 'color', Colors(3,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,'size_data',data_n{3,nvar});
    h3{2}.MarkerFaceAlpha=0.5; h3{2}.Marker=MarkersT{1};
    
    if nvar==1
        set(gca,'XLim', XlimPlot(nvar,:),'YTick',''); xlabel(plotNames{nvar});
    elseif nvar==2
        set(gca,'XLim', XlimPlot(nvar,:),'YTick',''); xlabel(plotNames{nvar});
    elseif nvar==3
        set(gca,'XLim', XlimPlot(nvar,:),'YTick',''); xlabel(plotNames{nvar});
    end
    box off
    format_fig;
end

%%
figure;
set(gcf,'Position',[680   657   733   321]);
subplot(1,2,1);
temp_topo=mean(thr_Wave,1);
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'on',[],0,1);
title('Mean Threshold (\muV)')
format_fig;
hb=colorbar;
set(hb,'Position',[0.92/2 0.6 0.035 0.25])

subplot(1,2,2);
pption1SD=[];
for k=1:63
    pption1SD(k)=mean(thr_Wave(:,k)>(mean(thr_Wave(:,k))-std(thr_Wave(:,k))) & thr_Wave(:,k)<(mean(thr_Wave(:,k))+std(thr_Wave(:,k))));
end

temp_topo=std(thr_Wave,1)./mean(thr_Wave,1);
simpleTopoPlot_ft(temp_topo', path_PsychFTlayout,'on',[],0,1);
title('SD Threshold (\muV)')
format_fig;
hb=colorbar;
set(hb,'Position',[0.92 0.6 0.035 0.25])
% ylabel(hb, '\muV')
export_fig([path_fig filesep 'LocalSleep_SWppties_Threshold.fig'])
export_fig([path_fig filesep 'LocalSleep_SWppties_Threshold.png'],'-r 300')

