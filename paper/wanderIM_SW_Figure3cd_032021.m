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

%%
figure; set(gcf,'Position',[ 562   669   325   316]);
% subplot(1,3,1)
temp_topo=grpstats(table_all.nWave,table_all.Elec); %in seconds but for all probes (60) so equivalent of in minutes
temp_topo([10 21])=NaN;
temp_topo_Dens=temp_topo'/20*60;
simpleTopoPlot_ft(temp_topo/20*60, path_PsychFTlayout,'on',[],0,1);
hb=colorbar; %caxis([0 9])
caxis([3 8])
set(hb,'YTick',[]);
title('SW density /min')
format_fig;
% export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo.fig'])
% export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo.eps'],'-r 300')

figure; set(gcf,'Position',[ 562   669   325   316]);
% subplot(1,3,1)
temp_topo=grpstats(table_all.P2P,table_all.Elec); %in seconds but for all probes (60) so equivalent of in minutes
temp_topo([10 21])=NaN;
temp_topo_P2P=temp_topo';
simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
hb=colorbar; caxis([25 60])
set(hb,'YTick',[]);

title('P2P amp \muV')
format_fig;
% export_fig([path_fig filesep 'LocalSleep_SWppties_P2PTopo.fig'])
% export_fig([path_fig filesep 'LocalSleep_SWppties_P2PTopo.eps'],'-r 300')

figure; set(gcf,'Position',[ 562   669   325   316]);
% subplot(1,3,1)
temp_topo=grpstats(table_all.DownSlope,table_all.Elec); %in seconds but for all probes (60) so equivalent of in minutes
temp_topo([10 21])=NaN;
temp_topo_DS=temp_topo';
simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
hb=colorbar; caxis([250 480])
set(hb,'YTick',[]);
title('Down Slope \muV.s^-^1')
format_fig;
% export_fig([path_fig filesep 'LocalSleep_SWppties_DownSlopeTopo.fig'])
% export_fig([path_fig filesep 'LocalSleep_SWppties_DownSlopeTopo.eps'],'-r 300')


figure; set(gcf,'Position',[ 562   669   325   316]);
% subplot(1,3,1)
temp_topo=grpstats(table_all.UpSlope,table_all.Elec); %in seconds but for all probes (60) so equivalent of in minutes
temp_topo([10 21])=NaN;
temp_topo_US=temp_topo';
simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
hb=colorbar; caxis([250 530])
set(hb,'YTick',[]);
title('Up Slope \muV.s^-^1')
format_fig;
% export_fig([path_fig filesep 'LocalSleep_SWppties_UpSlopeTopo.fig'])
% export_fig([path_fig filesep 'LocalSleep_SWppties_UpSlopeTopo.eps'],'-r 300')

%%
States={'ON','MW','MB'};
for nState=1:3
    figure; set(gcf,'Position',[ 562   669   325   316]);
    temp_topo=grpstats(table_all.nWave(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
    temp_topo([10 21])=NaN;
    temp_topo_Dens_byS(nState,:)=temp_topo';
    simpleTopoPlot_ft(temp_topo/20*60, path_PsychFTlayout,'on',[],0,1);
    hb=colorbar; caxis([3 8])
    set(hb,'YTick',[]);
    title(sprintf('D - %s',States{nState}))
    format_fig; %set(gca,'FontSize',26,'FontWeight','normal');
    %     export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo_PerState.fig'])
%     export_fig([path_fig filesep 'LocalSleep_SWppties_DensityTopo_' States{nState} '.eps'],'-r 300')
end


for nState=1:3
    figure; set(gcf,'Position',[ 562   669   325   316]);
    temp_topo=grpstats(table_all.P2P(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
    temp_topo([10 21])=NaN;
    temp_topo_P2P_byS(nState,:)=temp_topo';
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    hb=colorbar; caxis([25 60])
    set(hb,'YTick',[]);
    title(sprintf('A - %s',States{nState}))
    format_fig; %set(gca,'FontSize',22);
    
    
    %     export_fig([path_fig filesep 'LocalSleep_SWppties_P2PTopo_' '.fig'])
%     export_fig([path_fig filesep 'LocalSleep_SWppties_P2PTopo_' States{nState} '.eps'],'-r 300')
end

for nState=1:3
    figure; set(gcf,'Position',[ 562   669   325   316]);
    temp_topo=grpstats(table_all.DownSlope(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
    temp_topo([10 21])=NaN;
    temp_topo_DS_byS(nState,:)=temp_topo';
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    hb=colorbar; caxis([250 480])
    set(hb,'YTick',[]);
    title(sprintf('S - %s',States{nState}))
    format_fig; %set(gca,'FontSize',22);
    
    %     export_fig([path_fig filesep 'LocalSleep_SWppties_UpSlopeTopo_PerState.fig'])
%     export_fig([path_fig filesep 'LocalSleep_SWppties_DownSlopeTopo_' States{nState} '.eps'],'-r 300')
end

for nState=1:3
    figure; set(gcf,'Position',[ 562   669   325   316]);
    temp_topo=grpstats(table_all.UpSlope(table_all.State==num2str(nState)),table_all.Elec(table_all.State==num2str(nState))); %in seconds but for all probes (60) so equivalent of in minutes
    temp_topo([10 21])=NaN;
    temp_topo_US_byS(nState,:)=temp_topo';
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    hb=colorbar; caxis([250 530])
    set(hb,'YTick',[]);
    title(sprintf('S - %s',States{nState}))
    format_fig; %set(gca,'FontSize',22);
    
    %     export_fig([path_fig filesep 'LocalSleep_SWppties_UpSlopeTopo_PerState.fig'])
%     export_fig([path_fig filesep 'LocalSleep_SWppties_UpSlopeTopo_' States{nState} '.eps'],'-r 300')
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
% export_fig([path_fig filesep 'LocalSleep_SWppties_Threshold.fig'])
% export_fig([path_fig filesep 'LocalSleep_SWppties_Threshold.png'],'-r 300')

