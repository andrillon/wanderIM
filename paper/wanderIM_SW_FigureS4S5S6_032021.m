%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
% close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
% addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_ica'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'probe_infEEG_S3*.mat']);

%% Compute monte cartlo p
clus_alpha=0.025;
montecarlo_alpha=0.05/12/4;

addpath((path_fieldtrip)); ft_defaults;
cfg_neighb=[];
cfg_neighb.method = 'tri';
cfg_neighb.layout=path_PsychFTlayout;
neighbours = ft_prepare_neighbours(cfg_neighb);

for nWin=1:4
    load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsON_nWa_Win' num2str(nWin) '.mat'])
    load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsON_P2P_Win' num2str(nWin) '.mat'])
    load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsON_DwS_Win' num2str(nWin) '.mat'])
    load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsON_UpS_Win' num2str(nWin) '.mat'])
    load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MBvsON_nWa_Win' num2str(nWin) '.mat'])
    load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MBvsON_P2P_Win' num2str(nWin) '.mat'])
    load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MBvsON_DwS_Win' num2str(nWin) '.mat'])
    load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MBvsON_UpS_Win' num2str(nWin) '.mat'])
    load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsMB_nWa_Win' num2str(nWin) '.mat'])
    load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsMB_P2P_Win' num2str(nWin) '.mat'])
    load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsMB_DwS_Win' num2str(nWin) '.mat'])
    load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsMB_UpS_Win' num2str(nWin) '.mat'])
    
    
    [MWvsON_nWa_clus{nWin}]=get_clusterperm_lme(MWvsON_nWa_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
    [MWvsON_P2P_clus{nWin}]=get_clusterperm_lme(MWvsON_P2P_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
    [MWvsON_DwS_clus{nWin}]=get_clusterperm_lme(MWvsON_DwS_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
    [MWvsON_UpS_clus{nWin}]=get_clusterperm_lme(MWvsON_UpS_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
    
    [MBvsON_nWa_clus{nWin}]=get_clusterperm_lme(MBvsON_nWa_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
    [MBvsON_P2P_clus{nWin}]=get_clusterperm_lme(MBvsON_P2P_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
    [MBvsON_DwS_clus{nWin}]=get_clusterperm_lme(MBvsON_DwS_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
    [MBvsON_UpS_clus{nWin}]=get_clusterperm_lme(MBvsON_UpS_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
    
    [MWvsMB_nWa_clus{nWin}]=get_clusterperm_lme(MWvsMB_nWa_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
    [MWvsMB_P2P_clus{nWin}]=get_clusterperm_lme(MWvsMB_P2P_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
    [MWvsMB_DwS_clus{nWin}]=get_clusterperm_lme(MWvsMB_DwS_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
    [MWvsMB_UpS_clus{nWin}]=get_clusterperm_lme(MWvsMB_UpS_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
    
end
%%
% for nTask=1:2
Titles={'Density','P2P','DownSlope','UpSlope'};
AxesLim={[-5 5],[-5 5],[-5 5],[-5 5]};
cmap=cbrewer('div','RdBu',256); cmap=flipud(cmap);
for nP=1:4
    myplots(nP)=figure;
    set(myplots(nP),'Position',[5         637        1386         341]);
    tight_subplot(1,4);
end
for nP=1:4
    for nWin=1:4
        load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsON_nWa_Win' num2str(nWin) '.mat'])
        load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsON_P2P_Win' num2str(nWin) '.mat'])
        load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsON_DwS_Win' num2str(nWin) '.mat'])
        load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsON_UpS_Win' num2str(nWin) '.mat'])
        
        figure(myplots(nP));
        subplot(1,4,nWin);
        if nP==1
            temp_topo=MWvsON_nWa_est{1}(:,3);
            temp_clus=MWvsON_nWa_clus{nWin};
            
            DensWinTopo(nWin,:)=temp_topo';
        elseif nP==2
            temp_topo=MWvsON_P2P_est{1}(:,3);
            temp_clus=MWvsON_P2P_clus{nWin};
            
            
            P2PWinTopo(nWin,:)=temp_topo';
        elseif nP==3
            temp_topo=MWvsON_DwS_est{1}(:,3);
            temp_clus=MWvsON_DwS_clus{nWin};
            
            
            DSWinTopo(nWin,:)=temp_topo';
        elseif nP==4
            temp_topo=MWvsON_UpS_est{1}(:,3);
            temp_clus=MWvsON_UpS_clus{nWin};
            
            
            USWinTopo(nWin,:)=temp_topo';
        end
        
        simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
%         colorbar; 
        caxis(AxesLim{nP}); %caxis([-1 1]*max(abs(temp_topo)))
        colormap(cmap);
        
        load(path_PsychFTlayout);
        ft_plot_lay_me(layout,'chanindx',1:63,'pointsymbol','.','pointcolor',[1 1 1]*.5,'pointsize',64,'box','no','label','no')
        if ~isempty(temp_clus)
            for nclus=1:length(temp_clus)
                ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
                fprintf('... ... found %s cluster (%g) of %g electrodes (tval cluster=%g, Pmc=%g)\n',temp_clus{nclus}{1},nclus,length(temp_clus{nclus}{2}),temp_clus{nclus}{3},temp_clus{nclus}{4})
            end
        end
        title(['MW vs ON: ' Titles{nP}])
        %     export_fig([path_fig filesep 'LocalSleep_MWvsON_' Titles{nP} '.fig'])
        
    end
%             export_fig([path_fig filesep 'LocalSleep_MWvsON_' Titles{nP} '_Win.eps'],'-r 200')
end

%%
for nP=1:4
    myplots(nP)=figure;
    set(myplots(nP),'Position',[5         637        1386         341]);
    tight_subplot(1,4);
end
for nP=1:4
    for nWin=1:4
        load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MBvsON_nWa_Win' num2str(nWin) '.mat'])
        load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MBvsON_P2P_Win' num2str(nWin) '.mat'])
        load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MBvsON_DwS_Win' num2str(nWin) '.mat'])
        load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MBvsON_UpS_Win' num2str(nWin) '.mat'])
        
        figure(myplots(nP));
        subplot(1,4,nWin);
        if nP==1
            temp_topo=MBvsON_nWa_est{1}(:,3);
            temp_clus=MBvsON_nWa_clus{nWin};
                          
            
            DensWinTopo(nWin,:)=temp_topo';
        elseif nP==2
            temp_topo=MBvsON_P2P_est{1}(:,3);
            temp_clus=MBvsON_P2P_clus{nWin};
                       
            
            P2PWinTopo(nWin,:)=temp_topo';
        elseif nP==3
            temp_topo=MBvsON_DwS_est{1}(:,3);
            temp_clus=MBvsON_DwS_clus{nWin};
                  
            
            DSWinTopo(nWin,:)=temp_topo';
        elseif nP==4
            temp_topo=MBvsON_UpS_est{1}(:,3);
            temp_clus=MBvsON_UpS_clus{nWin};
                 
            
            USWinTopo(nWin,:)=temp_topo';   
        end
        
        simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
        caxis(AxesLim{nP}); %caxis([-1 1]*max(abs(temp_topo)))
        colormap(cmap);
        
        load(path_PsychFTlayout);
        ft_plot_lay_me(layout,'chanindx',1:63,'pointsymbol','.','pointcolor',[1 1 1]*.5,'pointsize',64,'box','no','label','no')
        if ~isempty(temp_clus)
            for nclus=1:length(temp_clus)
                ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
                fprintf('... ... found %s cluster (%g) of %g electrodes (tval cluster=%g, Pmc=%g)\n',temp_clus{nclus}{1},nclus,length(temp_clus{nclus}{2}),temp_clus{nclus}{3},temp_clus{nclus}{4})
            end
        end
        title(['MB vs ON: ' Titles{nP}])
        %     export_fig([path_fig filesep 'LocalSleep_MWvsON_' Titles{nP} '.fig'])
        %     export_fig([path_fig filesep 'LocalSleep_MWvsON_' Titles{nP} '.eps'],'-r 300')
        
    end
%                export_fig([path_fig filesep 'LocalSleep_MBvsON_' Titles{nP} '_Win.eps'],'-r 200')
   end

%%
for nP=1:4
    myplots(nP)=figure;
    set(myplots(nP),'Position',[5         637        1386         341]);
    tight_subplot(1,4);
end
for nP=1:4
    for nWin=1:4
        load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsMB_nWa_Win' num2str(nWin) '.mat'])
        load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsMB_P2P_Win' num2str(nWin) '.mat'])
        load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsMB_DwS_Win' num2str(nWin) '.mat'])
        load(['/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsMB_UpS_Win' num2str(nWin) '.mat'])
        
        figure(myplots(nP));
        subplot(1,4,nWin);
        if nP==1
            temp_topo=MWvsMB_nWa_est{1}(:,3);
            temp_clus=MWvsMB_nWa_clus{nWin};
        
                    DensWinTopo(nWin,:)=temp_topo';
elseif nP==2
            temp_topo=MWvsMB_P2P_est{1}(:,3);
            temp_clus=MWvsMB_P2P_clus{nWin};
       
                    P2PWinTopo(nWin,:)=temp_topo';
elseif nP==3
            temp_topo=MWvsMB_DwS_est{1}(:,3);
            temp_clus=MWvsMB_DwS_clus{nWin};
        
                    DSWinTopo(nWin,:)=temp_topo';
elseif nP==4
            temp_topo=MWvsMB_UpS_est{1}(:,3);
            temp_clus=MWvsMB_UpS_clus{nWin};
        
                    USWinTopo(nWin,:)=temp_topo';
end
        
        simpleTopoPlot_ft(-temp_topo, path_PsychFTlayout,'on',[],0,1);
        caxis(AxesLim{nP}); %caxis([-1 1]*max(abs(temp_topo)))
        colormap(cmap);
        
        load(path_PsychFTlayout);
        ft_plot_lay_me(layout,'chanindx',1:63,'pointsymbol','.','pointcolor',[1 1 1]*.5,'pointsize',64,'box','no','label','no')
        if ~isempty(temp_clus)
            for nclus=1:length(temp_clus)
                ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
                fprintf('... ... found %s cluster (%g) of %g electrodes (tval cluster=%g, Pmc=%g)\n',temp_clus{nclus}{1},nclus,length(temp_clus{nclus}{2}),temp_clus{nclus}{3},temp_clus{nclus}{4})
            end
        end
        title(['MB vs MW: ' Titles{nP}])
        %     export_fig([path_fig filesep 'LocalSleep_MWvsON_' Titles{nP} '.fig'])
        %     export_fig([path_fig filesep 'LocalSleep_MWvsON_' Titles{nP} '.eps'],'-r 300')
        
    end
%             export_fig([path_fig filesep 'LocalSleep_MBvsMW_' Titles{nP} '_Win.eps'],'-r 200')
end