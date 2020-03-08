%%
clear all
close all

run ../localdef_wanderIM;

addpath(genpath(lscpTools_path));
addpath(edf2mat_path);
eyet_path=[root_path filesep 'eyetracker'];
behav_path=[root_path filesep 'behav'];

files=dir([eyet_path filesep 'wanderIM_eyelink_S*_clean.mat']);

%% Loop on files
all_blinks_mat=[];
all_sacc_mat=[];
all_pup_mat=nan(20*6*10,25019);
all_PosX_mat=nan(20*6*10,25019);
all_PosY_mat=nan(20*6*10,25019);
all_eyet_probes=[];
ntotprobe=0;
for n=1:length(files)
    subID=files(n).name;
    bound=findstr(subID,'wanderIM_eyelink_S');
    subID=subID(length('wanderIM_eyelink_S')+(1:3));
    
    % load behavioural results
    behav_file=dir([behav_path filesep 'wanderIM_behavres_s' subID '_*.mat']);
    if isempty(behav_file)
continue;
    end
    load([behav_path filesep behav_file.name]);
    
    % load eye-tracker data
    savename=['wanderIM_eyelink_S' subID '_clean'];
    load([eyet_path filesep savename])
    fprintf('... %s loaded\n',subID)
    Fs=EL_headers.Fs;
    xTime=(-20*Fs:5*Fs)/Fs;

    probe_times=nan(6,10);
    probe_endtimes=nan(6,10);
    for npr=1:10
        probe_times(:,npr)=EL_events.Events.time(match_str(EL_events.Events.type,sprintf('P%g',npr)))';
        probe_endtimes(:,npr)=EL_events.Events.time(match_str(EL_events.Events.type,sprintf('EP%g',npr)))';
    end
    for nbl=1:6
        these_probes=probe_res(probe_res(:,4)==nbl,:);
        these_trials=test_res(test_res(:,1)==nbl,:);
        
        for npr=1:10
            fprintf('... Block %g Probe %g\n',nbl,npr)
            this_pr_tridx=these_probes(npr,6);
            ntotprobe=ntotprobe+1;
            if npr==1
                last_pr_tridx=0;
            else
                last_pr_tridx=these_probes(npr-1,6);
            end
            probe_details(1)=these_probes(npr,31); % look
            probe_details(2)=these_probes(npr,32); % state
            probe_details(3)=these_probes(npr,33); % origin
            probe_details(4)=these_probes(npr,34); % Awareness
            probe_details(5)=these_probes(npr,35); % Intention
            probe_details(6)=these_probes(npr,36); % Engagement
            probe_details(7)=these_probes(npr,37); % Performance
            probe_details(8)=these_probes(npr,38); % Alterness
            %  7: resp Q1 Looking (1: Yes / 2: No)
            %  8: resp Q2 Mind-State (1: On / 2: MW / 3: MB / 4: ?)
            %  9: resp Q3 Origine (1: room / 2: personal / 3: task)
            % 10: resp Q4 Awareness (1 (fully) to 4 (not at all))
            % 11: resp Q5 Intention (1 (fully) to 4 (not at all))
            % 12: resp Q6 Engagement (1 (not) to 4 (very))
            % 13: resp Q7 Performance (1 (bad) to 4 (good))
            % 14: resp Q8 Alterness (1 (alert) to 4 (sleepy))
            
            temp_testresgo=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,12)),:);
            tcorr_go=(temp_testresgo(end-19:end,12))';%/corr_go(n,these_probes(npr,5));
            temp_testresnogo=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,11)),:);
            tcorr_nogo=(temp_testresnogo(end-1:end,11))';%/corr_go(n,these_probes(npr,5));
            
            temp_corr_go=nanmean(tcorr_go);
            temp_corr_nogo=nanmean(tcorr_nogo);
            temp_rt_go=nanmean(temp_testresgo(:,10)-temp_testresgo(:,8));
            temp_corr_go2=tcorr_go;
            temp_corr_nogo2=tcorr_nogo;
            [dprime, crit]=calc_dprime(temp_corr_go2,temp_corr_nogo2==0);
            % Colum order
            temp_behav=[n nbl npr these_probes(npr,5) this_pr_tridx probe_details temp_corr_go temp_corr_nogo temp_rt_go  dprime crit];
            
            % Eye-Tracker - blinks
            this_probe_time=probe_times(nbl,npr);
            this_time_win=(-20*Fs:0)+this_probe_time;
            all_blinks_idx=ismember([EL_events.Blinks.start],this_time_win);
            temp_blinks=[EL_events.Blinks.start(all_blinks_idx)' EL_events.Blinks.duration(all_blinks_idx)'];
            all_blinks_mat=[all_blinks_mat ; [repmat(temp_behav,size(temp_blinks,1),1) temp_blinks]];
            
            % Eye-Tracker - saccades
            all_sacc_idx=ismember([EL_events.Sacc.start],this_time_win);
            temp_sacc=[EL_events.Sacc.start(all_sacc_idx)' EL_events.Sacc.duration(all_sacc_idx)' EL_events.Sacc.posX_start(all_sacc_idx)' EL_events.Sacc.posX_end(all_sacc_idx)' EL_events.Sacc.posY_start(all_sacc_idx)' EL_events.Sacc.posY_end(all_sacc_idx)' EL_events.Sacc.velo(all_sacc_idx)'];
            all_sacc_mat=[all_sacc_mat ; [repmat(temp_behav,size(temp_sacc,1),1) temp_sacc]];
            
            % Eye-Tracker - pupil
            this_probe_time=probe_times(nbl,npr);
            [~,idx_this_probe]=findclosest(EL_data.time,this_probe_time);
            temp_Pupil=EL_data.filt_pupilSize((-20*Fs:5*Fs)+idx_this_probe);
            all_pup_mat(ntotprobe,:)=[temp_behav temp_Pupil'];
            temp_Pupil_bef=nanmean(temp_Pupil(xTime>-4 & xTime<0));
            temp_Pupil_after=nanmean(temp_Pupil(xTime>0 & xTime<4))-temp_Pupil_bef;
            
%             % Eye-Tracker - pos X, Y
%             temp_PosX=EL_data.clean_posX((-20*Fs:5*Fs)+idx_this_probe);
%             all_PosX_mat(ntotprobe,:)=[temp_behav temp_PosX'];
%             temp_PosY=EL_data.clean_posY((-20*Fs:5*Fs)+idx_this_probe);
%             all_PosY_mat(ntotprobe,:)=[temp_behav temp_PosY'];
%             
            
            ampSac=sqrt((EL_events.Sacc.posX_end(all_sacc_idx)+EL_events.Sacc.posX_start(all_sacc_idx)).^2+(EL_events.Sacc.posY_end(all_sacc_idx)+EL_events.Sacc.posY_start(all_sacc_idx)).^2);
            all_eyet_probes=[all_eyet_probes ; [temp_behav sum(temp_blinks(:,end)<400) mean(temp_blinks(temp_blinks(:,end)<400,end)) size(temp_sacc,1) mean(ampSac) mean(temp_sacc(:,end)) temp_Pupil_bef temp_Pupil_after]];
        end
    end
end

%%
all_eyet_probes2=all_eyet_probes;
mysubs=unique(all_eyet_probes(:,1));
for ns=1:length(mysubs)
    for ntask=1:2
    all_eyet_probes2(all_eyet_probes2(:,1)==mysubs(ns) & all_eyet_probes2(:,4)==ntask,end-1)=nanzscore(all_eyet_probes2(all_eyet_probes2(:,1)==mysubs(ns) & all_eyet_probes2(:,4)==ntask,end-1));
    all_eyet_probes2(all_eyet_probes2(:,1)==mysubs(ns) & all_eyet_probes2(:,4)==ntask,end)=nanzscore(all_eyet_probes2(all_eyet_probes2(:,1)==mysubs(ns) & all_eyet_probes2(:,4)==ntask,end));
    end
end
%% figure Pupil
figure; set(gcf,'position',[4    12   542   762])
for task=1:2
    subplot(2,2,2*(task-1)+1);
    format_fig;
    for nstate=1:3
        tempplot=nanmean(all_pup_mat(all_pup_mat(:,7)==nstate & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
                tempplot=tempplot-nanmean(tempplot(xTime>0 & xTime<3));
        simpleTplot(xTime,tempplot,0,Colors(nstate,:),0,'-',1,0.5,500,1,2);
    end
    legend({'ON','MW','MB'})
    
    subplot(2,2,2*(task-1)+2);
    format_fig;
    for nstate=1:3
        tempplot=nanmean(all_pup_mat(all_pup_mat(:,7)==nstate & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
        tempplot=tempplot-nanmean(tempplot(xTime>-2 & xTime<0));
        simpleTplot(xTime,tempplot,0,Colors(nstate,:),0,'-',1,0.5,500,1,2);
        xlim([-2 5])
    end
    legend({'ON','MW','MB'})
end

figure
format_fig;
for task=1:2
    for nstate=1:3
        tempplot=(all_pup_mat(all_pup_mat(:,7)==nstate & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
        tempplot=nanmean(tempplot(:,xTime>-5 & xTime<0),2);
        if task==1
            simpleBarPlot(nstate+0.2*(task*2-3),tempplot,[Colors(nstate,:)],0.35,'k',[],2);
        elseif task==2
            simpleBarPlot(nstate+0.2*(task*2-3),tempplot,[1 1 1;Colors(nstate,:)],0.35,'k',[],2);
        end
    end
end
ylabel('Pupil before')
ylim([400 1400])
xlim([0.2 3.8])
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
xlim([0.2 3.8])

figure
format_fig;
for task=1:2
    for nstate=1:3
        tempplot=(all_pup_mat(all_pup_mat(:,7)==nstate & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
        tempplot=tempplot-nanmean(tempplot(:,xTime>-2 & xTime<0),2);
        tempplot=nanmean(tempplot(:,xTime>0 & xTime<4),2);
        if task==1
            simpleBarPlot(nstate+0.2*(task*2-3),tempplot,[Colors(nstate,:)],0.35,'k',[],2);
        elseif task==2
            simpleBarPlot(nstate+0.2*(task*2-3),tempplot,[1 1 1;Colors(nstate,:)],0.35,'k',[],2);
        end
    end
end
ylabel('Pupil after')
% ylim([400 1400])
xlim([0.2 3.8])
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
xlim([0.2 3.8])

%%
figure;
subplot(1,2,1)
format_fig;
for task=1:2
    for nstate=1:3
        tempplot=all_eyet_probes2(all_eyet_probes2(:,7)==nstate & all_eyet_probes2(:,4)==task,end-1);
        if task==1
            simpleBarPlot(nstate+0.2*(task*2-3),tempplot,[Colors(nstate,:)],0.35,'k',[],2);
        elseif task==2
            simpleBarPlot(nstate+0.2*(task*2-3),tempplot,[1 1 1;Colors(nstate,:)],0.35,'k',[],2);
        end
    end
end
ylabel('Pupil before')
% ylim([400 1400])
xlim([0.2 3.8])
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
xlim([0.2 3.8])

subplot(1,2,2)
format_fig;
for task=1:2
    for nstate=1:3
        tempplot=all_eyet_probes2(all_eyet_probes2(:,7)==nstate & all_eyet_probes2(:,4)==task,end);
        if task==1
            simpleBarPlot(nstate+0.2*(task*2-3),tempplot,[Colors(nstate,:)],0.35,'k',[],2);
        elseif task==2
            simpleBarPlot(nstate+0.2*(task*2-3),tempplot,[1 1 1;Colors(nstate,:)],0.35,'k',[],2);
        end
    end
end
ylabel('Pupil before')
% ylim([400 1400])
xlim([0.2 3.8])
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
xlim([0.2 3.8])

%% figure Pupil
figure; set(gcf,'position',[4    12   542   762])
xTime=(-20*Fs:5*Fs)/Fs;
mysval={[1],[2],[3]};
for task=1:2
    subplot(2,2,2*(task-1)+1);
    format_fig;
    for nstate=1:length(mysval)
        tempplot=nanmean(all_pup_mat(ismember(all_pup_mat(:,8),mysval{nstate}) & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
        %         tempplot=tempplot-nanmean(tempplot(xTime>0));
        simpleTplot(xTime,tempplot,0,Colors(nstate,:),0,'-',1,0.5,500,0,2);
    end
    legend({'Dist','MW','Interf'})
    
    subplot(2,2,2*(task-1)+2);
    format_fig;
    for nstate=1:length(mysval)
        tempplot=nanmean(all_pup_mat(ismember(all_pup_mat(:,8),mysval{nstate}) & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
        tempplot=tempplot-nanmean(tempplot(xTime>-2 & xTime<0));
        simpleTplot(xTime,tempplot,0,Colors(nstate,:),0,'-',1,0.5,500,0,2);
        xlim([-2 5])
    end
end

figure
format_fig;
for task=1:2
    for nstate=1:length(mysval)
        tempplot=(all_pup_mat(ismember(all_pup_mat(:,8),mysval{nstate}) & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
        tempplot=nanmean(tempplot(:,xTime>-5 & xTime<0),2);
        if task==1
            simpleBarPlot(nstate+0.2*(task*2-3),tempplot,[Colors(nstate,:)],0.35,'k',[],2);
        elseif task==2
            simpleBarPlot(nstate+0.2*(task*2-3),tempplot,[1 1 1;Colors(nstate,:)],0.35,'k',[],2);
        end
    end
end
ylabel('Pupil before')
ylim([400 1400])
xlim([0.2 length(mysval)+0.8])
set(gca,'XTick',1:3,'XTickLabel',{'Dist','MW','Interf'})
xlim([0.2 length(mysval)+0.8])

figure
format_fig;
for task=1:2
    for nstate=1:length(mysval)
        tempplot=(all_pup_mat(ismember(all_pup_mat(:,8),mysval{nstate}) & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
        tempplot=tempplot-nanmean(tempplot(:,xTime>-2 & xTime<0),2);
        tempplot=nanmean(tempplot(:,xTime>2 & xTime<4),2);
        if task==1
            simpleBarPlot(nstate+0.2*(task*2-3),tempplot,[Colors(nstate,:)],0.35,'k',[],2);
        elseif task==2
            simpleBarPlot(nstate+0.2*(task*2-3),tempplot,[1 1 1;Colors(nstate,:)],0.35,'k',[],2);
        end
    end
end
ylabel('Pupil after')
% ylim([400 1400])
xlim([0.2 length(mysval)+0.8])
set(gca,'XTick',1:3,'XTickLabel',{'Dist','MW','Interf'})
xlim([0.2 length(mysval)+0.8])


%% Heatmaps
figure
for task=1:2
    for nstate=1:3
        tempplotX=(all_PosX_mat(all_PosX_mat(:,7)==nstate & all_PosX_mat(:,4)==task,size(temp_behav,2)+1:end));
        tempplotY=(all_PosY_mat(all_PosY_mat(:,7)==nstate & all_PosY_mat(:,4)==task,size(temp_behav,2)+1:end));
        toplot=[];
        for ntr=1:size(tempplotX,1)
            [toplot(:,:,ntr),XEDGES,YEDGES] = histcounts2(tempplotX(ntr,:),tempplotY(ntr,:),[0:10:2500],[0:10:2500]);
            toplot(:,:,ntr)=toplot(:,:,ntr)/nansum(nansum(toplot(:,:,ntr),1),2)*100;
        end
        subplot(2,3,3*(task-1)+nstate);
        format_fig;
        imagesc(XEDGES,YEDGES,nanmean(toplot,3));
        ylim([1000 1500])
        xlim([500 1000])
        heatmaps{task,nstate}=nanmean(toplot,3);
    end
end

figure
for nstate=1:3
    for nstate2=nstate:3
        if nstate2==nstate
            subplot(3,3,3*(nstate-1)+nstate2);
            format_fig;
            imagesc(XEDGES,YEDGES,heatmaps{task,nstate});
            ylim([1000 1500])
            xlim([500 1000])
        else
            subplot(3,3,3*(nstate-1)+nstate2);
            format_fig;
            imagesc(XEDGES,YEDGES,heatmaps{task,nstate}-heatmaps{task,nstate2});
            ylim([1000 1500])
            xlim([500 1000])
        end
        
    end
end

%%
for task=1:2
    figure
    for nstate=1:3
        for nstate2=nstate:3
            if nstate2==nstate
                subplot(3,3,3*(nstate-1)+nstate2);
                format_fig;
                imagesc(XEDGES,YEDGES,heatmaps{task,nstate});
                ylim([1000 1500])
                xlim([500 1000])
            else
                subplot(3,3,3*(nstate-1)+nstate2);
                format_fig;
                imagesc(XEDGES,YEDGES,heatmaps{task,nstate}-heatmaps{task,nstate2});
                ylim([1000 1500])
                xlim([500 1000])
                caxis([-0.6 0.6])
            end
            
        end
    end
end
%% figure Bar
figure;
NamesPlot={'# blinks','dur blinks','# saccades'};
YLims=[5 25;80 140;5 25];
for nplot=1:3
    subplot(1,3,nplot);
    format_fig;
    for task=1:2
        for nstate=1:3
            tempplot=(all_eyet_probes(all_eyet_probes(:,7)==nstate & all_eyet_probes(:,4)==task,size(temp_behav,2)+1:end));
            if task==1
                simpleBarPlot(nstate+0.2*(task*2-3),tempplot(:,nplot),[Colors(nstate,:)],0.35,'k',[],2);
            elseif task==2
                simpleBarPlot(nstate+0.2*(task*2-3),tempplot(:,nplot),[1 1 1;Colors(nstate,:)],0.35,'k',[],2);
            end
        end
    end
    title(NamesPlot{nplot})
    ylim(YLims(nplot,:))
    xlim([0.2 3.8])
end


%% Figure behaviour (move elswhere)
% % Look
% figure; set(gcf,'Position',[1     8   282   797])
% varNames={'CorrGo', 'CorrNoGo', 'RTgo',  'd''', 'crit'};
% for nvar=1:5
%     subplot(5,1,nvar);
%     simpleCorPlotsetbin(all_eyet_probes(:,6),all_eyet_probes(:,13+nvar),1:2,{'o','b','b',96},[],0,1);
%     if nvar==5
%         xlabel('Look?')
%     set(gca,'xtick',1:2,'xticklabel',{'Yes','No'})
%     else
%     set(gca,'xtick',[])
%     end
%     ylabel(varNames{nvar})
%     title('')
% end

% Awareness
figure; set(gcf,'Position',[1     8   230   700])
varNames={'CorrGo', 'CorrNoGo',  'd''', 'crit', 'RTgo'};
varCol=[14 15 17 18 16];
for nvar=1:5
    subplot(5,1,nvar);
    simpleCorPlotsetbin(5-all_eyet_probes(:,9),all_eyet_probes(:,varCol(nvar)),1:4,{'o','b','b',96,20},[],0,1);
    if nvar==5
        xlabel('Awareness?')
        set(gca,'xtick',[1 4],'xticklabel',{'Low','High'})
    else
        set(gca,'xtick',[])
    end
    %     ylabel(varNames{nvar})
    title('')
end

% Intentionality
figure; set(gcf,'Position',[1     8   230   700])
varNames={'CorrGo', 'CorrNoGo',  'd''', 'crit', 'RTgo'};
varCol=[14 15 17 18 16];
for nvar=1:5
    subplot(5,1,nvar);
    simpleCorPlotsetbin(5-all_eyet_probes(:,10),all_eyet_probes(:,varCol(nvar)),1:4,{'o','b','b',96,20},[],0,1);
    if nvar==5
        xlabel('Intentionality?')
        set(gca,'xtick',[1 4],'xticklabel',{'Low','High'})
    else
        set(gca,'xtick',[])
    end
    %     ylabel(varNames{nvar})
    title('')
end
% Engagement
figure; set(gcf,'Position',[1     8   230   700])
varNames={'CorrGo', 'CorrNoGo',  'd''', 'crit', 'RTgo'};
varCol=[14 15 17 18 16];
for nvar=1:5
    subplot(5,1,nvar);
    simpleCorPlotsetbin(all_eyet_probes(:,11),all_eyet_probes(:,varCol(nvar)),1:4,{'o','b','b',96,20},[],0,1);
    if nvar==5
        xlabel('Engagement?')
        set(gca,'xtick',[1 4],'xticklabel',{'Low','High'})
    else
        set(gca,'xtick',[])
    end
    %     ylabel(varNames{nvar})
    title('')
end
% Performance
figure; set(gcf,'Position',[1     8   230   700])
varNames={'CorrGo', 'CorrNoGo',  'd''', 'crit', 'RTgo'};
varCol=[14 15 17 18 16];
for nvar=1:5
    subplot(5,1,nvar);
    simpleCorPlotsetbin(all_eyet_probes(:,12),all_eyet_probes(:,varCol(nvar)),1:4,{'o','b','b',96,20},[],0,1);
    if nvar==5
        xlabel('Performance?')
        set(gca,'xtick',[1 4],'xticklabel',{'Low','High'})
    else
        set(gca,'xtick',[])
    end
    %     ylabel(varNames{nvar})
    title('')
end
% Alterness
figure; set(gcf,'Position',[1     8   230   700])
varNames={'CorrGo', 'CorrNoGo',  'd''', 'crit', 'RTgo'};
varCol=[14 15 17 18 16];
for nvar=1:5
    subplot(5,1,nvar);
    simpleCorPlotsetbin(all_eyet_probes(:,13),all_eyet_probes(:,varCol(nvar)),1:4,{'o','b','b',96,20},[],0,1);
    if nvar==5
        xlabel('Sleepiness?')
        set(gca,'xtick',[1 4],'xticklabel',{'Low','High'})
    else
        set(gca,'xtick',[])
    end
    %     ylabel(varNames{nvar})
    title('')
    xlim([0.5 4.5])
end    % 14: resp Q8 Alterness (1 (alert) to 4 (sleepy))

%%
% # blinks / dur blinks / # saccades / amp saccades / vel saccades
figure; set(gcf,'Position',[1     8   230   700])
varNames={'# blinks', '# sacc',  'amp sacc', 'vel sacc'};
varCol=[20 21 22 23];
for nvar=1:4
    subplot(4,1,nvar);
    simpleCorPlotsetbin(5-all_eyet_probes(:,9),all_eyet_probes(:,varCol(nvar)),1:4,{'o','b','b',72,4},[],0,1);
    if nvar==4
        xlabel('Awareness?')
        set(gca,'xtick',[1 4],'xticklabel',{'Low','High'})
    else
        set(gca,'xtick',[])
    end
    %     ylabel(varNames{nvar})
    title('')
    xlim([0.5 4.5])
end    % 14: resp Q8 Alterness (1 (alert) to 4 (sleepy))


figure; set(gcf,'Position',[1     8   230   700])
varNames={'# blinks', '# sacc',  'amp sacc', 'vel sacc'};
varCol=[20 21 22 23];
for nvar=1:4
    subplot(4,1,nvar);
    simpleCorPlotsetbin(all_eyet_probes(:,11),all_eyet_probes(:,varCol(nvar)),1:4,{'o','b','b',72,4},[],0,1);
    if nvar==4
        xlabel('Engagement?')
        set(gca,'xtick',[1 4],'xticklabel',{'Low','High'})
    else
        set(gca,'xtick',[])
    end
    %     ylabel(varNames{nvar})
    title('')
    xlim([0.5 4.5])
end    % 14: resp Q8 Alterness (1 (alert) to 4 (sleepy))


figure; set(gcf,'Position',[1     8   230   700])
varNames={'# blinks', '# sacc',  'amp sacc', 'vel sacc'};
varCol=[20 21 22 23];
for nvar=1:4
    subplot(4,1,nvar);
    simpleCorPlotsetbin(all_eyet_probes(:,12),all_eyet_probes(:,varCol(nvar)),1:4,{'o','b','b',72,4},[],0,1);
    if nvar==4
        xlabel('Performance?')
        set(gca,'xtick',[1 4],'xticklabel',{'Low','High'})
    else
        set(gca,'xtick',[])
    end
    %     ylabel(varNames{nvar})
    title('')
    xlim([0.5 4.5])
end    % 14: resp Q8 Alterness (1 (alert) to 4 (sleepy))


figure; set(gcf,'Position',[1     8   230   700])
varNames={'# blinks', '# sacc',  'amp sacc', 'vel sacc'};
varCol=[20 21 22 23];
for nvar=1:4
    subplot(4,1,nvar);
    simpleCorPlotsetbin(all_eyet_probes(:,13),all_eyet_probes(:,varCol(nvar)),1:4,{'o','b','b',72,4},[],0,1);
    if nvar==4
        xlabel('Sleepiness?')
        set(gca,'xtick',[1 4],'xticklabel',{'Low','High'})
    else
        set(gca,'xtick',[])
    end
    %     ylabel(varNames{nvar})
    title('')
    xlim([0.5 4.5])
end    % 14: resp Q8 Alterness (1 (alert) to 4 (sleepy))

%% Pupil and taks parameters
myCols=[9 11 12 13];
myColsNames={'Awa','Eng','Perf','Sleep'};
myTaskNames={'F','D'};
figure;
for nvar=1:2
    for task=1:2
        subplot(2,2,2*(task-1)+nvar);  format_fig;
        for nstate=3:-1:1
            tempplot=(all_pup_mat(all_pup_mat(:,7)==nstate & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
            tempplot=nanmean(tempplot(:,xTime>-20 & xTime<0),2);
            tempplot2=all_eyet_probes(all_eyet_probes(:,7)==nstate & all_eyet_probes(:,4)==task,myCols(nvar));
            simpleCorPlotsetbin(tempplot2,tempplot,1:4,{'o',Colors(nstate,:),Colors(nstate,:),72,4},[],0,1);
        end
        title(sprintf('%s - %s',myTaskNames{task},myColsNames{nvar}))
        if task==1
            ylim([500 1500])
        elseif task==2
            ylim([800 2000]-200)
        end
        xlim([0.5 4.5])
        set(gca,'xtick',[1 4],'xticklabel',{'Low','High'})
    end
end


myCols=[9 11 12 13];
myColsNames={'Awareness','Engagement','Performance','Sleepiness'};
myTaskNames={'F','D'};
figure; set(gcf,'Position',[  1         325        1392         480])
for nvar=1:4
    subplot(1,4,nvar);  format_fig;
    for nstate=2:-1:1
        tempplot=(all_pup_mat(all_pup_mat(:,7)==nstate,size(temp_behav,2)+1:end));
        tempplot=nanmean(tempplot(:,xTime>-20 & xTime<0),2);
        tempplot2=all_eyet_probes(all_eyet_probes(:,7)==nstate,myCols(nvar));
        simpleCorPlotsetbin(tempplot2,tempplot,1:4,{'o',Colors(nstate,:),Colors(nstate,:),72,4},[],0,1);
    end
    title(sprintf('%s',myColsNames{nvar}))
    if task==1
        ylim([500 1500])
    elseif task==2
        ylim([800 2000]-200)
    end
    if nvar==1
        ylabel('Pupil size')
    end
    xlim([0.5 4.5])
    set(gca,'xtick',[1 4],'xticklabel',{'Low','High'})
end

%%
figure
format_fig;
% for task=1:2
for nstate=1:3
    if nstate~=2
        tempplot=(all_pup_mat(all_pup_mat(:,7)==nstate & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
    else
        tempplot=(all_pup_mat(all_pup_mat(:,7)==nstate & all_pup_mat(:,8)==nstate & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
    end
    tempplot=nanmean(tempplot(:,xTime>-5 & xTime<0),2);
    %         if task==1
    simpleBarPlot(nstate,tempplot,[Colors(nstate,:)],0.9,'k',[],2);
    %         elseif task==2
    %             simpleBarPlot(nstate+0.2*(task*2-3),tempplot,[1 1 1;Colors(nstate,:)],0.35,'k',[],2);
    %         end
end
% end
ylabel('Pupil before')
ylim([400 1400])
xlim([0.2 3.8])
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
xlim([0.2 3.8])

%%
% n nbl npr these_probes(npr,5) this_pr_tridx probe_details temp_corr_go
% temp_corr_nogo temp_rt_go  dprime crit sum(temp_blinks(:,end)<400) mean(temp_blinks(temp_blinks(:,end)<400,end)) size(temp_sacc,1) mean(ampSac) mean(temp_sacc(:,end))
table=array2table(all_eyet_probes2,'VariableNames',{'SubID','nBlock','nProbe','Task','ProbeIdx','Look','State','Ori','Awa','Int','Eng','Perf','Alert',...
    'CorrGo','CorrNoGo','RTgo','dp','crit','nBlinks','durBlinks','nSacc','ampSacc','velSacc','pupBef','pupAfter'});
table(table.State==4,:)=[];
table.SubID=categorical(table.SubID);
table.Task=categorical(table.Task);
table.State=categorical(table.State);
mdl1=fitlme(table,'nBlinks~nBlock*Task*State+(1|SubID)');
mdl2=fitlme(table,'durBlinks~nBlock*Task*State+(1|SubID)');
mdl3=fitlme(table,'nSacc~nBlock*Task*State+(1|SubID)');
mdl4=fitlme(table,'ampSacc~nBlock*Task*State+(1|SubID)');
mdl5=fitlme(table,'velSacc~nBlock*Task*State+(1|SubID)');
mdl6=fitlme(table,'pupBef~nBlock*Task*State+(1|SubID)');
mdl7=fitlme(table,'pupAfter~nBlock*Task*State+(1|SubID)');