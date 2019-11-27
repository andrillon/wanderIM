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
            
%             thesegotrials=find(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx & these_trials(:,5)~=3);
%             thesenogotrials=find(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx & these_trials(:,5)==3);
%             temp_testres=these_trials([thesegotrials(end-17:end) ; thesenogotrials(end-1:end)],:);
            
            temp_testresgo=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,12)),:);
            tcorr_go=(temp_testresgo(end-17:end,12))';%/corr_go(n,these_probes(npr,5));
            temp_testresnogo=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,11)),:);
            tcorr_nogo=(temp_testresnogo(end-1:end,11))';%/corr_go(n,these_probes(npr,5));
            
            temp_corr_go=nanmean(tcorr_go);
            temp_corr_nogo=nanmean(tcorr_nogo);
            temp_rt_go=nan;
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
            temp_Pupil=EL_data.clean_pupilSize((-20*Fs:5*Fs)+idx_this_probe);
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
    subplot(2,1,task);
    format_fig;
    for nstate=1:3
        tempplot=(all_pup_mat(all_pup_mat(:,7)==nstate & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
                tempplot=tempplot-repmat(nanmean(tempplot(:,xTime>0.5 & xTime<2.5),2),[1 size(tempplot,2)]);
        simpleTplot(xTime,tempplot,0,Colors(nstate,:),0,'-',0.5,1,[],0,2);
    end
    ylim([-250 120])
%     legend({'ON','MW','MB'})
    
%     subplot(2,2,2*(task-1)+2);
%     format_fig;
%     for nstate=1:3
%         tempplot=nanmean(all_pup_mat(all_pup_mat(:,7)==nstate & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
%         tempplot=tempplot-nanmean(tempplot(xTime>-2 & xTime<0));
%         simpleTplot(xTime,tempplot,0,Colors(nstate,:),0,'-',1,0.5,500,0,2);
%         xlim([-2 5])
%     end
%     legend({'ON','MW','MB'})
end

%%
myS=unique(all_pup_mat(:,1));

figure; set(gcf,'position',[680   696   427   282])
format_fig;
for task=1:2
    for nstate=1:3
        tempplot=squeeze(all_pup_mat(all_pup_mat(:,7)==nstate & all_pup_mat(:,4)==task,size(temp_behav,2)+1:end));
        tempplot=nanmean(tempplot(:,xTime>-5 & xTime<0),2)-nanmean(tempplot(:,xTime>0.5 & xTime<2.5),2);
        tempS=squeeze(all_pup_mat(all_pup_mat(:,7)==nstate & all_pup_mat(:,4)==task,1));
        tempbyS=[];
        for nS=1:length(myS)
            tempbyS(nS,1)=nanmean(tempplot(tempS==myS(nS)));
            tempbyS(nS,2)=nansum(tempS==myS(nS));
        end
         mean_temp(nstate)=nansum(tempbyS(:,1).*tempbyS(:,2))./sum(tempbyS(:,2));
        sem_temp(nstate)=std(tempbyS(~isnan(tempbyS(:,1)),1),tempbyS(~isnan(tempbyS(:,1)),2))/sqrt(sum(~isnan(tempbyS(:,1)))-1);
       
       line(nstate*[1 1]+(2*task-3)*0.1,[-1 1]*sem_temp(nstate)+mean_temp(nstate),'Color',Colors(nstate,:),'LineWidth',3)
    end
    hold on;
    plot((1:3)+(2*task-3)*0.1,mean_temp,'Color',[1 1 1]*0.5,'LineWidth',2)
    for nstate=1:3
        if task==1
        scatter(nstate+(2*task-3)*0.1,mean_temp(nstate),'MarkerFaceColor',Colors(nstate,:),'MarkerEdgeColor',Colors(nstate,:),'SizeData',144,'LineWidth',3);
        else
        scatter(nstate+(2*task-3)*0.1,mean_temp(nstate),'MarkerFaceColor','w','MarkerEdgeColor',Colors(nstate,:),'SizeData',144,'LineWidth',3);
        end
    end
end
ylabel('Pupil before-after')
% ylim([400 1400])
xlim([0.2 3.8])
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
xlim([0.2 3.8])

tempplot=(all_pup_mat(:,size(temp_behav,2)+1:end));
tempplot=nanmean(tempplot(:,xTime>-5 & xTime<0),2)-nanmean(tempplot(:,xTime>0.5 & xTime<2.5),2);
fortbl=[all_pup_mat(:,1:size(temp_behav,2)) tempplot];
%             temp_behav=[n nbl npr these_probes(npr,5) this_pr_tridx probe_details temp_corr_go temp_corr_nogo temp_rt_go  dprime crit];
fortbl=fortbl(fortbl(:,7)<4,:);
tbl=array2table(fortbl,'VariableNames',{'SubID','nBl','nPr','Task','Pr','look','MS','origin','awareness','intention','engagement','performance','alterness','GO','NOGO','RT','dp','crit','PUP'});
tbl.MS=categorical(tbl.MS);
tbl.Task=categorical(tbl.Task);
tbl.MS=reordercats(tbl.MS,[2 1 3]);
mdl_0= fitlme(tbl,'PUP~1+(1|SubID)');
mdl_1= fitlme(tbl,'PUP~Task+(1|SubID)');
mdl_2= fitlme(tbl,'PUP~Task+MS+(1|SubID)');
mdl_3= fitlme(tbl,'PUP~Task*MS+(1|SubID)');

