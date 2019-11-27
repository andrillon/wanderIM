%%
clear all
% close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))

% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path_adhd filesep 'behav/'];
data_path_ori=[root_path filesep 'behav/'];
files=dir([data_path filesep 'MWADHD_behavres_s4*.mat']);
files_ori=dir([data_path_ori filesep 'wanderIM_behavres_s3*.mat']);

state_colours=[0 146 146 ; 182 109 255; 219 209 0; 200 200 200]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];

GroupStates={[1],[2],[3,4]};
%load([data_path filesep 'CARS_quest'])
%%
all_test_res_ctrl=[];
for n=1:20 %length(files_ori)
    % load
    load([data_path_ori filesep files_ori(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s duration of task %3.1f min\n',SubID,(all_GrandEnd-all_GrandStart)/60)
age_ctr(n)=str2num(SubjectInfo.Age);
gender_ctr(n)=strcmp(SubjectInfo.Gender,'F');
   %%% cleaning too fast RTs
    RTs=[test_res(:,10)-test_res(:,8)];
    %     findFalseStarts=find(test_res(:,10)-test_res(:,8)<0.1);
    %     test_res(findFalseStarts,[10 11 12])=NaN;
    %     warning('correcting for false starts')
    nbt=1;
    tp_nogos=test_res(~isnan(test_res(:,11)),11);
    tp_gos=test_res(~isnan(test_res(:,12)),12);
    [dprime_test_ctrl(n,nbt), crit_test_ctrl(n,nbt)]=calc_dprime((tp_gos==1),(tp_nogos==0));
    corr_go_ctrl(n,nbt)=nanmean(tp_gos);
    corr_nogo_ctrl(n,nbt)=nanmean(tp_nogos);
    
    rt_gos_ctrl(n,nbt)=nanmedian(RTs(~isnan(test_res(:,12))));
    rt_nogos_ctrl(n,nbt)=nanmedian(RTs(~isnan(test_res(:,11))));
    all_test_res_ctrl=[all_test_res_ctrl ; [n*ones(size(test_res,1),1) test_res(:,[1 2 4 11 12])]];
        temp=RTs(~isnan(test_res(:,12)));
    rt_distrib_go_ctrl(n,:)=histc(temp',0:0.05:1);
    temp=RTs(~isnan(test_res(:,11)));
    rt_distrib_nogo_ctrl(n,:)=histc(temp',0:0.05:1);
    
    temp=probe_res(:,[5 31:38]);
    clear pption pption2
    for nstate=1:length(GroupStates)
        pption(nstate)=mean(ismember(temp(:,3),GroupStates{nstate}));
        pption2(nstate)=mean(temp(temp(:,3)==2,4)==nstate);
    end
    pption_MS_ctrl(n,:)=pption;
    pption_Ori_ctrl(n,:)=pption2;
    
        tp_probes=probe_res(:,31:38);
    for nstate=1:length(GroupStates)
        
        if sum(tp_probes(:,2)==(nstate))~=0
            
            mean_byprobe_awa_ctrl(n,nstate)=5-nanmean(tp_probes(ismember(tp_probes(:,2),GroupStates{nstate}),4));
            mean_byprobe_wil_ctrl(n,nstate)=5-nanmean(tp_probes(ismember(tp_probes(:,2),GroupStates{nstate}),5));
            mean_byprobe_eng_ctrl(n,nstate)=nanmean(tp_probes(ismember(tp_probes(:,2),GroupStates{nstate}),6));
            mean_byprobe_prf_ctrl(n,nstate)=nanmean(tp_probes(ismember(tp_probes(:,2),GroupStates{nstate}),7));
            mean_byprobe_vig_ctrl(n,nstate)=5-nanmean(tp_probes(ismember(tp_probes(:,2),GroupStates{nstate}),8));
        else
            mean_byprobe_awa_ctrl(n,nstate)=nan;
            mean_byprobe_wil_ctrl(n,nstate)=nan;
            mean_byprobe_eng_ctrl(n,nstate)=nan;
            mean_byprobe_prf_ctrl(n,nstate)=nan;
            mean_byprobe_vig_ctrl(n,nstate)=nan;
        end
    end
    mean_allprobes_awa_ctrl(n,nstate)=5-nanmean(tp_probes(:,4));
    mean_allprobes_wil_ctrl(n,nstate)=5-nanmean(tp_probes(:,5));
    mean_allprobes_eng_ctrl(n,nstate)=nanmean(tp_probes(:,6));
    mean_allprobes_prf_ctrl(n,nstate)=nanmean(tp_probes(:,7));
    mean_allprobes_vig_ctrl(n,nstate)=5-nanmean(tp_probes(:,8));
end


%
all_test_resprobes_perblock=[];
all_test_res=[];
all_probes_timing=[];
all_probes_mat=[];
RTs=[];
all_probes_mat_fullgo=[];
all_probes_mat_fullnogo=[];
TreatmentSuspension=[401	1; 402	1; 403	1; 404	1; 405	1; 406	1; 407	0; 408	0];
nc=0;
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
    if TreatmentSuspension(TreatmentSuspension(:,1)==str2num(SubID),2)==0
    fprintf('... %s DISCARDED\n',SubID)
        continue;
    end
    nc=nc+1;
    age_adhd(nc)=str2num(SubjectInfo.Age);
gender_adhd(nc)=strcmp(SubjectInfo.Gender,'F');
%    CARS_flag(n)=CARS_bool(CARS_bool(:,1)==str2num(SubID),2);
    % SART
    %  1: num block
    %  2: block cond (1: Faces / 2: Squares)
    %  3: image set
    %  4: num trial
    %  5: seq trial
    %  6: target
    %  7: resp
    %  8: stim onset
    %  9: stim pre
    % 10: resp onset
    % 11: nogo
    % 12: go
    
    %%% cleaning too fast RTs
    RTs=[test_res(:,10)-test_res(:,8)];
    %     findFalseStarts=find(test_res(:,10)-test_res(:,8)<0.1);
    %     test_res(findFalseStarts,[10 11 12])=NaN;
    %     warning('correcting for false starts')
    nbt=1;
    tp_nogos=test_res(~isnan(test_res(:,11)),11);
    tp_gos=test_res(~isnan(test_res(:,12)),12);
    [dprime_test(nc,nbt), crit_test(nc,nbt)]=calc_dprime((tp_gos==1),(tp_nogos==0));
    corr_go(nc,nbt)=nanmean(tp_gos);
    corr_nogo(nc,nbt)=nanmean(tp_nogos);
    
    rt_gos(nc,nbt)=nanmedian(RTs(~isnan(test_res(:,12))));
    rt_nogos(nc,nbt)=nanmedian(RTs(~isnan(test_res(:,11))));
    all_test_res=[all_test_res ; [n*ones(size(test_res,1),1) test_res(:,[1 2 4 11 12])]];
    
    temp=RTs(~isnan(test_res(:,12)));
    rt_distrib_go(nc,:)=histc(temp',0:0.05:1);
    temp=RTs(~isnan(test_res(:,11)));
    rt_distrib_nogo(nc,:)=histc(temp',0:0.05:1);
    % probe
    %  1: num probe within block
    %  2: probe onset (expected)
    %  3: probe onset (actual)
    %  4: num block
    %  5: block type (1: FACE / 2: SQUARE)
    %  6: num trial
    %  7: resp Q1 Looking (1: Yes / 2: No)
    %  8: resp Q2 Mind-State (1: On / 2: MW / 3: MB / 4: ?)
    %  9: resp Q3 Origine (1: room / 2: personal / 3: task)
    % 10: resp Q4 Awareness (1 (fully) to 4 (not at all))
    % 11: resp Q5 Intention (1 (fully) to 4 (not at all))
    % 12: resp Q6 Engagement (1 (not) to 4 (very))
    % 13: resp Q7 Performance (1 (bad) to 4 (good))
    % 14: resp Q8 Alterness (1 (alert) to 4 (sleepy))
    % 15-22: onset R1-8
    % 23-30: onset Q1-8
    % 31-38: resp number Q1-8
    temp=probe_res(:,[5 31:38]);
    clear pption pption2
    for nstate=1:length(GroupStates)
        pption(nstate)=mean(ismember(temp(:,3),GroupStates{nstate}));
        pption2(nstate)=mean(temp(temp(:,3)==2,4)==nstate);
    end
    pption_MS(nc,:)=pption;
    pption_Ori(nc,:)=pption2;
    for nb=1:max(probe_res(:,4))
        temp=probe_res(probe_res(:,4)==nb,[5 31:38]);
        temp2=test_res(test_res(:,1)==nb,[11 12]);
        all_test_resprobes_perblock=[all_test_resprobes_perblock ; [n*ones(size(temp,1),1) nb*ones(size(temp,1),1) temp]];
        all_probes_timing=[all_probes_timing ; diff(probe_res(probe_res(:,3)==nb,2))];
    end
    tp_probes=probe_res(:,31:38);
    mean_bytask_prf(nc)=nanmean(tp_probes(:,7));
    
    
    for nstate=1:length(GroupStates)
        mean_count_mindstates(nc,nstate)=sum(ismember(tp_probes(:,2),GroupStates{nstate}));
        mean_count_ori(nc,nstate)=sum(tp_probes(:,3)==(nstate));
        
        if sum(ismember(tp_probes(:,2),GroupStates{nstate}))~=0
            
            mean_byprobe_awa(nc,nstate)=5-nanmean(tp_probes(ismember(tp_probes(:,2),GroupStates{nstate}),4));
            mean_byprobe_wil(nc,nstate)=5-nanmean(tp_probes(ismember(tp_probes(:,2),GroupStates{nstate}),5));
            mean_byprobe_eng(nc,nstate)=nanmean(tp_probes(ismember(tp_probes(:,2),GroupStates{nstate}),6));
            mean_byprobe_prf(nc,nstate)=nanmean(tp_probes(ismember(tp_probes(:,2),GroupStates{nstate}),7));
            mean_byprobe_vig(nc,nstate)=5-nanmean(tp_probes(ismember(tp_probes(:,2),GroupStates{nstate}),8));
        else
            mean_byprobe_awa(nc,nstate)=nan;
            mean_byprobe_wil(nc,nstate)=nan;
            mean_byprobe_eng(nc,nstate)=nan;
            mean_byprobe_prf(nc,nstate)=nan;
            mean_byprobe_vig(nc,nstate)=nan;
        end
    end
    mean_allprobes_awa(nc,nstate)=5-nanmean(tp_probes(:,4));
    mean_allprobes_wil(nc,nstate)=5-nanmean(tp_probes(:,5));
    mean_allprobes_eng(nc,nstate)=nanmean(tp_probes(:,6));
    mean_allprobes_prf(nc,nstate)=nanmean(tp_probes(:,7));
    mean_allprobes_vig(nc,nstate)=5-nanmean(tp_probes(:,8));
    
    %%% performance between probes
    for nbl=1:max(probe_res(:,4))
        these_probes=probe_res(probe_res(:,4)==nbl,:);
        these_trials=test_res(test_res(:,1)==nbl,:);
        for npr=1:max(these_probes(:,1))
            this_pr_tridx=these_probes(npr,6);
            %             if npr==1
            %                 last_pr_tridx=0;
            %             else
            %                 last_pr_tridx=these_probes(npr-1,6);
            %             end
            last_pr_tridx=this_pr_tridx-21;
            number_STD(nc,nbl,npr)=sum(these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,5)~=3);
            number_DEV(nc,nbl,npr)=sum(these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,5)==3);
            
            temp_testres=these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,:);
            %             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            tcorr_go=nanmean(temp_testres(:,12));%/corr_go(n,these_probes(npr,5));
            tcorr_nogo=nanmean(temp_testres(:,11));%/corr_nogo(n,these_probes(npr,5));
            num_go=sum(~isnan(temp_testres(:,12)));
            num_nogo=sum(~isnan(temp_testres(:,11)));
            rt_go=nanmean(temp_testres(~isnan(temp_testres(:,12)),10)-temp_testres(~isnan(temp_testres(:,12)),8));
            rt_nogo=nanmean(temp_testres(~isnan(temp_testres(:,11)),10)-temp_testres(~isnan(temp_testres(:,11)),8));
            

            
            tdp=calc_dprime((temp_testres(~isnan(temp_testres(:,12)),12)==1),(temp_testres(~isnan(temp_testres(:,11)),11)==0));
            all_probes_mat=[all_probes_mat ; [n nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_go tcorr_nogo num_go num_nogo tdp (these_probes(npr,37)) rt_go rt_nogo]];
            
            temp_testres=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,12)),:);
            tcorr_go=(temp_testres(end-19:end,12))';%/corr_go(n,these_probes(npr,5));
             temp_testres=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,11)),:);
            tcorr_nogo=(temp_testres(end-1:end,11))';%/corr_go(n,these_probes(npr,5));
            all_probes_mat_fullgo=[all_probes_mat_fullgo ; [n nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_go]];
            all_probes_mat_fullnogo=[all_probes_mat_fullnogo ; [n nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_nogo]];
        end
    end
end



% 
%% Behaviour - Performance on SART
figure; set(gcf,'Position',[529   943   615   750]);
subplot(3,2,1); format_fig;
simple_violinplot(dprime_test_ctrl,[],1,'ViolinColor',[0.7 0 0.1],'Width',0.45,'BoxWidth',0.04);
simple_violinplot(dprime_test,[],2,'ViolinColor',[0.1 0 0.7],'Width',0.45,'BoxWidth',0.04);
set(gca,'XTick',1:2,'XTickLabel',{'CTRL','ADHD'})
ylabel('d-prime')
xlim([0.2 2.8])
ylim([0 2.5])
[pV,h,stats]=ranksum(dprime_test_ctrl,dprime_test);
fprintf('... ADHD vs Control ... dprime: pv=%g U=%g\n',pV,stats.ranksum)

subplot(3,2,2); format_fig;
simple_violinplot(crit_test_ctrl,[],1,'ViolinColor',[0.7 0 0.1],'Width',0.45,'BoxWidth',0.04);
simple_violinplot(crit_test,[],2,'ViolinColor',[0.1 0 0.7],'Width',0.45,'BoxWidth',0.04);
set(gca,'XTick',1:2,'XTickLabel',{'CTRL','ADHD'})
ylabel('criterion')
xlim([0.2 2.8])
[pV,h,stats]=ranksum(crit_test_ctrl,crit_test);
fprintf('... ADHD vs Control ... crit: pv=%g U=%g\n',pV,stats.ranksum)

subplot(3,2,3); format_fig;
simple_violinplot(100*corr_go_ctrl,[],1,'ViolinColor',[0.7 0 0.1],'Width',0.45,'BoxWidth',0.04);
simple_violinplot(100*corr_go,[],2,'ViolinColor',[0.1 0 0.7],'Width',0.45,'BoxWidth',0.04);
set(gca,'XTick',1:2,'XTickLabel',{'CTRL','ADHD'})
ylabel('Perf. GO')
xlim([0.2 2.8])
ylim([90 100])
[pV,h,stats]=ranksum(corr_go_ctrl,corr_go);
fprintf('... ADHD vs Control ... corr go: pv=%g U=%g\n',pV,stats.ranksum)

subplot(3,2,4); format_fig;
simple_violinplot(100*corr_nogo_ctrl,[],1,'ViolinColor',[0.7 0 0.1],'Width',0.45,'BoxWidth',0.04);
simple_violinplot(100*corr_nogo,[],2,'ViolinColor',[0.1 0 0.7],'Width',0.45,'BoxWidth',0.04);
set(gca,'XTick',1:2,'XTickLabel',{'CTRL','ADHD'})
ylabel('Perf. NO-GO')
xlim([0.2 2.8])
ylim([0 100])
[pV,h,stats]=ranksum(corr_nogo_ctrl,corr_nogo);
fprintf('... ADHD vs Control ... corr nogo: pv=%g U=%g\n',pV,stats.ranksum)

subplot(3,2,5); format_fig;
simple_violinplot(rt_gos_ctrl,[],1,'ViolinColor',[0.7 0 0.1],'Width',0.45,'BoxWidth',0.04);
simple_violinplot(rt_gos,[],2,'ViolinColor',[0.1 0 0.7],'Width',0.45,'BoxWidth',0.04);
set(gca,'XTick',1:2,'XTickLabel',{'CTRL','ADHD'})
ylabel('RT. GO')
xlim([0.2 2.8])
ylim([0.4 0.6])
[pV,h,stats]=ranksum(rt_gos_ctrl,rt_gos);
fprintf('... ADHD vs Control ... RT GO: pv=%g U=%g\n',pV,stats.ranksum)

subplot(3,2,6); format_fig;
simple_violinplot(rt_nogos_ctrl,[],1,'ViolinColor',[0.7 0 0.1],'Width',0.45,'BoxWidth',0.04);
simple_violinplot(rt_nogos,[],2,'ViolinColor',[0.1 0 0.7],'Width',0.45,'BoxWidth',0.04);
set(gca,'XTick',1:2,'XTickLabel',{'CTRL','ADHD'})
ylabel('RT. NO-GO')
xlim([0.2 2.8])
ylim([0.4 0.6])
[pV,h,stats]=ranksum(rt_nogos_ctrl,rt_nogos);
fprintf('... ADHD vs Control ... RT NOGO: pv=%g U=%g\n',pV,stats.ranksum)

%% Probes - % mind-state
pption_MS2=pption_MS;
pption_MS_ctrl2=pption_MS_ctrl;
% pption_MS2(:,3)=pption_MS(:,3)+pption_MS(:,4);
% pption_MS_ctrl2(:,3)=pption_MS_ctrl(:,3)+pption_MS_ctrl(:,4);
figure;  set(gcf,'Position',[-29        1275        1218         412]);
subplot(1,2,1); format_fig;
for nst=1:size(pption_MS_ctrl,2)
    hb(1)=simpleBarPlot(nst-0.2,100*squeeze(pption_MS_ctrl(:,nst)),[1 1 1;state_colours(nst,:)],0.35,'k',[],4);
    hb(2)=simpleBarPlot(nst+0.2,100*squeeze(pption_MS(:,nst)),[state_colours(nst,:)],0.35,'k',{2 100*squeeze(pption_MS_ctrl2(:,nst)) 0.05},4);
end
set(gca,'XTick',1:3,'XTickLabel',{'On','MW','MB'})
ylabel('% probes')
title('Mind-State')
xlim([0.2 size(pption_MS_ctrl,2)+0.8])
legend(hb,{'Control','ADHD'})

subplot(1,2,2); format_fig;
for nst=1:3
    simpleBarPlot(nst-0.2,100*squeeze(pption_Ori_ctrl(:,nst)),[1 1 1;state_colours(nst,:)],0.35,'k',[],4);
    simpleBarPlot(nst+0.2,100*squeeze(pption_Ori(:,nst)),[state_colours(nst,:)],0.35,'k',[],4);
end
set(gca,'XTick',1:3,'XTickLabel',{'room','perso','task'})
ylabel('% probes')
title('Origin of MW')
xlim([0.2 3.8])

%%
figure; format_fig;
for nsta=1:length(GroupStates)
    subplot(2,3,1)
     tempctrl=squeeze(mean(mean_byprobe_awa_ctrl(:,nsta),2));
   temp=squeeze(mean(mean_byprobe_awa(:,nsta),2));
    simpleBarPlot(nsta-0.2,tempctrl,[state_colours(nsta,:)],0.35,'k',[],4);
    simpleBarPlot(nsta+0.2,temp,[1 1 1;state_colours(nsta,:)],0.35,'k',{2 tempctrl 0.05},4);
    xlim([0.2 length(GroupStates)+0.8])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title('Awareness')
    
    subplot(2,3,2)
    tempctrl=squeeze(mean(mean_byprobe_wil_ctrl(:,nsta),2));
    temp=squeeze(mean(mean_byprobe_wil(:,nsta),2));
    simpleBarPlot(nsta-0.2,tempctrl,[state_colours(nsta,:)],0.35,'k',[],4);
    simpleBarPlot(nsta+0.2,temp,[1 1 1;state_colours(nsta,:)],0.35,'k',{2 tempctrl 0.05},4);
    xlim([0.2 length(GroupStates)+0.8])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title('Intention')
    
    subplot(2,3,3)
    tempctrl=squeeze(mean(mean_byprobe_eng_ctrl(:,nsta),2));
    temp=squeeze(mean(mean_byprobe_eng(:,nsta),2));
    simpleBarPlot(nsta-0.2,tempctrl,[state_colours(nsta,:)],0.35,'k',[],4);
    simpleBarPlot(nsta+0.2,temp,[1 1 1;state_colours(nsta,:)],0.35,'k',{2 tempctrl 0.05},4);
    xlim([0.2 length(GroupStates)+0.8])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title('Engagement')
    
    subplot(2,3,4)
    tempctrl=squeeze(mean(mean_byprobe_prf_ctrl(:,nsta),2));
    temp=squeeze(mean(mean_byprobe_prf(:,nsta),2));
    simpleBarPlot(nsta-0.2,tempctrl,[state_colours(nsta,:)],0.35,'k',[],4);
    simpleBarPlot(nsta+0.2,temp,[1 1 1;state_colours(nsta,:)],0.35,'k',{2 tempctrl 0.05},4);
    xlim([0.2 length(GroupStates)+0.8])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title('Perf')
    
    subplot(2,3,5)
    tempctrl=squeeze(mean(mean_byprobe_vig_ctrl(:,nsta),2));
    temp=squeeze(mean(mean_byprobe_vig(:,nsta),2));
    simpleBarPlot(nsta-0.2,tempctrl,[state_colours(nsta,:)],0.35,'k',[],4);
    simpleBarPlot(nsta+0.2,temp,[1 1 1;state_colours(nsta,:)],0.35,'k',{2 tempctrl 0.05},4);

    xlim([0.2 length(GroupStates)+0.8])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title('Vigilance')
 
%      subplot(2,3,5)
%     temp=squeeze(mean(mean_byprobe_(:,nsta),2));
%     simpleBarPlot(nsta,temp,state_colours(nsta,:),0.9,'k');
%     xlim([0.2 3.8])
%     set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
%     title('Not Looking')
end

%%
figure; format_fig;
for nval=1:5
    switch nval
        case 1
            tempctrl=squeeze(mean(mean_allprobes_awa_ctrl(:,nsta),2));
            temp=squeeze(mean(mean_allprobes_awa(:,nsta),2));
        case 2
            tempctrl=squeeze(mean(mean_allprobes_wil_ctrl(:,nsta),2));
            temp=squeeze(mean(mean_allprobes_wil(:,nsta),2));
        case 3
            tempctrl=squeeze(mean(mean_allprobes_eng_ctrl(:,nsta),2));
            temp=squeeze(mean(mean_allprobes_eng(:,nsta),2));
        case 4
            tempctrl=squeeze(mean(mean_allprobes_prf_ctrl(:,nsta),2));
            temp=squeeze(mean(mean_allprobes_prf(:,nsta),2));
        case 5
            tempctrl=squeeze(mean(mean_allprobes_vig_ctrl(:,nsta),2));
            temp=squeeze(mean(mean_allprobes_vig(:,nsta),2));
    end
     simpleBarPlot(nval-0.2,tempctrl,[0 0 0],0.35,'r',[],4);
    simpleBarPlot(nval+0.2,temp,[1 1 1;0 0 0],0.35,'r',{2 tempctrl 0.05},4);
end
xlim([0.2 5.8])
    set(gca,'XTick',1:5,'XTickLabel',{'AWAR','WILL','ENGT','PERF','VIGI'})
    
%%
nn = 1:100;
pwrout = sampsizepwr('t2',[mean(pption_MS_ctrl(:,3)) std(pption_MS_ctrl(:,3))],mean(pption_MS(:,3)),[],nn);

figure;
plot(nn,pwrout,'b-')
title('Power versus Sample Size')
xlabel('Sample Size')
ylabel('Power')

%%
nn = 1:100;
pwrout = sampsizepwr('t2',[mean(corr_nogo) std(corr_nogo)],mean(corr_nogo_ctrl),[],nn);

figure;
plot(nn,pwrout,'b-')
title('Power versus Sample Size')
xlabel('Sample Size')
ylabel('Power')