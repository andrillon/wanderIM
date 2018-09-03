%%
clear all
close all

run localdef_wanderIM

addpath(genpath(lscpTools_path))

% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];

%load([data_path filesep 'CARS_quest'])
%%
all_test_resprobes_perblock=[];
all_test_res=[];
all_probes_timing=[];
all_probes_mat=[];
RTs=[];
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
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
    RTs=[RTs; test_res(:,10)-test_res(:,8)];
    %     findFalseStarts=find(test_res(:,10)-test_res(:,8)<0.1);
    %     test_res(findFalseStarts,[10 11 12])=NaN;
    %     warning('correcting for false starts')
    for nbt=1:2
        tp_nogos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,11)),11);
        tp_gos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,12)),12);
        [dprime_test(n,nbt), crit_test(n,nbt)]=calc_dprime((tp_gos==1),(tp_nogos==0));
        corr_go(n,nbt)=nanmean(tp_gos);
        corr_nogo(n,nbt)=nanmean(tp_nogos);
        
        rt_gos(n,nbt)=nanmean(RTs(test_res(:,2)==nbt & ~isnan(test_res(:,12))));
        rt_nogos(n,nbt)=nanmean(RTs(test_res(:,2)==nbt & ~isnan(test_res(:,11))));
    end
    all_test_res=[all_test_res ; [n*ones(size(test_res,1),1) test_res(:,[1 2 4 11 12])]];
    
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
    for nbt=1:2
        temp=probe_res(probe_res(:,5)==nbt,[5 31:38]);
        clear pption pption2
        for nstate=1:4
            pption(nstate)=mean(temp(:,3)==nstate);
            pption2(nstate)=mean(temp(temp(:,3)==2,4)==nstate);
        end
        pption_MS(n,nbt,:)=pption;
        pption_Ori(n,nbt,:)=pption2;
    end
    for nb=1:max(probe_res(:,3))
        temp=probe_res(probe_res(:,4)==nb,[5 31:38]);
        temp2=test_res(test_res(:,1)==nb,[11 12]);
        %             pption=[];
        %             for nstate=1:4
        %                 pption(nstate)=mean(temp(:,3)==nstate);
        %             end
        %             all_test_resprobes_perblock=[all_test_resprobes_perblock ; [n 1 unique(temp(:,1)) pption mean(temp(:,[5:8])) nanmean(temp2)]];
        all_test_resprobes_perblock=[all_test_resprobes_perblock ; [n*ones(size(temp,1),1) nb*ones(size(temp,1),1) temp]];
        all_probes_timing=[all_probes_timing ; diff(probe_res(probe_res(:,3)==nb,2))];
    end
    for nbt=1:2
        tp_probes=probe_res(probe_res(:,5)==nbt,31:38);
        mean_bytask_prf(n,nbt)=nanmean(tp_probes(:,7));
        
        
        for nstate=1:3
            mean_count_mindstates(n,nbt,nstate)=sum(tp_probes(:,2)==(nstate));
            mean_count_ori(n,nbt,nstate)=sum(tp_probes(:,3)==(nstate));
            
            if sum(tp_probes(:,2)==(nstate))~=0
                
                mean_byprobe_awa(n,nbt,nstate)=5-nanmean(tp_probes(tp_probes(:,2)==(nstate),4));
                mean_byprobe_wil(n,nbt,nstate)=5-nanmean(tp_probes(tp_probes(:,2)==(nstate),5));
                mean_byprobe_eng(n,nbt,nstate)=nanmean(tp_probes(tp_probes(:,2)==(nstate),6));
                mean_byprobe_prf(n,nbt,nstate)=nanmean(tp_probes(tp_probes(:,2)==(nstate),7));
                mean_byprobe_vig(n,nbt,nstate)=5-nanmean(tp_probes(tp_probes(:,2)==(nstate),8));
            else
                mean_byprobe_awa(n,nbt,nstate)=nan;
                mean_byprobe_wil(n,nbt,nstate)=nan;
                mean_byprobe_eng(n,nbt,nstate)=nan;
                mean_byprobe_prf(n,nbt,nstate)=nan;
                mean_byprobe_vig(n,nbt,nstate)=nan;
            end
        end
    end
    
    %%% performance between probes
    for nbl=1:6
        these_probes=probe_res(probe_res(:,4)==nbl,:);
        these_trials=test_res(test_res(:,1)==nbl,:);
        for npr=1:10
            this_pr_tridx=these_probes(npr,6);
            %             if npr==1
            %                 last_pr_tridx=0;
            %             else
            %                 last_pr_tridx=these_probes(npr-1,6);
            %             end
            last_pr_tridx=20;
            number_STD(n,nbl,npr)=sum(these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,5)~=3);
            number_DEV(n,nbl,npr)=sum(these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,5)==3);
            
            temp_testres=these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,:);
            %             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            tcorr_go=nanmean(temp_testres(:,12));
            tcorr_nogo=nanmean(temp_testres(:,11));
            num_go=sum(~isnan(temp_testres(:,12)));
            num_nogo=sum(~isnan(temp_testres(:,11)));
            rt_go=nanmean(temp_testres(~isnan(temp_testres(:,12)),10)-temp_testres(~isnan(temp_testres(:,12)),8));
            rt_nogo=nanmean(temp_testres(~isnan(temp_testres(:,11)),10)-temp_testres(~isnan(temp_testres(:,11)),8));
            
            tdp=calc_dprime((temp_testres(~isnan(temp_testres(:,12)),12)==1),(temp_testres(~isnan(temp_testres(:,11)),11)==0));
            all_probes_mat=[all_probes_mat ; [n nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_go tcorr_nogo num_go num_nogo tdp (these_probes(npr,37)) rt_go rt_nogo]];
        end
    end
end

%% Average probe results across subjects
mysub=unique(all_probes_mat(:,1));
for ns=1:length(mysub)
    tp_probes_mat=all_probes_mat(all_probes_mat(:,1)==mysub(ns),:);
    
    % by taks
    for nt=1:2
        % by probe type
        for npr=1:3
            tp_probes_thistype=tp_probes_mat(tp_probes_mat(:,4)==npr & tp_probes_mat(:,3)==nt,:);
            corr_GO_byprobe_bysubj(ns,nt,npr)=nanmean(tp_probes_thistype(:,6));
            corr_NOGO_byprobe_bysubj(ns,nt,npr)=nanmean(tp_probes_thistype(:,7));
            
            rt_GO_byprobe_bysubj(ns,nt,npr)=nanmean(tp_probes_thistype(:,12));
            rt_NOGO_byprobe_bysubj(ns,nt,npr)=nanmean(tp_probes_thistype(:,13));
            
            n_GO_byprobe_bysubj(ns,nt,npr)=length(tp_probes_thistype(:,6));
            n_NOGO_byprobe_bysubj(ns,nt,npr)=length(tp_probes_thistype(:,7));
        end
    end
end

%% Average probe results across probes
% by taks
for nt=1:2
    % by probe type
    for npr=1:3
        tp_probes_thistype=all_probes_mat(all_probes_mat(:,4)==npr & all_probes_mat(:,3)==nt,:);
        
        corr_GO_byprobe{npr,nt}=(tp_probes_thistype(:,6));
        corr_NOGO_byprobe{npr,nt}=(tp_probes_thistype(:,7));
        
        rt_GO_byprobe{npr,nt}=(tp_probes_thistype(:,12));
        rt_NOGO_byprobe{npr,nt}=(tp_probes_thistype(:,13));
    end
end
%% Behaviour - Performance on SART - BY MIND-STATE AVGD ACROSS PROBES & PARTICIPANTS
% 
% figure;  set(gcf,'Position',[-29        1275        1218         412]);
% subplot(2,2,1); format_fig;
% for nst=1:3
%     hb(1)=simpleBarPlot(nst-0.2,100*squeeze(corr_GO_byprobe_bysubj(:,1,nst)),state_colours(nst,:),0.35,'k');
%     hb(2)=simpleBarPlot(nst+0.2,100*squeeze(corr_GO_byprobe_bysubj(:,2,nst)),[1 1 1 ; state_colours(nst,:)],0.35,'k');
% end
% set(gca,'XTick',1:3,'XTickLabel',{'On','MW','MB'})
% ylabel('Corr Go')
% title('Correctness by Mind-State')
% xlim([0.2 3.8])
% ylim([50 100])
% legend(hb,{'Face','Square'})
% 
% subplot(2,2,2); format_fig;
% for nst=1:3
%     hb(1)=simpleBarPlot(nst-0.2,100*squeeze(corr_NOGO_byprobe_bysubj(:,1,nst)),state_colours(nst,:),0.35,'k');
%     hb(2)=simpleBarPlot(nst+0.2,100*squeeze(corr_NOGO_byprobe_bysubj(:,2,nst)),[1 1 1 ; state_colours(nst,:)],0.35,'k');
% end
% set(gca,'XTick',1:3,'XTickLabel',{'On','MW','MB'})
% ylabel('Corr No-Go')
% title('Correctness by Mind-State')
% xlim([0.2 3.8])
% ylim([50 100])
% legend(hb,{'Face','Square'})
% 
% subplot(2,2,3); format_fig;
% for nst=1:3
%     hb(1)=simpleBarPlot(nst-0.2,squeeze(rt_GO_byprobe_bysubj(:,1,nst)),state_colours(nst,:),0.35,'k');
%     hb(2)=simpleBarPlot(nst+0.2,squeeze(rt_GO_byprobe_bysubj(:,2,nst)),[1 1 1 ; state_colours(nst,:)],0.35,'k');
% end
% set(gca,'XTick',1:3,'XTickLabel',{'On','MW','MB'})
% ylabel('RT Go')
% title('Reaction Time by Mind-State')
% xlim([0.2 3.8])
% ylim([0.2 0.6])
% legend(hb,{'Face','Square'})
% 
% subplot(2,2,4); format_fig;
% for nst=1:3
%     hb(1)=simpleBarPlot(nst-0.2,squeeze(rt_NOGO_byprobe_bysubj(:,1,nst)),state_colours(nst,:),0.35,'k');
%     hb(2)=simpleBarPlot(nst+0.2,squeeze(rt_NOGO_byprobe_bysubj(:,2,nst)),[1 1 1 ; state_colours(nst,:)],0.35,'k');
% end
% set(gca,'XTick',1:3,'XTickLabel',{'On','MW','MB'})
% ylabel('RT No-Go')
% title('Reaction Time by Mind-State')
% xlim([0.2 3.8])
% ylim([0.2 0.6])
% legend(hb,{'Face','Square'})


%cond_colours


% figure;  set(gcf,'Position',[-29        1275        1218         412]);
% subplot(2,2,1); format_fig;
% for nst=1:3
%     hb(1)=simpleBarPlot(nst-0.2,100*squeeze(corr_GO_byprobe_bysubj(:,1,nst)),state_colours(nst,:),0.35,'k');
%     hb(2)=simpleBarPlot(nst+0.2,100*squeeze(corr_GO_byprobe_bysubj(:,2,nst)),[1 1 1 ; state_colours(nst,:)],0.35,'k');
% end
% set(gca,'XTick',1:3,'XTickLabel',{'On','MW','MB'})
% ylabel('Corr Go')
% title('Correctness by Mind-State')
% xlim([0.2 3.8])
% ylim([50 100])
% legend(hb,{'Face','Square'})

figure;  set(gcf,'Position',[-29        1275        1218         412]);
subplot(2,2,1); format_fig;
hb(1)=simpleBarPlot(0.6,100*squeeze(corr_GO_byprobe_bysubj(:,1,1)),state_colours(1,:),0.35,'k');
hb(2)=simpleBarPlot(1.0,100*squeeze(corr_GO_byprobe_bysubj(:,1,2)),state_colours(2,:),0.35,'k');
hb(3)=simpleBarPlot(1.4,100*squeeze(corr_GO_byprobe_bysubj(:,1,3)),state_colours(3,:),0.35,'k');
hb(1)=simpleBarPlot(2.0,100*squeeze(corr_GO_byprobe_bysubj(:,2,1)),[1 1 1 ; state_colours(1,:)],0.35,'k');
hb(2)=simpleBarPlot(2.4,100*squeeze(corr_GO_byprobe_bysubj(:,2,2)),[1 1 1 ; state_colours(2,:)],0.35,'k');
hb(3)=simpleBarPlot(2.8,100*squeeze(corr_GO_byprobe_bysubj(:,2,3)),[1 1 1 ; state_colours(3,:)],0.35,'k');
set(gca,'XTick',1:2,'XTickLabel',{'Face','Square'})
ylabel('Corr Go')
title('Correctness by Task')
xlim([0.2 4.4])
ylim([90 100])
legend(hb,{'ON','MW', 'MB'})

subplot(2,2,2); format_fig;
hb(1)=simpleBarPlot(0.6,100*squeeze(corr_NOGO_byprobe_bysubj(:,1,1)),state_colours(1,:),0.35,'k');
hb(2)=simpleBarPlot(1.0,100*squeeze(corr_NOGO_byprobe_bysubj(:,1,2)),state_colours(2,:),0.35,'k');
hb(3)=simpleBarPlot(1.4,100*squeeze(corr_NOGO_byprobe_bysubj(:,1,3)),state_colours(3,:),0.35,'k');
hb(1)=simpleBarPlot(2.0,100*squeeze(corr_NOGO_byprobe_bysubj(:,2,1)),[1 1 1 ; state_colours(1,:)],0.35,'k');
hb(2)=simpleBarPlot(2.4,100*squeeze(corr_NOGO_byprobe_bysubj(:,2,2)),[1 1 1 ; state_colours(2,:)],0.35,'k');
hb(3)=simpleBarPlot(2.8,100*squeeze(corr_NOGO_byprobe_bysubj(:,2,3)),[1 1 1 ; state_colours(3,:)],0.35,'k');
set(gca,'XTick',1:2,'XTickLabel',{'Face','Square'})
ylabel('Corr No-Go')
title('Correctness by Task')
xlim([0.2 4.4])
ylim([50 80])
legend(hb,{'ON','MW', 'MB'})

subplot(2,2,3); format_fig;
hb(1)=simpleBarPlot(0.6,100*squeeze(rt_GO_byprobe_bysubj(:,1,1)),state_colours(1,:),0.35,'k');
hb(2)=simpleBarPlot(1.0,100*squeeze(rt_GO_byprobe_bysubj(:,1,2)),state_colours(2,:),0.35,'k');
hb(3)=simpleBarPlot(1.4,100*squeeze(rt_GO_byprobe_bysubj(:,1,3)),state_colours(3,:),0.35,'k');
hb(1)=simpleBarPlot(2.0,100*squeeze(rt_GO_byprobe_bysubj(:,2,1)),[1 1 1 ; state_colours(1,:)],0.35,'k');
hb(2)=simpleBarPlot(2.4,100*squeeze(rt_GO_byprobe_bysubj(:,2,2)),[1 1 1 ; state_colours(2,:)],0.35,'k');
hb(3)=simpleBarPlot(2.8,100*squeeze(rt_GO_byprobe_bysubj(:,2,3)),[1 1 1 ; state_colours(3,:)],0.35,'k');
set(gca,'XTick',1:2,'XTickLabel',{'Face','Square'})
ylabel('RT Go')
title('Reaction Time by Task')
xlim([0.2 4.4])
ylim([40 60])
legend(hb,{'ON','MW', 'MB'})

subplot(2,2,4); format_fig;
hb(1)=simpleBarPlot(0.6,100*squeeze(rt_NOGO_byprobe_bysubj(:,1,1)),state_colours(1,:),0.35,'k');
hb(2)=simpleBarPlot(1.0,100*squeeze(rt_NOGO_byprobe_bysubj(:,1,2)),state_colours(2,:),0.35,'k');
hb(3)=simpleBarPlot(1.4,100*squeeze(rt_NOGO_byprobe_bysubj(:,1,3)),state_colours(3,:),0.35,'k');
hb(1)=simpleBarPlot(2.0,100*squeeze(rt_NOGO_byprobe_bysubj(:,2,1)),[1 1 1 ; state_colours(1,:)],0.35,'k');
hb(2)=simpleBarPlot(2.4,100*squeeze(rt_NOGO_byprobe_bysubj(:,2,2)),[1 1 1 ; state_colours(2,:)],0.35,'k');
hb(3)=simpleBarPlot(2.8,100*squeeze(rt_NOGO_byprobe_bysubj(:,2,3)),[1 1 1 ; state_colours(3,:)],0.35,'k');
set(gca,'XTick',1:2,'XTickLabel',{'Face','Square'})
ylabel('RT No-Go')
title('Reaction Time by Task')
xlim([0.2 4.4])
ylim([30 50])
legend(hb,{'ON','MW', 'MB'})


% figure; set(gcf,'Position',[529   943   615   750]);
% subplot(2,2,1); format_fig;
% simpleBarPlot(1,100*corr_GO_byprobe_bysubj(:,1,1),state_colours(:,1),0.35,'k');
% simpleBarPlot(2,100*corr_GO_byprobe_bysubj(:,1,2),state_colours(:,2),0.35,'k');
% simpleBarPlot(3,100*corr_GO_byprobe_bysubj(:,1,3),state_colours(:,3),0.35,'k');
% set(gca,'XTick',1:2,'XTickLabel',{'ON','MW','MB'})
% ylabel('Perf. GO')
% xlim([0.2 6])
% ylim([40 100])
% 
% subplot(2,2,2); format_fig;
% simpleBarPlot(1,100*corr_NOGO_byprobe_bysubj(:,1),cond_colours(1,:),0.9,'k');
% simpleBarPlot(2,100*corr_NOGO_byprobe_bysubj(:,2),cond_colours(2,:),0.9,'k');
% set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
% ylabel('Perf. NO-GO')
% xlim([0.2 2.8])
% ylim([40 100])
% 
% subplot(2,2,3); format_fig;
% simpleBarPlot(1,rt_GO_byprobe_bysubj(:,1),cond_colours(1,:),0.9,'k');
% simpleBarPlot(2,rt_GO_byprobe_bysubj(:,2),cond_colours(2,:),0.9,'k');
% set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
% ylabel('Perf. GO')
% xlim([0.2 2.8])
% ylim([0 1])
% 
% subplot(2,2,4); format_fig;
% simpleBarPlot(1,rt_NOGO_byprobe_bysubj(:,1),cond_colours(1,:),0.9,'k');
% simpleBarPlot(2,rt_NOGO_byprobe_bysubj(:,2),cond_colours(2,:),0.9,'k');
% set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
% ylabel('Perf. NO-GO')
% xlim([0.2 2.8])
% ylim([0 1])



%% Behaviour - Performance on SART

figure; set(gcf,'Position',[529   943   615   750]);
subplot(3,2,1); format_fig;
simpleBarPlot(1,dprime_test(:,1),cond_colours(1,:),0.9,'k');
simpleBarPlot(2,dprime_test(:,2),cond_colours(2,:),0.9,'k');
set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
ylabel('d-prime')
xlim([0.2 2.8])

subplot(3,2,2); format_fig;
simpleBarPlot(1,crit_test(:,1),cond_colours(1,:),0.9,'k');
simpleBarPlot(2,crit_test(:,2),cond_colours(2,:),0.9,'k');
set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
ylabel('criterion')
xlim([0.2 2.8])

subplot(3,2,3); format_fig;
simpleBarPlot(1,100*corr_go(:,1),cond_colours(1,:),0.9,'k');
simpleBarPlot(2,100*corr_go(:,2),cond_colours(2,:),0.9,'k');
set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
ylabel('Perf. GO')
xlim([0.2 2.8])
ylim([90 100])

subplot(3,2,4); format_fig;
simpleBarPlot(1,100*corr_nogo(:,1),cond_colours(1,:),0.9,'k');
simpleBarPlot(2,100*corr_nogo(:,2),cond_colours(2,:),0.9,'k');
set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
ylabel('Perf. NO-GO')
xlim([0.2 2.8])
ylim([50 80])

subplot(3,2,5); format_fig;
simpleBarPlot(1,rt_gos(:,1),cond_colours(1,:),0.9,'k');
simpleBarPlot(2,rt_gos(:,2),cond_colours(2,:),0.9,'k');
set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
ylabel('RT. GO')
xlim([0.2 2.8])
ylim([0.5 0.6])

subplot(3,2,6); format_fig;
simpleBarPlot(1,rt_nogos(:,1),cond_colours(1,:),0.9,'k');
simpleBarPlot(2,rt_nogos(:,2),cond_colours(2,:),0.9,'k');
set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
ylabel('RT. NO-GO')
xlim([0.2 2.8])
ylim([0.5 0.6])
% 
% figure; set(gcf,'Position',[529   943   615   750]);
% subplot(2,2,1); format_fig;
% simpleBarPlot(1,dprime_test(:,1),cond_colours(1,:),0.9,'k');
% simpleBarPlot(2,dprime_test(:,2),cond_colours(2,:),0.9,'k');
% set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
% ylabel('d-prime')
% xlim([0.2 2.8])
% 
% subplot(2,2,2); format_fig;
% simpleBarPlot(1,crit_test(:,1),cond_colours(1,:),0.9,'k');
% simpleBarPlot(2,crit_test(:,2),cond_colours(2,:),0.9,'k');
% set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
% ylabel('criterion')
% xlim([0.2 2.8])
% 
% subplot(2,2,3); format_fig;
% simpleBarPlot(1,100*corr_go(:,1),cond_colours(1,:),0.9,'k');
% simpleBarPlot(2,100*corr_go(:,2),cond_colours(2,:),0.9,'k');
% set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
% ylabel('Perf. GO')
% xlim([0.2 2.8])
% ylim([40 100])
% 
% subplot(2,2,4); format_fig;
% simpleBarPlot(1,100*corr_nogo(:,1),cond_colours(1,:),0.9,'k');
% simpleBarPlot(2,100*corr_nogo(:,2),cond_colours(2,:),0.9,'k');
% set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
% ylabel('Perf. NO-GO')
% xlim([0.2 2.8])
% ylim([40 100])

%% Probes - % mind-state
figure;  set(gcf,'Position',[-29        1275        1218         412]);
subplot(1,2,1); format_fig;
for nst=1:3
    hb(1)=simpleBarPlot(nst-0.2,100*squeeze(pption_MS(:,1,nst)),state_colours(nst,:),0.35,'k');
    hb(2)=simpleBarPlot(nst+0.2,100*squeeze(pption_MS(:,2,nst)),[1 1 1 ; state_colours(nst,:)],0.35,'k');
end
set(gca,'XTick',1:3,'XTickLabel',{'On','MW','MB'})
ylabel('% probes')
title('Mind-State')
xlim([0.2 3.8])
legend(hb,{'Face','Square'})

subplot(1,2,2); format_fig;
for nst=1:3
    simpleBarPlot(nst-0.2,100*squeeze(pption_Ori(:,1,nst)),state_colours(nst,:),0.35,'k');
    simpleBarPlot(nst+0.2,100*squeeze(pption_Ori(:,2,nst)),[1 1 1 ; state_colours(nst,:)],0.35,'k');
end
set(gca,'XTick',1:3,'XTickLabel',{'room','perso','task'})
ylabel('% probes')
title('Origin of MW')
xlim([0.2 3.8])

%%
figure; format_fig;
for nsta=1:3
    subplot(2,3,1)
    temp=squeeze(mean(mean_byprobe_awa(:,:,nsta),2));
    simpleBarPlot(nsta,temp,state_colours(nsta,:),0.9,'k');
    xlim([0.2 3.8])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title('Awareness')
    
    subplot(2,3,2)
    temp=squeeze(mean(mean_byprobe_wil(:,:,nsta),2));
    simpleBarPlot(nsta,temp,state_colours(nsta,:),0.9,'k');
    xlim([0.2 3.8])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title('Intention')
    
    subplot(2,3,3)
    temp=squeeze(mean(mean_byprobe_eng(:,:,nsta),2));
    simpleBarPlot(nsta,temp,state_colours(nsta,:),0.9,'k');
    xlim([0.2 3.8])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title('Engagement')
    
    subplot(2,3,4)
    temp=squeeze(mean(mean_byprobe_prf(:,:,nsta),2));
    simpleBarPlot(nsta,temp,state_colours(nsta,:),0.9,'k');
    xlim([0.2 3.8])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title('Perf')
    
    subplot(2,3,5)
    temp=squeeze(mean(mean_byprobe_vig(:,:,nsta),2));
    simpleBarPlot(nsta,temp,state_colours(nsta,:),0.9,'k');
    xlim([0.2 3.8])
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title('Vigilance')
    
    %      subplot(2,3,6)
    %     temp=squeeze(mean(mean_byprobe_look(:,:,nsta),2));
    %     simpleBarPlot(nsta,temp,state_colours(nsta,:),0.9,'k');
    %     xlim([0.2 3.8])
    %     set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    %     title('Looking')
    
end

%%
figure; format_fig;
radar_data=[nanmean(squeeze(mean(mean_byprobe_awa(:,:,:),2))); nanmean(squeeze(mean(mean_byprobe_wil(:,:,:),2))); nanmean(squeeze(mean(mean_byprobe_eng(:,:,:),2))); nanmean(squeeze(mean(mean_byprobe_prf(:,:,:),2))); nanmean(squeeze(mean(mean_byprobe_vig(:,:,:),2)))];
%radar_colors={'r', 'g', 'b'};
radar_colors={state_colours(1, :), state_colours(2, :), state_colours(3,:)};
radarplot(transpose(radar_data), {'Engagement','Awareness','Intention','Vigilance','Performance'}, radar_colors, radar_colors);
legend('on-tasks', 'mind-wandering', 'mind-blanking');

%%
tbl=array2table(all_test_resprobes_perblock,'VariableNames',{'SubID','nBlock','Task','Look','State','Ori','Awa','Wil','Enga','Perf','Vig'});
tbl.SubID=categorical(tbl.SubID);
tbl.Task=categorical(tbl.Task);
tbl.Look=categorical(tbl.Look);
tbl.State=categorical(tbl.State);
tbl.Ori=categorical(tbl.Ori);

tbl.On=zeros(size(all_test_resprobes_perblock,1),1);
tbl.On(tbl.State=="1")=1;
tbl.MBMW=nan(size(all_test_resprobes_perblock,1),1);
tbl.MBMW(tbl.State=="2")=1;
tbl.MBMW(tbl.State=="3")=0;
% tbl.On=categorical(tbl.On);

lme_0= fitlme(tbl,'Perf~1+(1|SubID)');
lme_1= fitlme(tbl,'Perf~State+(1|SubID)');

%%
tbl2=array2table(all_test_res,'VariableNames',{'SubID','nBlock','Task','nTrial','Go','NoGo'});
tbl2.SubID=categorical(tbl2.SubID);
tbl2.Task=categorical(tbl2.Task);

lme_1= fitlme(tbl2,'Go~Task+nBlock+nTrial+(1|SubID)');

lme_0= fitlme(tbl2,'Go~1+(1|SubID)');
lme_1= fitlme(tbl2,'Go~Task+(1|SubID)');

%%
tblPr=array2table(all_probes_mat,'VariableNames',{'SubID','nBlock','Task','State','nProbe','Go','NoGo','nGo','nNoGo','dprime','PerfEst'});
tblPr(tblPr.State==4,:)=[];
tblPr.SubID=categorical(tblPr.SubID);
tblPr.State=categorical(tblPr.State);
tblPr.Task=categorical(tblPr.Task);

lme_NoGo= fitlme(tblPr,'NoGo~Task+nBlock+nProbe+State+(1|SubID)');
lme_Go= fitlme(tblPr,'Go~Task+nBlock+nProbe+State+(1|SubID)');

subtblPr=tblPr(ismember(tblPr.State,{'2','3'}),:);
subtblPr.State=categorical(subtblPr.State);
subtblPr.State=removecats(subtblPr.State);
lme2_NoGo= fitlme(subtblPr,'NoGo~Task+nBlock+nProbe+State+(1|SubID)');
lme2_Go= fitlme(subtblPr,'Go~Task+nBlock+nProbe+State+(1|SubID)');

addpath(genpath(gramm_path))

clear g

g(1,1)=gramm('x',double(tblPr.Task),'y',tblPr.Go,'color',double(tblPr.State),'subset',double(tblPr.State)~=4);
g(1,2)=copy(g(1));
g(2,1)=copy(g(1));
g(2,2)=copy(g(1));

%Jittered scatter plot
g(1,1).geom_jitter('width',0.4,'height',0);
g(1,1).set_title('geom_jitter()');
g(1,1).axe_property('XLim',[0 3],'XTick',[ 1 2],'XTickLabel',{'FACE','DIGIT'});

%Averages with confidence interval
g(1,2).stat_summary('geom',{'bar','black_errorbar'});
g(1,2).set_title('stat_summary()');
g(1,2).axe_property('XLim',[0 3],'XTick',[ 1 2],'XTickLabel',{'FACE','DIGIT'});

%Boxplots
g(2,1).stat_boxplot();
g(2,1).set_title('stat_boxplot()');
g(2,1).axe_property('XLim',[0 3],'XTick',[ 1 2],'XTickLabel',{'FACE','DIGIT'});

%Violin plots
g(2,2).stat_violin('fill','transparent');
g(2,2).set_title('stat_violin()');
g(2,2).axe_property('XLim',[0 3],'XTick',[ 1 2],'XTickLabel',{'FACE','DIGIT'});

%These functions can be called on arrays of gramm objects
g.set_names('x','Task','y','NoGO','color','State');
g.set_title('Visualization of Y~X relationships with X as categorical variable');

gf = copy(g);

figure('Position',[100 100 800 550]);
g.draw();

gf.set_title('Visualization of Y~X relationships with X as categorical variable and flipped coordinates');

%%
Y1 = splitapply(@mean,tblPr.Go(tblPr.Task=='1',:),double(tblPr.SubID(tblPr.Task=='1',:)));
Y2 = splitapply(@mean,tblPr.Go(tblPr.Task=='2',:),double(tblPr.SubID(tblPr.Task=='2',:)));
figure; subplot(1,3,1)
simpleCorPlot(Y1,Y2,[],'pearson',0);

Y1 = splitapply(@mean,tblPr.NoGo(tblPr.Task=='1',:),double(tblPr.SubID(tblPr.Task=='1',:)));
Y2 = splitapply(@mean,tblPr.NoGo(tblPr.Task=='2',:),double(tblPr.SubID(tblPr.Task=='2',:)));
subplot(1,3,2)
simpleCorPlot(Y1,Y2,[],'pearson',0);

Y1 = splitapply(@mean,tblPr.dprime(tblPr.Task=='1',:),double(tblPr.SubID(tblPr.Task=='1',:)));
Y2 = splitapply(@mean,tblPr.dprime(tblPr.Task=='2',:),double(tblPr.SubID(tblPr.Task=='2',:)));
subplot(1,3,3)
simpleCorPlot(Y1,Y2,[],'pearson',0);

%%
R1 = splitapply(@corr,all_probes_mat(all_probes_mat(:,3)==1,11),all_probes_mat(all_probes_mat(:,3)==1,10),all_probes_mat(all_probes_mat(:,3)==1,1));
R2 = splitapply(@corr,all_probes_mat(all_probes_mat(:,3)==2,11),all_probes_mat(all_probes_mat(:,3)==2,10),all_probes_mat(all_probes_mat(:,3)==2,1));
