%%
clear all
% close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))

% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];

load([data_path filesep 'WanderIM_ASRSscores_15Feb19'])
%%
all_test_resprobes_perblock=[];
all_test_res=[];
all_probes_timing=[];
all_probes_mat=[];
RTs=[];
all_probes_mat_fullgo=[];
all_probes_mat_fullnogo=[];
nc=0;
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
    if sum(ASRS_scores(:,1)==str2num(SubID))==0
        fprintf('... missing ASRS value\n')
        continue
    end
    nc=nc+1;
   ASRS_scores2(nc,:)=ASRS_scores(ASRS_scores(:,1)==str2num(SubID),:);
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
    for nbt=1:2
        tp_nogos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,11)),11);
        tp_gos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,12)),12);
        [dprime_test(nc,nbt), crit_test(nc,nbt)]=calc_dprime((tp_gos==1),(tp_nogos==0));
        corr_go(nc,nbt)=nanmean(tp_gos);
        corr_nogo(nc,nbt)=nanmean(tp_nogos);
        
        rt_gos(nc,nbt)=nanmean(RTs(test_res(:,2)==nbt & ~isnan(test_res(:,12))));
        rt_nogos(nc,nbt)=nanmean(RTs(test_res(:,2)==nbt & ~isnan(test_res(:,11))));
    end
    all_test_res=[all_test_res ; [nc*ones(size(test_res,1),1) test_res(:,[1 2 4 11 12])]];
    
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
        pption_MS(nc,nbt,:)=pption;
        pption_Ori(nc,nbt,:)=pption2;
    end
    for nb=1:max(probe_res(:,4))
        temp=probe_res(probe_res(:,4)==nb,[5 31:38]);
        temp2=test_res(test_res(:,1)==nb,[11 12]);
        %             pption=[];
        %             for nstate=1:4
        %                 pption(nstate)=mean(temp(:,3)==nstate);
        %             end
        %             all_test_resprobes_perblock=[all_test_resprobes_perblock ; [nc 1 unique(temp(:,1)) pption mean(temp(:,[5:8])) nanmean(temp2)]];
        all_test_resprobes_perblock=[all_test_resprobes_perblock ; [nc*ones(size(temp,1),1) nb*ones(size(temp,1),1) temp]];
        all_probes_timing=[all_probes_timing ; diff(probe_res(probe_res(:,3)==nb,2))];
    end
    for nbt=1:2
        tp_probes=probe_res(probe_res(:,5)==nbt,31:38);
        mean_bytask_prf(nc,nbt)=nanmean(tp_probes(:,7));
        
        mean_byprobe2_awa(nc,nbt)=5-nanmean(tp_probes(:,4));
        mean_byprobe2_wil(nc,nbt)=5-nanmean(tp_probes(:,5));
        mean_byprobe2_eng(nc,nbt)=nanmean(tp_probes(:,6));
        mean_byprobe2_prf(nc,nbt)=nanmean(tp_probes(:,7));
        mean_byprobe2_vig(nc,nbt)=5-nanmean(tp_probes(:,8));
        
        for nstate=1:3
            mean_count_mindstates(nc,nbt,nstate)=sum(tp_probes(:,2)==(nstate));
            mean_count_ori(nc,nbt,nstate)=sum(tp_probes(:,3)==(nstate));
            
            if sum(tp_probes(:,2)==(nstate))~=0
                
                mean_byprobe_awa(nc,nbt,nstate)=5-nanmean(tp_probes(tp_probes(:,2)==(nstate),4));
                mean_byprobe_wil(nc,nbt,nstate)=5-nanmean(tp_probes(tp_probes(:,2)==(nstate),5));
                mean_byprobe_eng(nc,nbt,nstate)=nanmean(tp_probes(tp_probes(:,2)==(nstate),6));
                mean_byprobe_prf(nc,nbt,nstate)=nanmean(tp_probes(tp_probes(:,2)==(nstate),7));
                mean_byprobe_vig(nc,nbt,nstate)=5-nanmean(tp_probes(tp_probes(:,2)==(nstate),8));
            else
                mean_byprobe_awa(nc,nbt,nstate)=nan;
                mean_byprobe_wil(nc,nbt,nstate)=nan;
                mean_byprobe_eng(nc,nbt,nstate)=nan;
                mean_byprobe_prf(nc,nbt,nstate)=nan;
                mean_byprobe_vig(nc,nbt,nstate)=nan;
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
            last_pr_tridx=this_pr_tridx-21;
            number_STD(nc,nbl,npr)=sum(these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,5)~=3);
            number_DEV(nc,nbl,npr)=sum(these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,5)==3);
            
            temp_testres=these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,:);
            %             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            tcorr_go=nanmean(temp_testres(:,12));%/corr_go(nc,these_probes(npr,5));
            tcorr_nogo=nanmean(temp_testres(:,11));%/corr_nogo(nc,these_probes(npr,5));
            num_go=sum(~isnan(temp_testres(:,12)));
            num_nogo=sum(~isnan(temp_testres(:,11)));
            rt_go=nanmean(temp_testres(~isnan(temp_testres(:,12)),10)-temp_testres(~isnan(temp_testres(:,12)),8));
            rt_nogo=nanmean(temp_testres(~isnan(temp_testres(:,11)),10)-temp_testres(~isnan(temp_testres(:,11)),8));
            
            tdp=calc_dprime((temp_testres(~isnan(temp_testres(:,12)),12)==1),(temp_testres(~isnan(temp_testres(:,11)),11)==0));
            all_probes_mat=[all_probes_mat ; [nc nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_go tcorr_nogo num_go num_nogo tdp (these_probes(npr,37)) rt_go rt_nogo]];
            
            temp_testres=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,12)),:);
            tcorr_go=(temp_testres(end-19:end,12))';%/corr_go(nc,these_probes(npr,5));
             temp_testres=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,11)),:);
            tcorr_nogo=(temp_testres(end-1:end,11))';%/corr_go(nc,these_probes(npr,5));
            all_probes_mat_fullgo=[all_probes_mat_fullgo ; [nc nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_go]];
            all_probes_mat_fullnogo=[all_probes_mat_fullnogo ; [nc nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_nogo]];
        end
    end
end

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
ylim([0.4 0.6])

subplot(3,2,6); format_fig;
simpleBarPlot(1,rt_nogos(:,1),cond_colours(1,:),0.9,'k');
simpleBarPlot(2,rt_nogos(:,2),cond_colours(2,:),0.9,'k');
set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
ylabel('RT. NO-GO')
xlim([0.2 2.8])
ylim([0.4 0.6])
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
figure;
format_fig;
StateNames={'ON','MW','MB'};
for nst=1:3
    subplot(1,3,nst);
    simpleCorPlot(ASRS_scores2(:,4),mean(pption_MS(:,:,nst),2),{'o',state_colours(nst,:),state_colours(nst,:),72},'pearson',0);
    
    format_fig;
    % title(StateNames{nst})
    ylabel(sprintf('%% of %s',StateNames{nst}))
    xlabel('ASRS score (A+B)')
end

% %%
% figure;
% format_fig;
% thisColu=4;
% StateNames={'ON','MW','MB'};
% for nst=1:3
%     subplot(1,3,nst);
%     if nst==1
%     simpleCorPlot(ASRS_scores2(:,thisColu),mean(dprime_test,2),{'o',state_colours(nst,:),state_colours(nst,:),72},'pearson',0);
%     elseif nst==2
%     simpleCorPlot(ASRS_scores2(:,thisColu),mean(crit_test,2),{'o',state_colours(nst,:),state_colours(nst,:),72},'pearson',0);
%     elseif nst==3
%     simpleCorPlot(ASRS_scores2(:,thisColu),mean(rt_gos,2),{'o',state_colours(nst,:),state_colours(nst,:),72},'pearson',0);
%     end
%     format_fig;
%     % title(StateNames{nst})
%     ylabel(sprintf('% of %s',StateNames{nst}))
%     xlabel('ASRS score (A+B)')
% end