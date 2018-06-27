%%
clear all
close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))

% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];
%%
all_test_resprobes_perblock=[];
all_test_res=[];
all_probes_timing=[];
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    
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
    for nbt=1:2
        tp_nogos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,11)),11);
        tp_gos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,12)),12);
        [dprime_test(n,nbt), crit_test(n,nbt)]=calc_dprime((tp_gos==1),(tp_nogos==0));
        corr_go(n,nbt)=nanmean(tp_gos);
        corr_nogo(n,nbt)=nanmean(tp_nogos);
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
            pption2(nstate)=mean(temp(:,4)==nstate);
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
end

%% Behaviour - Performance on SART
figure; set(gcf,'Position',[529   943   615   750]);
subplot(2,2,1); format_fig;
simpleBarPlot(1,dprime_test(:,1),cond_colours(1,:),0.9,'k');
simpleBarPlot(2,dprime_test(:,2),cond_colours(2,:),0.9,'k');
set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
ylabel('d-prime')
xlim([0.2 2.8])

subplot(2,2,2); format_fig;
simpleBarPlot(1,crit_test(:,1),cond_colours(1,:),0.9,'k');
simpleBarPlot(2,crit_test(:,2),cond_colours(2,:),0.9,'k');
set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
ylabel('criterion')
xlim([0.2 2.8])

subplot(2,2,3); format_fig;
simpleBarPlot(1,100*corr_go(:,1),cond_colours(1,:),0.9,'k');
simpleBarPlot(2,100*corr_go(:,2),cond_colours(2,:),0.9,'k');
set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
ylabel('Perf. GO')
xlim([0.2 2.8])
ylim([40 100])

subplot(2,2,4); format_fig;
simpleBarPlot(1,100*corr_nogo(:,1),cond_colours(1,:),0.9,'k');
simpleBarPlot(2,100*corr_nogo(:,2),cond_colours(2,:),0.9,'k');
set(gca,'XTick',1:2,'XTickLabel',{'FACE','DIGIT'})
ylabel('Perf. NO-GO')
xlim([0.2 2.8])
ylim([40 100])

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
