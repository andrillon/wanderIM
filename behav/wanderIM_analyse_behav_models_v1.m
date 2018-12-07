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

%load([data_path filesep 'CARS_quest'])
%%
all_test_resprobes_perblock=[];
all_test_res=[];
all_probes_timing=[];
all_probes_mat=[];
RTs=[];
all_probes_mat_fullgo=[];
all_probes_mat_fullnogo=[];
all_probes_mat_fullrtgo=[];
all_probes_mat_fullrtnogo=[];
all_probes_mat2=[];
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
    RTs=[test_res(:,10)-test_res(:,8)];
    %     findFalseStarts=find(test_res(:,10)-test_res(:,8)<0.1);
    %     test_res(findFalseStarts,[10 11 12])=NaN;
    %     warning('correcting for false starts')
    for nbt=1:2
        tp_nogos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,11)),11);
        tp_gos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,12)),12);
        [dprime_test(n,nbt), crit_test(n,nbt)]=calc_dprime2((tp_gos==1),(tp_nogos==0));
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
            if npr==1
                last_pr_tridx=0;
            else
                last_pr_tridx=these_probes(npr-1,6);
            end
            %             last_pr_tridx=20;
            number_STD(n,nbl,npr)=sum(these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,5)~=3);
            number_DEV(n,nbl,npr)=sum(these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,5)==3);
            
            temp_testres=these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,:);
            %             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            tcorr_go=nanmean(temp_testres(:,12));%/corr_go(n,these_probes(npr,5));
            tcorr_nogo=nanmean(temp_testres(:,11));%/corr_nogo(n,these_probes(npr,5));
            num_go=sum(~isnan(temp_testres(:,12)));
            num_nogo=sum(~isnan(temp_testres(:,11)));
            rt_go=nanmean(temp_testres(~isnan(temp_testres(:,12)),10)-temp_testres(~isnan(temp_testres(:,12)),8));
            rt_nogo=nanmean(temp_testres(~isnan(temp_testres(:,11)),10)-temp_testres(~isnan(temp_testres(:,11)),8));
            
            tdp=calc_dprime2((temp_testres(~isnan(temp_testres(:,12)),12)==1),(temp_testres(~isnan(temp_testres(:,11)),11)==0));
            all_probes_mat=[all_probes_mat ; [n nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_go tcorr_nogo num_go num_nogo tdp (these_probes(npr,37)) rt_go rt_nogo]];
            
            temp_testres=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,12)),:);
            tcorr_go=(temp_testres(end-15:end,12))';%/corr_go(n,these_probes(npr,5));
            rtcorr_go=(temp_testres(~isnan(temp_testres(end-15:end,12)),10)-temp_testres(~isnan(temp_testres(end-15:end,12)),8))';%/corr_go(n,these_probes(npr,5));
            temp_testres=these_trials(these_trials(:,4)<this_pr_tridx & ~isnan(these_trials(:,11)),:);
            tcorr_nogo=(temp_testres(end-1:end,11))';%/corr_go(n,these_probes(npr,5));
            rtcorr_nogo=(temp_testres(~isnan(temp_testres(end-1:end,11)),10)-temp_testres(~isnan(temp_testres(end-1:end,11)),8))';%/corr_go(n,these_probes(npr,5));
            all_probes_mat_fullgo=[all_probes_mat_fullgo ; [n nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_go]];
            all_probes_mat_fullnogo=[all_probes_mat_fullnogo ; [n nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_nogo]];
            all_probes_mat_fullrtgo=[all_probes_mat_fullrtgo ; [n nbl these_probes(npr,5) these_probes(npr,32) npr rtcorr_go]];
            all_probes_mat_fullrtnogo=[all_probes_mat_fullrtnogo ; [n nbl these_probes(npr,5) these_probes(npr,32) npr rtcorr_nogo]];
            
            [dprime, crit]=calc_dprime2(tcorr_go==1,tcorr_nogo==0);
            all_probes_mat2=[all_probes_mat2 ; [n nbl these_probes(npr,5) these_probes(npr,32) npr nanmean(tcorr_go) nanmean(tcorr_nogo) num_go num_nogo tdp (these_probes(npr,37)) rt_go rt_nogo dprime crit]];
        end
    end
end

%%
forANOVA=[];
forANOVAgr=[];
forANOVA2=[];
for n=1:length(files)
    for ntask=1:2
        for nSta=1:3
            this_sub_go=all_probes_mat_fullgo(all_probes_mat_fullgo(:,1)==n & all_probes_mat_fullgo(:,3)==ntask & all_probes_mat_fullgo(:,4)==nSta,6:end);
            this_sub_nogo=all_probes_mat_fullnogo(all_probes_mat_fullnogo(:,1)==n & all_probes_mat_fullnogo(:,3)==ntask & all_probes_mat_fullnogo(:,4)==nSta,6:end);
            this_sub_go=reshape(this_sub_go,1,numel(this_sub_go));
            this_sub_nogo=reshape(this_sub_nogo,1,numel(this_sub_nogo));
            
            this_sub_rtgo=all_probes_mat_fullrtgo(all_probes_mat_fullrtgo(:,1)==n & all_probes_mat_fullrtgo(:,3)==ntask & all_probes_mat_fullrtgo(:,4)==nSta,6:end);
            this_sub_rtnogo=all_probes_mat_fullrtnogo(all_probes_mat_fullrtnogo(:,1)==n & all_probes_mat_fullrtnogo(:,3)==ntask & all_probes_mat_fullrtnogo(:,4)==nSta,6:end);
            this_sub_rtgo=reshape(this_sub_rtgo,1,numel(this_sub_rtgo));
            this_sub_rtnogo=reshape(this_sub_rtnogo,1,numel(this_sub_rtnogo));
            
            numdata_bysubj_state(n,ntask,nSta,:)=[length(this_sub_go) length(this_sub_nogo)];
            if length(this_sub_nogo)>=4
                [dprime_bysubj_state(n,ntask,nSta), crit_bysubj_state(n,ntask,nSta)] = calc_dprime2(this_sub_go,this_sub_nogo==0);
                forANOVA=[forANOVA ; [dprime_bysubj_state(n,ntask,nSta), crit_bysubj_state(n,ntask,nSta)]];
                forANOVAgr=[forANOVAgr ; [n ntask nSta]];
                
                RTGO_bysubj_state(n,ntask,nSta)=nanmean(this_sub_rtgo);
                RTNOGO_bysubj_state(n,ntask,nSta)=nanmean(this_sub_rtnogo);
                forANOVA2=[forANOVA2 ; [RTGO_bysubj_state(n,ntask,nSta), RTNOGO_bysubj_state(n,ntask,nSta)]];
                
            else
                dprime_bysubj_state(n,ntask,nSta)=NaN;
                crit_bysubj_state(n,ntask,nSta)=NaN;
                
                RTGO_bysubj_state(n,ntask,nSta)=NaN;
                RTNOGO_bysubj_state(n,ntask,nSta)=NaN;
            end
        end
    end
end

%%
tbl=array2table([forANOVA forANOVA2 forANOVAgr],'VariableNames',{'dp','crit','RTGO','RTNOGO','SubID','Task','State'});
tbl.SubID=categorical(tbl.SubID);
tbl.Task=categorical(tbl.Task);
tbl.State=categorical(tbl.State);


mdl_0= fitlme(tbl,'dp~1+(1|SubID)');
mdl_1= fitlme(tbl,'dp~Task+(1|SubID)');
mdl_2= fitlme(tbl,'dp~Task+State+(1|SubID)');
mdl_3= fitlme(tbl,'dp~Task*State+(1|SubID)');

mdl2_0= fitlme(tbl,'crit~1+(1|SubID)');
mdl2_1= fitlme(tbl,'crit~Task+(1|SubID)');
mdl2_2= fitlme(tbl,'crit~Task+State+(1|SubID)');
mdl2_3= fitlme(tbl,'crit~Task*State+(1|SubID)');

mdl3_0= fitlme(tbl,'RTGO~1+(1|SubID)');
mdl3_1= fitlme(tbl,'RTGO~Task+(1|SubID)');
mdl3_2= fitlme(tbl,'RTGO~Task+State+(1|SubID)');
mdl3_3= fitlme(tbl,'RTGO~Task*State+(1|SubID)');

mdl4_0= fitlme(tbl,'RTNOGO~1+(1|SubID)');
mdl4_1= fitlme(tbl,'RTNOGO~Task+(1|SubID)');
mdl4_2= fitlme(tbl,'RTNOGO~Task+State+(1|SubID)');
mdl4_3= fitlme(tbl,'RTNOGO~Task*State+(1|SubID)');

tbl.State=reordercats(tbl.State,{'2','1','3'});
mdlb_0= fitlme(tbl,'dp~1+(1|SubID)');
mdlb_1= fitlme(tbl,'dp~Task+(1|SubID)');
mdlb_2= fitlme(tbl,'dp~Task+State+(1|SubID)');
mdlb_3= fitlme(tbl,'dp~Task*State+(1|SubID)');

mdl2b_0= fitlme(tbl,'crit~1+(1|SubID)');
mdl2b_1= fitlme(tbl,'crit~Task+(1|SubID)');
mdl2b_2= fitlme(tbl,'crit~Task+State+(1|SubID)');
mdl2b_3= fitlme(tbl,'crit~Task*State+(1|SubID)');

%%
figure;
subplot(1,2,1); format_fig;
for nSta=1:3
    simpleBarPlot(nSta-0.2,squeeze(dprime_bysubj_state(:,1,nSta)),Colors(nSta,:),0.35,'k',[],3);
    simpleBarPlot(nSta+0.2,squeeze(dprime_bysubj_state(:,2,nSta)),[1 1 1;Colors(nSta,:)],0.35,'k',[],3);
end
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
ylabel('d''')

subplot(1,2,2); format_fig;
for nSta=1:3
    simpleBarPlot(nSta-0.2,squeeze(crit_bysubj_state(:,1,nSta)),Colors(nSta,:),0.35,'k',[],3);
    simpleBarPlot(nSta+0.2,squeeze(crit_bysubj_state(:,2,nSta)),[1 1 1;Colors(nSta,:)],0.35,'k',[],3);
end
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
ylabel('crit')

% figure;
% subplot(1,2,1); format_fig;
% for nSta=1:3
%     simpleBarPlot(nSta-0.2,squeeze(RTGO_bysubj_state(:,1,nSta)),Colors(nSta,:),0.35,'k',[],3);
%     simpleBarPlot(nSta+0.2,squeeze(RTGO_bysubj_state(:,2,nSta)),[1 1 1;Colors(nSta,:)],0.35,'k',[],3);
% end
% set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
%
% subplot(1,2,2); format_fig;
% for nSta=1:3
%     simpleBarPlot(nSta-0.2,squeeze(RTNOGO_bysubj_state(:,1,nSta)),Colors(nSta,:),0.35,'k',[],3);
%     simpleBarPlot(nSta+0.2,squeeze(RTNOGO_bysubj_state(:,2,nSta)),[1 1 1;Colors(nSta,:)],0.35,'k',[],3);
% end
% set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})

%%
figure;
subplot(1,2,1); format_fig;
for nSta=1:3
    %     simpleBarPlot(nSta-0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==1,30)),Colors(nSta,:),0.35,'k',[],3);
    simple_PiratePlot(nSta-0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==1,14)),0.15,0.5,'p',Colors(nSta,:),'y',[],{'SizeData',36});
    %     simpleBarPlot(nSta+0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,30)),[1 1 1;Colors(nSta,:)],0.35,'k',[],3);
    simple_PiratePlot(nSta+0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,14)),0.15,0.5,'p',Colors(nSta,:),'y',[],{'SizeData',36});
end
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
ylabel('d''')
xlim([0.2 3.8])

subplot(1,2,2); format_fig;
for nSta=1:3
    %     simpleBarPlot(nSta-0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==1,31)),Colors(nSta,:),0.35,'k',[],3);
    simple_PiratePlot(nSta-0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==1,15)),0.15,0.5,'p',Colors(nSta,:),'y',[],{'SizeData',36});
    %     simpleBarPlot(nSta+0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,31)),[1 1 1;Colors(nSta,:)],0.35,'k',[],3);
    simple_PiratePlot(nSta+0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,15)),0.15,0.5,'p',Colors(nSta,:),'y',[],{'SizeData',36});
end
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
ylabel('crit')
xlim([0.2 3.8])


%% By subj - Norm on ON
figure;
subplot(2,2,1); format_fig;
for nSta=1:3
    simpleBarPlot(nSta-0.2,100*squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==1,6)),Colors(nSta,:),0.35,'k',[],3);
    simpleBarPlot(nSta+0.2,100*squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,6)),[1 1 1;Colors(nSta,:)],0.35,'k',[],3);
end
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
ylabel('Corr GO')
xlim([0.2 3.8])
ylim([90 100])

subplot(2,2,2); format_fig;
for nSta=1:3
    simpleBarPlot(nSta-0.2,100*squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==1,7)),Colors(nSta,:),0.35,'k',[],3);
    simpleBarPlot(nSta+0.2,100*squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,7)),[1 1 1;Colors(nSta,:)],0.35,'k',[],3);
end
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
ylabel('Corr NOGO')
xlim([0.2 3.8])
ylim([40 90])

subplot(2,2,3); format_fig;
for nSta=1:3
    simpleBarPlot(nSta-0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==1,14)),Colors(nSta,:),0.35,'k',[],3);
    simpleBarPlot(nSta+0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,14)),[1 1 1;Colors(nSta,:)],0.35,'k',[],3);
    %      simple_PiratePlot(nSta+0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,31)),0.15,0.5,'p',Colors(nSta,:),'y',[],{'SizeData',36});
end
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
ylabel('d''')
xlim([0.2 3.8])

subplot(2,2,4); format_fig;
for nSta=1:3
    simpleBarPlot(nSta-0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==1,15)),Colors(nSta,:),0.35,'k',[],3);
    simpleBarPlot(nSta+0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,15)),[1 1 1;Colors(nSta,:)],0.35,'k',[],3);
    %      simple_PiratePlot(nSta+0.2,squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,31)),0.15,0.5,'p',Colors(nSta,:),'y',[],{'SizeData',36});
end
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
ylabel('norm-crit')
xlim([0.2 3.8])

%% By subj - Norm on ON
figure; set(gcf,'Position',[ 680     1   604   804])
plotNames={'GO','NOGO','d''','crit'};
plotCol=[6 7 14 15];
plotCorr=[100 100 1 1];
plotLims=[92 100; 40 90; 0 3; -1 -0.4];
for nplot=1:4
    subplot(4,2,2*nplot-1); format_fig;
    for nSta=1:3
        simpleBarPlot(nSta,plotCorr(nplot)*squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==1,plotCol(nplot))),Colors(nSta,:),0.9,'k',[],3);
        
    end
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    ylabel(plotNames{nplot})
    xlim([0.2 3.8])
    ylim(plotLims(nplot,:))
    
    subplot(4,2,2*nplot); format_fig;
    for nSta=1:3
        simpleBarPlot(nSta,plotCorr(nplot)*squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,plotCol(nplot))),[1 1 1;Colors(nSta,:)],0.9,'k',[],3);
    end
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    ylabel(plotNames{nplot})
    xlim([0.2 3.8])
    ylim(plotLims(nplot,:))
end

%%
%% By subj - Norm on ON
figure; set(gcf,'Position',[ 1         608        1239         370])
plotNames={'GO','NOGO','d''','crit','RT'};
plotCol=[6 7 14 15 12];
plotCorr=[100 100 1 1 1];
plotLims=[92 100; 45 85; 1 3;-1 -0.5; 0.45 0.65];
for nplot=1:5
%     subplot(5,2,2*nplot-1); format_fig; hold on;
    subplot(2,5,nplot); format_fig; hold on;
    mean_temp=[]; sem_temp=[];
    for nSta=1:3
        temp=plotCorr(nplot)*squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==1,plotCol(nplot)));
        mean_temp(nSta)=mean(temp);
        sem_temp(nSta)=sem(temp);
        line(nSta*[1 1],[-1 1]*sem_temp(nSta)+mean_temp(nSta),'Color',Colors(nSta,:),'LineWidth',3)
    end
    plot(1:3,mean_temp,'Color',[1 1 1]*0.5,'LineWidth',2)
    for nSta=1:3
        scatter(nSta,mean_temp(nSta),'MarkerFaceColor',Colors(nSta,:),'MarkerEdgeColor',Colors(nSta,:),'SizeData',72,'LineWidth',3);
    end
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    ylabel(plotNames{nplot})
    xlim([0.2 3.8])
    ylim(plotLims(nplot,:))
    
%     subplot(5,2,2*nplot); format_fig; hold on;
    subplot(2,5,5+nplot); format_fig; hold on;
   mean_temp=[]; sem_temp=[];
    for nSta=1:3
        temp=plotCorr(nplot)*squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,plotCol(nplot)));
        mean_temp(nSta)=mean(temp);
        sem_temp(nSta)=sem(temp);
        line(nSta*[1 1],[-1 1]*sem_temp(nSta)+mean_temp(nSta),'Color',Colors(nSta,:),'LineWidth',3)
    end
    plot(1:3,mean_temp,'Color',[1 1 1]*0.5,'LineWidth',2)
    for nSta=1:3
        scatter(nSta,mean_temp(nSta),'Marker','s','MarkerFaceColor',Colors(nSta,:),'MarkerEdgeColor',Colors(nSta,:),'SizeData',72,'LineWidth',3);
    end
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    ylabel(plotNames{nplot})
    xlim([0.2 3.8])
    ylim(plotLims(nplot,:))
end

%% Try to predict On vs Off
% n nbl these_probes(npr,5) these_probes(npr,32) npr nanmean(tcorr_go) nanmean(tcorr_nogo) num_go num_nogo tdp (these_probes(npr,37)) rt_go rt_nogo dprime crit
OnOFF=double(all_probes_mat2(:,4)==1);
OnOFF(all_probes_mat2(:,4)==4)=NaN;
MWMB=nan(size(OnOFF,1),1);
MWMB(all_probes_mat2(:,4)==2)=1;
MWMB(all_probes_mat2(:,4)==3)=0;
tbl=array2table([all_probes_mat2 OnOFF MWMB],'VariableNames',{'SubID','nBlock','Task','State','nProbe','CorrGo','CorrNoGo','nGo','nNoGo','dp2','perfrating','RTGO','RTNOGO','dp','crit','ONOFF','MWMB'});
tbl.SubID=categorical(tbl.SubID);
tbl.Task=categorical(tbl.Task);
tbl.State=categorical(tbl.State);
writetable(tbl,'/Users/tand0009/Data/WanderIM/behav/wanderIM_binary_table_ONOFF_MWMB.txt')

mdl_ONOFF_full= fitglme(tbl,'ONOFF~nBlock+nProbe+(CorrGo+CorrNoGo+dp+crit+RTGO)+(1|SubID)','Distribution','binomial');
mdl_ONOFF_full2= fitglme(tbl,'ONOFF~nBlock+nProbe+Task*(CorrGo+CorrNoGo+dp+crit+RTGO)+(1|SubID)','Distribution','binomial');

mdl_MWMB_full= fitglme(tbl,'MWMB~nBlock+nProbe+CorrGo+CorrNoGo+dp+crit+RTGO+(1|SubID)','Distribution','binomial');
mdl_MWMB_full2= fitglme(tbl,'MWMB~nBlock+nProbe+Task*(CorrGo+CorrNoGo+dp+crit+RTGO)+(1|SubID)','Distribution','binomial');