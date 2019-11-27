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
all_probes_mat=[];
all_probes_mat_ExGauss=[];
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
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
    
    %%% across blocks
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
            last_pr_tridx=this_pr_tridx-20;
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
            
            tempRT=temp_testres(:,10)-temp_testres(:,8);
            [X,fVal,exitFlag,solverOutput] = exgauss_fit(tempRT');
            
            
            [tdp tcrit]=calc_dprime((temp_testres(~isnan(temp_testres(:,12)),12)==1),(temp_testres(~isnan(temp_testres(:,11)),11)==0));
            all_probes_mat=[all_probes_mat ; [n nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_go tcorr_nogo num_go num_nogo tdp tcrit (these_probes(npr,38)) rt_go rt_nogo]];
            all_probes_mat_ExGauss=[all_probes_mat_ExGauss ; X];
        end
    end
end

%% Splitting trials based on tiredness
all_probes_mat(all_probes_mat(:,3)==1,:)=[];
uniSubs=unique(all_probes_mat(:,1));
normRT=[];
for nS=1:length(uniSubs)
    sub_res=all_probes_mat(all_probes_mat(:,1)==uniSubs(nS),:);
    sub_res2=all_probes_mat_ExGauss(all_probes_mat(:,1)==uniSubs(nS),:);
    normRT=[normRT ; zscore(sub_res(:,end-1))];
    fprintf('... ... alert %g - drowsy %g\n',sum(sub_res(:,12)<3),sum(sub_res(:,12)>=3))
    tired_tr=find(sub_res(:,12)>=3);
    alert_tr=find(sub_res(:,12)<3);
    if sum(sub_res(:,12)<3)<3 || sum(sub_res(:,12)>=3)<3
        dprime_vig(nS,:)=nan(1,2);
        crit_vig(nS,:)=nan(1,2);
        RTgo_vig(nS,:)=nan(1,2);
        CorrGo_vig(nS,:)=nan(1,2);
        CorrNoGo_vig(nS,:)=nan(1,2);
        
        ExGauss_vig(nS,:,:)=nan(1,3,2);
    else
        dprime_vig(nS,1)=mean(sub_res(alert_tr,end-4));
        dprime_vig(nS,2)=mean(sub_res(tired_tr,end-4));
        
        crit_vig(nS,1)=mean(sub_res(alert_tr,end-3));
        crit_vig(nS,2)=mean(sub_res(tired_tr,end-3));
        
        RTgo_vig(nS,1)=mean(sub_res(alert_tr,end-1));
        RTgo_vig(nS,2)=mean(sub_res(tired_tr,end-1));
        
        CorrGo_vig(nS,1)=mean(sub_res(alert_tr,end-8));
        CorrGo_vig(nS,2)=mean(sub_res(tired_tr,end-8));
        
        CorrNoGo_vig(nS,1)=mean(sub_res(alert_tr,end-7));
        CorrNoGo_vig(nS,2)=mean(sub_res(tired_tr,end-7));
        
        ExGauss_vig(nS,:,1)=nanmean(sub_res2(alert_tr,:));
        ExGauss_vig(nS,:,2)=nanmean(sub_res2(tired_tr,:));
    end
end

%%
figure; set(gcf,'Position',[680   330   799   648])
subplot(2,3,1)
format_fig;
simpleBarPlot(1,dprime_vig(:,1),'b',0.9,'k',[],3);
simpleBarPlot(2,dprime_vig(:,2),'r',0.9,'k',{0 dprime_vig(:,1) 0.05},3);
ylabel('d''')
xlim([0.2 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'Alert','Drowsy'})

subplot(2,3,2)
format_fig;
simpleBarPlot(1,crit_vig(:,1),'b',0.9,'k',[],3);
simpleBarPlot(2,crit_vig(:,2),'r',0.9,'k',{0 crit_vig(:,1) 0.05},3);
ylabel('crit')
xlim([0.2 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'Alert','Drowsy'})
set(gca,'YDir','reverse')

subplot(2,3,3)
format_fig;
simpleBarPlot(1,RTgo_vig(:,1),'b',0.9,'k',[],3);
simpleBarPlot(2,RTgo_vig(:,2),'r',0.9,'k',{0 RTgo_vig(:,1) 0.05},3);
ylabel('RT (s)')
xlim([0.2 2.8])
ylim([0.4 0.6])
set(gca,'XTick',1:2,'XTickLabel',{'Alert','Drowsy'})

subplot(2,3,4)
format_fig;
simpleBarPlot(1,100*CorrGo_vig(:,1),'b',0.9,'k',[],3);
simpleBarPlot(2,100*CorrGo_vig(:,2),'r',0.9,'k',{0 100*CorrGo_vig(:,1) 0.05},3);
ylabel('Hit Rate (%)')
xlim([0.2 2.8])
ylim([95 100])
set(gca,'XTick',1:2,'XTickLabel',{'Alert','Drowsy'})

subplot(2,3,5)
format_fig;
simpleBarPlot(1,100-100*CorrNoGo_vig(:,1),'b',0.9,'k',[],3);
simpleBarPlot(2,100-100*CorrNoGo_vig(:,2),'r',0.9,'k',{0 100*CorrNoGo_vig(:,1) 0.05},3);
ylabel('False Alarm (%)')
xlim([0.2 2.8])
ylim([10 45])
set(gca,'XTick',1:2,'XTickLabel',{'Alert','Drowsy'})

%%
figure; set(gcf,'Position',[680   330   799   648])
for nplot=1:3
    subplot(1,3,nplot)
    format_fig;
    simpleBarPlot(1,ExGauss_vig(:,nplot,1),'b',0.9,'k',[],3);
    simpleBarPlot(2,ExGauss_vig(:,nplot,2),'r',0.9,'k',{0 ExGauss_vig(:,nplot,1) 0.05},3);
%     ylabel('d''')
    xlim([0.2 2.8])
    set(gca,'XTick',1:2,'XTickLabel',{'Alert','Drowsy'})
end

%%
figure;
simpleCorPlotsetbin(all_probes_mat(:,12),all_probes_mat(:,end-1),1:4,{'o','k','k',288,3},[],0);
xlim([0.2 4.8])
ylim([0.45 0.6])
set(gca,'XTick',1:4);
xlabel('Tiredness scale')
ylabel('Mean GO RT (s)')
format_fig;
% figure;
% simpleCorPlotsetbin(all_probes_mat(:,12),normRT,1:4,[],[],0);