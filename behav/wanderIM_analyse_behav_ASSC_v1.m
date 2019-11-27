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
all_probes_mat3=[];
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
            %                         last_pr_tridx=this_pr_tridx-20;
            thesegotrials=find(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx & these_trials(:,5)~=3);
            thesenogotrials=find(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx & these_trials(:,5)==3);
            number_STD(n,nbl,npr)=length(thesegotrials(end-17:end));
            number_DEV(n,nbl,npr)=length(thesenogotrials(end-1:end));
            
            temp_testres=these_trials([thesegotrials(end-17:end) ; thesenogotrials(end-1:end)],:);
            %             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            tcorr_go=nanmean(temp_testres(:,12));%/corr_go(n,these_probes(npr,5));
            tcorr_nogo=nanmean(temp_testres(:,11));%/corr_nogo(n,these_probes(npr,5));
            num_go=sum(~isnan(temp_testres(:,12)));
            num_nogo=sum(~isnan(temp_testres(:,11)));
            rt_go=nanmean(temp_testres(~isnan(temp_testres(:,12)),10)-temp_testres(~isnan(temp_testres(:,12)),8));
            rt_nogo=nanmean(temp_testres(~isnan(temp_testres(:,11)),10)-temp_testres(~isnan(temp_testres(:,11)),8));
            
            tdp=calc_dprime2((temp_testres(~isnan(temp_testres(:,12)),12)==1),(temp_testres(~isnan(temp_testres(:,11)),11)==0));
            all_probes_mat=[all_probes_mat ; [n nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_go tcorr_nogo num_go num_nogo tdp (these_probes(npr,37)) rt_go rt_nogo]];
            
            temp_testres=these_trials([thesegotrials(end-17:end) ; thesenogotrials(end-1:end)],:);
            tcorr_go=(temp_testres(~isnan(temp_testres(:,12)),12))';%/corr_go(n,these_probes(npr,5));
            rtcorr_go=(temp_testres(~isnan(temp_testres(:,12)),10)-temp_testres(~isnan(temp_testres(:,12)),8))';%/corr_go(n,these_probes(npr,5));
            tcorr_nogo=(temp_testres(~isnan(temp_testres(:,11)),11))';%/corr_go(n,these_probes(npr,5));
            rtcorr_nogo=(temp_testres(~isnan(temp_testres(:,11)),10)-temp_testres(~isnan(temp_testres(:,11)),8))';%/corr_go(n,these_probes(npr,5));
            all_probes_mat_fullgo=[all_probes_mat_fullgo ; [n nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr tcorr_go]];
            all_probes_mat_fullnogo=[all_probes_mat_fullnogo ; [n nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr tcorr_nogo]];
            all_probes_mat_fullrtgo=[all_probes_mat_fullrtgo ; [n nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr rtcorr_go]];
            all_probes_mat_fullrtnogo=[all_probes_mat_fullrtnogo ; [n nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr rtcorr_nogo]];
            
            [dprime, crit]=calc_dprime2(tcorr_go==1,tcorr_nogo==0);
            all_probes_mat2=[all_probes_mat2 ; [n nbl these_probes(npr,5) these_probes(npr,32) npr nanmean(tcorr_go) nanmean(tcorr_nogo) num_go num_nogo tdp (these_probes(npr,38)) rt_go rt_nogo dprime crit]];
            all_probes_mat3=[all_probes_mat3 ; [n nbl npr these_probes(npr,5) these_probes(npr,32) npr nanmean(tcorr_go) nanmean(tcorr_nogo) num_go num_nogo tdp (these_probes(npr,38)) rt_go rt_nogo dprime crit these_probes(npr,31:38)]];
        end
    end
end

%%
allpption_alongtask=[];
allpption_alongtask2=[];
forANOVA=[];
forANOVAgr=[];
forANOVA2=[];
for n=1:length(files)
    for ntask=1:2
        for nSta=1:3
            %             if nSta==2
            %                 this_sub_go=all_probes_mat_fullgo(all_probes_mat_fullgo(:,1)==n & all_probes_mat_fullgo(:,3)==ntask & all_probes_mat_fullgo(:,4)==nSta & all_probes_mat_fullgo(:,5)==nSta,end);
            %                 this_sub_nogo=all_probes_mat_fullnogo(all_probes_mat_fullnogo(:,1)==n & all_probes_mat_fullnogo(:,3)==ntask & all_probes_mat_fullnogo(:,4)==nSta & all_probes_mat_fullgo(:,5)==nSta,end);
            %                 this_sub_go=reshape(this_sub_go,1,numel(this_sub_go));
            %                 this_sub_nogo=reshape(this_sub_nogo,1,numel(this_sub_nogo));
            %
            %                 this_sub_rtgo=all_probes_mat_fullrtgo(all_probes_mat_fullrtgo(:,1)==n & all_probes_mat_fullrtgo(:,3)==ntask & all_probes_mat_fullrtgo(:,4)==nSta & all_probes_mat_fullgo(:,5)==nSta,end);
            %                 this_sub_rtnogo=all_probes_mat_fullrtnogo(all_probes_mat_fullrtnogo(:,1)==n & all_probes_mat_fullrtnogo(:,3)==ntask & all_probes_mat_fullrtnogo(:,4)==nSta & all_probes_mat_fullgo(:,5)==nSta,end);
            %                 this_sub_rtgo=reshape(this_sub_rtgo,1,numel(this_sub_rtgo));
            %                 this_sub_rtnogo=reshape(this_sub_rtnogo,1,numel(this_sub_rtnogo));
            %             else
            this_sub_go=all_probes_mat_fullgo(all_probes_mat_fullgo(:,1)==n & all_probes_mat_fullgo(:,3)==ntask & all_probes_mat_fullgo(:,4)==nSta,end);
            this_sub_nogo=all_probes_mat_fullnogo(all_probes_mat_fullnogo(:,1)==n & all_probes_mat_fullnogo(:,3)==ntask & all_probes_mat_fullnogo(:,4)==nSta,end);
            this_sub_go=reshape(this_sub_go,1,numel(this_sub_go));
            this_sub_nogo=reshape(this_sub_nogo,1,numel(this_sub_nogo));
            
            this_sub_rtgo=all_probes_mat_fullrtgo(all_probes_mat_fullrtgo(:,1)==n & all_probes_mat_fullrtgo(:,3)==ntask & all_probes_mat_fullrtgo(:,4)==nSta,end);
            this_sub_rtnogo=all_probes_mat_fullrtnogo(all_probes_mat_fullrtnogo(:,1)==n & all_probes_mat_fullrtnogo(:,3)==ntask & all_probes_mat_fullrtnogo(:,4)==nSta,end);
            this_sub_rtgo=reshape(this_sub_rtgo,1,numel(this_sub_rtgo));
            this_sub_rtnogo=reshape(this_sub_rtnogo,1,numel(this_sub_rtnogo));
            %             end
            
            numdata_bysubj_state(n,ntask,nSta,:)=[length(this_sub_go) length(this_sub_nogo)];
            if length(this_sub_nogo)>=4
                [dprime_bysubj_state(n,ntask,nSta), crit_bysubj_state(n,ntask,nSta)] = calc_dprime2(this_sub_go,this_sub_nogo==0);
                go_bysubj_state(n,ntask,nSta)=nanmean(this_sub_go);
                nogo_bysubj_state(n,ntask,nSta)=nanmean(this_sub_nogo);
                
                RTGO_bysubj_state(n,ntask,nSta)=nanmedian(this_sub_rtgo);
                RTNOGO_bysubj_state(n,ntask,nSta)=nanmedian(this_sub_rtnogo);
                
                forANOVA=[forANOVA ; [dprime_bysubj_state(n,ntask,nSta), crit_bysubj_state(n,ntask,nSta)]];
                forANOVAgr=[forANOVAgr ; [n ntask nSta]];
                forANOVA2=[forANOVA2 ; [RTGO_bysubj_state(n,ntask,nSta), RTNOGO_bysubj_state(n,ntask,nSta)]];
            else
                %                 disp('missing')
                dprime_bysubj_state(n,ntask,nSta)=NaN;
                crit_bysubj_state(n,ntask,nSta)=NaN;
                go_bysubj_state(n,ntask,nSta)=NaN;
                nogo_bysubj_state(n,ntask,nSta)=NaN;
                RTGO_bysubj_state(n,ntask,nSta)=NaN;
                RTNOGO_bysubj_state(n,ntask,nSta)=NaN;
            end
        end
        
        % phenomenon characteristics
        for nSta=1:3
            this_sub_charact=all_probes_mat3(all_probes_mat3(:,1)==n & all_probes_mat3(:,4)==ntask & all_probes_mat3(:,5)==nSta,end-7:end);
            probes_quest(n,ntask,nSta,:)=nanmean(this_sub_charact,1);
            probes_quest_n(n,ntask,nSta)=size(this_sub_charact,1);
        end
        % phenomenon characteristics
        for nSta2=1:3
            %             if nSta2==2
            %             this_sub_charact=all_probes_mat3(all_probes_mat3(:,1)==n & all_probes_mat3(:,4)==ntask & all_probes_mat3(:,5)==2 & all_probes_mat3(:,19)==nSta2,end-7:end);
            %             else
            this_sub_charact=all_probes_mat3(all_probes_mat3(:,1)==n & all_probes_mat3(:,4)==ntask & all_probes_mat3(:,5)==nSta2,end-7:end);
            %             end
            probes_quest2(n,ntask,nSta2,:)=nanmean(this_sub_charact,1);
            probes_quest2_n(n,ntask,nSta2)=size(this_sub_charact,1);
        end
    end
    % along task
    partidx=[1:4;3:6;5:8;7:10];
    for nblock=1:6
        this_sub_part=all_probes_mat3(all_probes_mat3(:,1)==n & all_probes_mat3(:,2)==nblock,end-6);
        this_sub_part2=all_probes_mat3(all_probes_mat3(:,1)==n & all_probes_mat3(:,2)==nblock,end);
        this_sub_part3=all_probes_mat3(all_probes_mat3(:,1)==n & all_probes_mat3(:,2)==nblock,19);
        for nstate=1:4
            %                 if nstate~=2
            pption_alongtask(n,nblock,nstate)=mean(this_sub_part==nstate);
            %                 else
            %                 pption_alongtask(n,nblock,nstate)=mean(this_sub_part==nstate & this_sub_part3==nstate);
            %                 end
            allpption_alongtask=[allpption_alongtask ; [n nblock nstate mean(this_sub_part==nstate)]];
            allpption_alongtask2=[allpption_alongtask2 ; [n nblock nstate mean(this_sub_part2)]];
        end
    end
end

%% Proportion of probes and characteristics
%%% Pption
thisColors=[Colors ; [1 1 1]*0.4];
figure;  set(gcf,'Position',[-29        1275        1218/2         412]);
% subplot(1,2,1);
format_fig;
all_pption_MS=pption_MS; %cat(3,pption_MS(:,:,[1]),pption_MS(:,:,[2]).*pption_Ori(:,:,2),pption_MS(:,:,[3]));
for nst=1:3
    hb(1)=simpleBarPlot(nst-0.2,100*squeeze(all_pption_MS(:,1,nst)),Colors(nst,:),0.35,'k');
    hb(2)=simpleBarPlot(nst+0.2,100*squeeze(all_pption_MS(:,2,nst)),[1 1 1 ; Colors(nst,:)],0.35,'k');
end
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
ylabel('% probes')
title('Mind-State')
xlim([0.2 3.8])
% legend(hb,{'Face','Square'})
set(gca,'fontsize',20)
%
% subplot(1,2,2); format_fig;
% for nst=1:3
%     simpleBarPlot(nst-0.2,100*squeeze(pption_Ori(:,1,nst)),Colors5(nst+3,:),0.35,'k');
%     simpleBarPlot(nst+0.2,100*squeeze(pption_Ori(:,2,nst)),[1 1 1 ; Colors5(nst+3,:)],0.35,'k');
% end
% set(gca,'XTick',1:3,'XTickLabel',{'environ^n^t','personal','task'})
% ylabel('% probes')
% title('Origin of OFF-state')
% xlim([0.2 3.8])
% set(gca,'fontsize',20)

%% Spider
figure; hold on; format_fig;
set(gcf,'Position',[100   100   420   420]);
set(gca,'XColor','w','YColor','w');
for nstep=1:4
    viscircles([0 0], nstep,'Color','k','LineStyle','--','LineWidth',1);
end
myQ=[4 5 6 7 8];
for nQ=1:length(myQ)
    for nstate=1:3
        if ismember(myQ(nQ),[4 5 8])
            thisQ=5-squeeze(nanmean(probes_quest(:,:,nstate,myQ(nQ)),2));
        else
            thisQ=squeeze(nanmean(probes_quest(:,:,nstate,myQ(nQ)),2));
        end
        toplot(nstate,nQ)=nanmean(thisQ);
        tostat{nstate,nQ}=(thisQ);
    end
end
% plot the axis
th = (2*pi/length(myQ))*(ones(2,1)*(length(myQ):-1:1));
% Axis start and end
r = [0;4]*ones(1,length(myQ));
% Conversion to cartesian coordinates to plot using regular plot.
[x,y] = pol2cart(th, r);
hLine = line(x, y,...
    'LineWidth', 1.5,...
    'Color', [1, 1, 1]*0.5  );
for nstate=1:3
    % Radius
    R = toplot(nstate,:); % - (1)./((4-1)) + 0.1;
    R = [R'; R(1)];
    Th = (2*pi/length(myQ)) * ((length(myQ):-1:0)'*ones(1,1));
    % polar(Th, R)
    [X, Y] = pol2cart(Th, R);
    h = plot(X, Y,'Color',thisColors(nstate,:),'LineWidth',4);
    
end
for nQ=1:5
    pVano(nQ)=anova1([tostat{1,nQ} tostat{2,nQ} tostat{3,nQ}],[],'off');
    [~, pVpost(nQ)]=ttest(tostat{2,nQ},tostat{3,nQ});
end
% %  add labels
% minV = 0;
% maxV = 4;
% mylabels={'Awareness','Intention','Engagement','Performance','Alertness'};
% for j = 1:length(myQ)
%     % Generate the axis label
%     msg = mylabels{j};
%     [mx, my] = pol2cart( th(1, j), 4.1);
%     text(mx, my, msg);
% end
%% Dynamics withi task
figure;
format_fig; set(gcf,'Position',[440   321   694   477])
for nstate=1:4
    templot=100*squeeze(pption_alongtask(:,:,nstate));
    errorbar(1:size(templot,2),nanmean(templot,1),sem(templot),'Color',thisColors(nstate,:),'LineWidth',2);
    hold on;
end
xlabel('Block #');
title('Across block dynamics')
ylabel('% of probes')
format_fig;
xlim([0.5 6.5])
ylim([0 70])


%% ON/MW/MB
myS=unique(all_probes_mat2(:,1));
% figure; set(gcf,'Position',[ 1         608        1239         370])
% plotNames={'GO','NOGO','d''','crit','RT'};
% plotCol=[6 7 14 15 12];
% plotCorr=[100 100 1 1 1];
% plotLims=[92 100; 45 85; 1 3;-0.9 -0.6; 0.55 0.65; 0.48 0.55];
figure; set(gcf,'Position',[ 1         608        831         370])
plotNames={'d''','crit','RT'};
plotCol=[14 15 12];
plotCorr=[1 1 1];
plotLims=[1 3;-0.95 -0.6; 0.55 0.65; 0.48 0.55];
for nplot=[3 2 1]
    %     subplot(5,2,2*nplot-1); format_fig; hold on;
    subplot(2,3,nplot); format_fig; hold on;
    mean_temp=[]; sem_temp=[];
    for nSta=1:3
        temp=plotCorr(nplot)*squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==1,plotCol(nplot)));
        tempS=squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==1,1));
        tempbyS=[];
        for nS=1:length(myS)
            tempbyS(nS,1)=nanmean(temp(tempS==myS(nS)));
            tempbyS(nS,2)=nansum(tempS==myS(nS));
        end
        mean_temp(nSta)=nansum(tempbyS(:,1).*tempbyS(:,2))./sum(tempbyS(:,2));
        sem_temp(nSta)=std(tempbyS(~isnan(tempbyS(:,1)),1),tempbyS(~isnan(tempbyS(:,1)),2))/sqrt(sum(~isnan(tempbyS(:,1)))-1);
        line(nSta*[1 1],[-1 1]*sem_temp(nSta)+mean_temp(nSta),'Color',Colors(nSta,:),'LineWidth',3)
    end
    plot(1:3,mean_temp,'Color',[1 1 1]*0.5,'LineWidth',2)
    for nSta=1:3
        scatter(nSta,mean_temp(nSta),'MarkerFaceColor',Colors(nSta,:),'MarkerEdgeColor',Colors(nSta,:),'SizeData',72,'LineWidth',3);
    end
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    title(plotNames{nplot})
    xlim([0.2 3.8])
    ylim(plotLims(nplot,:))
    set(gca,'fontsize',20)
    
    %     subplot(5,2,2*nplot); format_fig; hold on;
    subplot(2,3,3+nplot); format_fig; hold on;
    mean_temp=[]; sem_temp=[];
    for nSta=1:3
        temp=plotCorr(nplot)*squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,plotCol(nplot)));
        tempS=squeeze(all_probes_mat2(all_probes_mat2(:,4)==nSta & all_probes_mat2(:,3)==2,1));
        tempbyS=[];
        for nS=1:length(myS)
            tempbyS(nS,1)=nanmean(temp(tempS==myS(nS)));
            tempbyS(nS,2)=nansum(tempS==myS(nS));
        end
        mean_temp(nSta)=nansum(tempbyS(:,1).*tempbyS(:,2))./sum(tempbyS(:,2));
        sem_temp(nSta)=std(tempbyS(~isnan(tempbyS(:,1)),1),tempbyS(~isnan(tempbyS(:,1)),2))/sqrt(sum(~isnan(tempbyS(:,1)))-1);
        line(nSta*[1 1],[-1 1]*sem_temp(nSta)+mean_temp(nSta),'Color',Colors(nSta,:),'LineWidth',3)
    end
    plot(1:3,mean_temp,'Color',[1 1 1]*0.5,'LineWidth',2)
    for nSta=1:3
        scatter(nSta,mean_temp(nSta),'Marker','s','MarkerFaceColor',Colors(nSta,:),'MarkerEdgeColor',Colors(nSta,:),'SizeData',72,'LineWidth',3);
    end
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
    %     ylabel(plotNames{nplot})
    xlim([0.2 3.8])
    if nplot==3
        ylim(plotLims(nplot+1,:))
    else
        ylim(plotLims(nplot,:))
    end
    set(gca,'fontsize',20)
end

%%
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

%%
temp=all_probes_mat2(all_probes_mat2(:,4)<4,:);
tbl=array2table(temp,'VariableNames',{'SubID','nBl','Task','MS','nPr','GO','NOGO','nGO','nNOGO','tdp','vig','rtGO','rtNOGO','dp','crit'});
tbl.MS=categorical(tbl.MS);
tbl.Task=categorical(tbl.Task);
tbl.MS=reordercats(tbl.MS,[2 1 3]);

mdl_0= fitlme(tbl,'crit~1+(1|SubID)');
mdl_1= fitlme(tbl,'crit~MS+(1|SubID)'); %wining model
% mdl_2= fitlme(tbl,'crit~Task+MS+(1|SubID)');
% mdl_3= fitlme(tbl,'crit~Task*MS+(1|SubID)');

mdl2_0= fitlme(tbl,'rtGO~1+(1|SubID)');
mdl2_1= fitlme(tbl,'rtGO~MS+(1|SubID)');
mdl2_2= fitlme(tbl,'rtGO~Task+MS+(1|SubID)'); %wining model
% mdl2_3= fitlme(tbl,'rtGO~Task*MS+(1|SubID)');

mdl3_0= fitlme(tbl,'dp~1+(1|SubID)');
mdl3_1= fitlme(tbl,'dp~MS+(1|SubID)');
mdl3_2= fitlme(tbl,'dp~Task+MS+(1|SubID)'); %wining model
% mdl3_3= fitlme(tbl,'dp~Task*MS+(1|SubID)');

mdl4_0= fitlme(tbl,'GO~1+(1|SubID)');
mdl4_1= fitlme(tbl,'GO~MS+(1|SubID)');
mdl4_2= fitlme(tbl,'GO~Task+MS+(1|SubID)'); %wining model
% mdl4_3= fitlme(tbl,'dp~Task*MS+(1|SubID)');


mdl5_0= fitlme(tbl,'NOGO~1+(1|SubID)');
mdl5_1= fitlme(tbl,'NOGO~MS+(1|SubID)');
mdl5_2= fitlme(tbl,'NOGO~Task+MS+(1|SubID)'); %wining model
% mdl5_3= fitlme(tbl,'NOGO~Task*MS+(1|SubID)');

%%
mdl6_0= fitlme(tbl,'vig~1+(1|SubID)');
mdl6_1= fitlme(tbl,'vig~MS+(1|SubID)');
mdl6_2= fitlme(tbl,'vig~Task+MS+(1|SubID)'); %wining model
mdl6_3= fitlme(tbl,'vig~Task*MS+(1|SubID)'); %wining model
figure; format_fig;     hold on; set(gcf,'Position',[858   424   382   381])
myS=unique(all_probes_mat2(:,1));
for nT=1:2
    mean_temp=[];
    sem_temp=[];
    for nCond=1:3
        if nT==1
            tempC=Colors(nCond,:);
        else
            tempC=[ [1 1 1]; Colors(nCond,:)];
        end
        temp=all_probes_mat2(all_probes_mat2(:,3)==nT & all_probes_mat2(:,4)==nCond,11);
        tempS=squeeze(all_probes_mat2(all_probes_mat2(:,3)==nT & all_probes_mat2(:,4)==nCond,1));
        tempbyS=[];
        for nS=1:length(myS)
            tempbyS(nS,1)=nanmean(temp(tempS==myS(nS)));
            tempbyS(nS,2)=nansum(tempS==myS(nS));
        end
        mean_temp(nCond)=nansum(tempbyS(:,1).*tempbyS(:,2))./sum(tempbyS(:,2));
        sem_temp(nCond)=std(tempbyS(~isnan(tempbyS(:,1)),1),tempbyS(~isnan(tempbyS(:,1)),2))/sqrt(sum(~isnan(tempbyS(:,1)))-1);
        
        line(nCond*[1 1]+(2*nT-3)*0.1,[-1 1]*sem_temp(nCond)+mean_temp(nCond),'Color',Colors(nCond,:),'LineWidth',3)
    end
        plot((1:3)+(2*nT-3)*0.1,mean_temp,'Color',[1 1 1]*0.5,'LineWidth',2)
    for nCond=1:3
                if nT==1
                scatter(nCond+(2*nT-3)*0.1,mean_temp(nCond),'MarkerFaceColor',Colors(nCond,:),'MarkerEdgeColor',Colors(nCond,:),'SizeData',288,'LineWidth',3);
                else
                scatter(nCond+(2*nT-3)*0.1,mean_temp(nCond),'MarkerFaceColor','w','MarkerEdgeColor',Colors(nCond,:),'SizeData',288,'LineWidth',3);
                end
    end
    
end
set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'})
ylabel('Subj. Fatigue')
xlim([0.5 3.5])
ylim([1.5 3.5])
set(gca,'fontsize',22)