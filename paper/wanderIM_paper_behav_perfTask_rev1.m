%%
clear all
close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))
addpath(genpath(path_export))
addpath(genpath(path_RainCloudPlot))
addpath(path_export);
% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

% state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
% cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];

% load([pwd filesep 'paper_SubID'])
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
nSc=0;
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    probe_res(probe_res(:,32)==4,32)=3;
    SubID=SubjectInfo.subID;
    if ~ismember(SubID,GoodSubID)
        continue;
    end
    nSc=nSc+1;
    fprintf('... %s (%g)\n',SubID,nSc)
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
    %         findFalseStarts=find(test_res(:,10)-test_res(:,8)<0.1);
    %         test_res(findFalseStarts,[10 11 12])=NaN;
    % test_res(RTs<0.3,:)=[];
    %     warning('correcting for false starts')
    for nbt=1:2
        tp_nogos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,11)),11);
        tp_gos=test_res(test_res(:,2)==nbt & ~isnan(test_res(:,12)),12);
        [dprime_test(nSc,nbt), crit_test(nSc,nbt)]=calc_dprime2((tp_gos==1),(tp_nogos==0));
        corr_go(nSc,nbt)=nanmean(tp_gos);
        corr_nogo(nSc,nbt)=nanmean(tp_nogos);
        
        rt_gos(nSc,nbt)=nanmean(RTs(test_res(:,2)==nbt & ~isnan(test_res(:,12))));
        rt_nogos(nSc,nbt)=nanmean(RTs(test_res(:,2)==nbt & ~isnan(test_res(:,11))));
    end
    all_test_res=[all_test_res ; [nSc*ones(size(test_res,1),1) test_res(:,[1 2 4 11 12])]];
    
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
        pption_MS(nSc,nbt,:)=pption;
        pption_Ori(nSc,nbt,:)=pption2;
    end
    for nb=1:max(probe_res(:,3))
        temp=probe_res(probe_res(:,4)==nb,[5 31:38]);
        temp2=test_res(test_res(:,1)==nb,[11 12]);
        %             pption=[];
        %             for nstate=1:4
        %                 pption(nstate)=mean(temp(:,3)==nstate);
        %             end
        %             all_test_resprobes_perblock=[all_test_resprobes_perblock ; [nSc 1 unique(temp(:,1)) pption mean(temp(:,[5:8])) nanmean(temp2)]];
        all_test_resprobes_perblock=[all_test_resprobes_perblock ; [nSc*ones(size(temp,1),1) nb*ones(size(temp,1),1) temp]];
        all_probes_timing=[all_probes_timing ; diff(probe_res(probe_res(:,3)==nb,2))];
    end
    for nbt=1:2
        tp_probes=probe_res(probe_res(:,5)==nbt,31:38);
        mean_bytask_prf(nSc,nbt)=nanmean(tp_probes(:,7));
        
        
        for nstate=1:3
            mean_count_mindstates(nSc,nbt,nstate)=sum(tp_probes(:,2)==(nstate));
            mean_count_ori(nSc,nbt,nstate)=sum(tp_probes(:,3)==(nstate));
            
            if sum(tp_probes(:,2)==(nstate))~=0
                
                mean_byprobe_awa(nSc,nbt,nstate)=5-nanmean(tp_probes(tp_probes(:,2)==(nstate),4));
                mean_byprobe_wil(nSc,nbt,nstate)=5-nanmean(tp_probes(tp_probes(:,2)==(nstate),5));
                mean_byprobe_eng(nSc,nbt,nstate)=nanmean(tp_probes(tp_probes(:,2)==(nstate),6));
                mean_byprobe_prf(nSc,nbt,nstate)=nanmean(tp_probes(tp_probes(:,2)==(nstate),7));
                mean_byprobe_vig(nSc,nbt,nstate)=5-nanmean(tp_probes(tp_probes(:,2)==(nstate),8));
            else
                mean_byprobe_awa(nSc,nbt,nstate)=nan;
                mean_byprobe_wil(nSc,nbt,nstate)=nan;
                mean_byprobe_eng(nSc,nbt,nstate)=nan;
                mean_byprobe_prf(nSc,nbt,nstate)=nan;
                mean_byprobe_vig(nSc,nbt,nstate)=nan;
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
            number_STD(nSc,nbl,npr)=length(thesegotrials(end-17:end));
            number_DEV(nSc,nbl,npr)=length(thesenogotrials(end-1:end));
            
            temp_testres=these_trials([thesegotrials(end-17:end) ; thesenogotrials(end-1:end)],:);
            %             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            tcorr_go=nanmean(temp_testres(:,12));%/corr_go(nSc,these_probes(npr,5));
            tcorr_nogo=nanmean(temp_testres(:,11));%/corr_nogo(nSc,these_probes(npr,5));
            num_go=sum(~isnan(temp_testres(:,12)));
            num_nogo=sum(~isnan(temp_testres(:,11)));
            rt_go=(temp_testres(~isnan(temp_testres(:,12)),10)-temp_testres(~isnan(temp_testres(:,12)),8));
            rt_go=nanmean(rt_go(rt_go>0.3));
            rt_nogo=(temp_testres(~isnan(temp_testres(:,11)),10)-temp_testres(~isnan(temp_testres(:,11)),8));
            rt_nogo=nanmean(rt_nogo(rt_nogo>0.3));
            
            tdp=calc_dprime2((temp_testres(~isnan(temp_testres(:,12)),12)==1),(temp_testres(~isnan(temp_testres(:,11)),11)==0));
            all_probes_mat=[all_probes_mat ; [nSc nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_go tcorr_nogo num_go num_nogo tdp (these_probes(npr,37)) rt_go rt_nogo]];
            
            temp_testres=these_trials([thesegotrials(end-17:end) ; thesenogotrials(end-1:end)],:);
            tcorr_go=(temp_testres(~isnan(temp_testres(:,12)),12))';%/corr_go(nSc,these_probes(npr,5));
            rtcorr_go=(temp_testres(~isnan(temp_testres(:,12)),10)-temp_testres(~isnan(temp_testres(:,12)),8))';%/corr_go(nSc,these_probes(npr,5));
            tcorr_nogo=(temp_testres(~isnan(temp_testres(:,11)),11))';%/corr_go(nSc,these_probes(npr,5));
            rtcorr_nogo=(temp_testres(~isnan(temp_testres(:,11)),10)-temp_testres(~isnan(temp_testres(:,11)),8))';%/corr_go(nSc,these_probes(npr,5));
            all_probes_mat_fullgo=[all_probes_mat_fullgo ; [nSc nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr tcorr_go]];
            all_probes_mat_fullnogo=[all_probes_mat_fullnogo ; [nSc nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr tcorr_nogo]];
            all_probes_mat_fullrtgo=[all_probes_mat_fullrtgo ; [nSc nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr rtcorr_go]];
            all_probes_mat_fullrtnogo=[all_probes_mat_fullrtnogo ; [nSc nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr rtcorr_nogo]];
            
            [dprime, crit]=calc_dprime2(tcorr_go==1,tcorr_nogo==0);
            all_probes_mat2=[all_probes_mat2 ; [nSc nbl these_probes(npr,5) these_probes(npr,32) npr nanmean(tcorr_go) nanmean(tcorr_nogo) num_go num_nogo tdp (these_probes(npr,38)) rt_go rt_nogo dprime crit]];
            all_probes_mat3=[all_probes_mat3 ; [nSc nbl npr these_probes(npr,5) these_probes(npr,32) npr nanmean(tcorr_go) nanmean(tcorr_nogo) num_go num_nogo tdp (these_probes(npr,38)) rt_go rt_nogo dprime crit these_probes(npr,31:38)]];
        end
    end
end

%% Normalise by ON state
myS=unique(all_probes_mat2(:,1));
all_probes_mat2norm=all_probes_mat2;
% all_probes_mat2norm(all_probes_mat2norm(:,4)==4,4)=3;
for nS=1:length(myS)
    for nTask=1:2
        all_probes_mat2norm(all_probes_mat2norm(:,1)==myS(nS) & all_probes_mat2norm(:,3)==nTask,[6 7 10 12:15])=all_probes_mat2norm(all_probes_mat2norm(:,1)==myS(nS) & all_probes_mat2norm(:,3)==nTask,[6 7 10 12:15])./repmat(...
            nanmean(all_probes_mat2norm(all_probes_mat2norm(:,1)==myS(nS) & all_probes_mat2norm(:,3)==nTask & all_probes_mat2norm(:,4)==1,[6 7 10 12:15]),1),sum(all_probes_mat2norm(:,1)==myS(nS) & all_probes_mat2norm(:,3)==nTask),1);
    end
end


%% RT by probe
plotNames={'RT'};
plotCol=[12];
thismat=all_probes_mat2;
plotLims=[0.3 0.9];

for nPlot=1 %:3
    data=[];
    for i = 1:3
        for j = 1:2
            temp=squeeze(thismat(thismat(:,4)==i & thismat(:,3)==j,plotCol(nPlot)));
            %             tempS=squeeze(thismat(thismat(:,4)==i & thismat(:,3)==j,1));
            %             tempbyS=[]; myS=unique(tempS);
            %             for nS=1:length(myS)
            %                 tempbyS(nS)=nanmean(temp(tempS==myS(nS)));
            %             end
            %             data{i, j} = tempbyS;
            data{i, j} = temp;
        end
    end
    
    figure; set(gcf,'Position',[ 440   315   800   360]);
    for ntask=1:2
        subplot(1,2,ntask);
        h1 = raincloud_plot(data{1,ntask}, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
            'box_col_match', 0); h1{2}.SizeData=72; h1{2}.MarkerFaceAlpha=0.5;
        h2 = raincloud_plot(data{2,ntask}, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
        h2{2}.SizeData=72; h2{2}.MarkerFaceAlpha=0.5;
        h3 = raincloud_plot(data{3,ntask}, 'box_on', 1, 'color', Colors(3,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
        h3{2}.SizeData=72; h3{2}.MarkerFaceAlpha=0.5;
        if ntask==1
            set(gca,'XLim', plotLims(nPlot,:), 'YLim', [-3.5 6],'YTick',''); xlabel(plotNames{nPlot});
        else
            set(gca,'XLim', plotLims(nPlot,:), 'YLim', [-3.5 6],'YTick',''); xlabel(plotNames{nPlot});
        end
        box off
        format_fig;
    end
end
%export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Behav_RT_perProbeAndState.fig'])
%export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Behav_RT_perProbeAndState.eps'],'-r 300')

%% RT by probe
plotNames={'RT'};
plotCol=[12];
thismat=all_probes_mat2;
plotLims=[0.3 0.85];
MarkersT={'d','o'};
for nPlot=1 %:3
    data=[];
    for i = 1:3
        for j = 1:2
            temp=squeeze(thismat(thismat(:,4)==i & thismat(:,3)==j,plotCol(nPlot)));
            tempS=squeeze(thismat(thismat(:,4)==i & thismat(:,3)==j,1));
            tempbyS=[]; myS=unique(tempS);
            tempbyS_n=[];
            for nS=1:length(myS)
                tempbyS(nS)=nanmean(temp(tempS==myS(nS)));
                tempbyS_n(nS)=sum((tempS==myS(nS)));
            end
            data{i, j} = tempbyS;
            data_n{i, j} = tempbyS_n;
            
            myS=unique(thismat(:,1));
            datattest_n{i, j}=nan(1,length(myS));
            datattest{i, j}=nan(1,length(myS));
            for nS=1:length(myS)
                datattest{i, j}(nS)=nanmean(temp(tempS==myS(nS)));
                datattest_n{i, j}(nS)=sum((tempS==myS(nS)));
            end
            
            datattest{i, j}(datattest_n{i, j}<3)=NaN;
            %             data{i, j} = temp;
        end
        
        
        datattest{i, 3} = nanmean([datattest{i, 1} ; datattest{i, 2}],1);
    end
    
    figure; set(gcf,'Position',[ 440   315   800   360]);
    for ntask=1:2
        subplot(1,2,ntask);
        h1 = raincloud_plot(data{1,ntask}, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
            'box_col_match', 0,'size_data',data_n{1,ntask});
        h1{2}.MarkerFaceAlpha=0.5; h1{2}.Marker=MarkersT{ntask};
        %         h1{7}.MarkerFaceAlpha=0.5; h1{7}.Marker=MarkersT{ntask};  h1{7}.MarkerFaceColor=[1 1 1]*0.5;   h1{7}.MarkerEdgeColor=Colors(1,:);    h1{7}.SizeData=144;
        h2 = raincloud_plot(data{2,ntask}, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'size_data',data_n{2,ntask});
        h2{2}.MarkerFaceAlpha=0.5; h2{2}.Marker=MarkersT{ntask};
        h3 = raincloud_plot(data{3,ntask}, 'box_on', 1, 'color', Colors(3,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0,'size_data',data_n{3,ntask});
        h3{2}.MarkerFaceAlpha=0.5; h3{2}.Marker=MarkersT{ntask};
        
        if ntask==1
            set(gca,'XLim', plotLims(nPlot,:), 'YLim', [-4.5 7],'YTick',''); xlabel(plotNames{nPlot});
        else
            set(gca,'XLim', plotLims(nPlot,:), 'YLim', [-4.5 7],'YTick',''); xlabel(plotNames{nPlot});
        end
        box off
        format_fig;
    end
end
export_fig([path_fig 'Behav_RT_perSubjectAndState.fig'])
export_fig([path_fig 'Behav_RT_perSubjectAndState.eps'],'-r 300')


[h, pV, ~, stats]=ttest([datattest{2,3}],[datattest{1,3}]);
fprintf('post-hoc ttests (av task): MW vs ON t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)
[h, pV, ~, stats]=ttest([datattest{3,3}],[datattest{1,3}]);
fprintf('post-hoc ttests (av task): MB vs ON t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)
[h, pV, ~, stats]=ttest([datattest{3,3}],[datattest{2,3}]);
fprintf('post-hoc ttests (av task): MB vs MW t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)
%%
for nTask=1:2
    temp=all_probes_mat2(all_probes_mat2(:,3)==nTask,:);
    tbl=array2table(temp,'VariableNames',{'SubID','nBl','Task','MS','nPr','GO','NOGO','nGO','nNOGO','tdp','vig','rtGO','rtNOGO','dp','crit'});
    tbl.MS=categorical(tbl.MS);
    tbl.Task=categorical(tbl.Task);
    tbl.MS=reordercats(tbl.MS,[2 1 3]);
    %
    % mdl_0= fitlme(tbl,'crit~1+(1|SubID)');
    % mdl_1= fitlme(tbl,'crit~MS+(1|SubID)'); %wining model
    % mdl_2= fitlme(tbl,'crit~Task+MS+(1|SubID)');
    % mdl_3= fitlme(tbl,'crit~Task*MS+(1|SubID)');
    
    if nTask==1
        mdl2_0.F= fitlme(tbl,'rtGO~1+(1|SubID)');
        mdl2_1.F= fitlme(tbl,'rtGO~MS+(1|SubID)');
    else
        mdl2_0.D= fitlme(tbl,'rtGO~1+(1|SubID)');
        mdl2_1.D= fitlme(tbl,'rtGO~MS+(1|SubID)');
        
    end
    % mdl2_2{nTask}= fitlme(tbl,'rtGO~Task+MS+(1|SubID)'); %wining model
    % mdl2_3{nTask}= fitlme(tbl,'rtGO~Task*MS+(1|SubID)');
end


temp=all_probes_mat2; %(all_probes_mat2(:,5)==nTask,:);
tbl=array2table(temp,'VariableNames',{'SubID','nBl','Task','MS','nPr','GO','NOGO','nGO','nNOGO','tdp','vig','rtGO','rtNOGO','dp','crit'});
tbl.MS=categorical(tbl.MS);
tbl.Task=categorical(tbl.Task);
tbl.MS=reordercats(tbl.MS,[1 2 3]);

mdl_B_1_0= fitlme(tbl,'rtGO~1+(1|SubID)');
mdl_B_1_1= fitlme(tbl,'rtGO~Task+(1|SubID)');
mdl_B_1_2= fitlme(tbl,'rtGO~Task+MS+(1|SubID)'); % winning model
mdl_B_1_3= fitlme(tbl,'rtGO~Task*MS+(1|SubID)');

mdl_B_2_0= fitlme(tbl,'GO~1+(1|SubID)');
mdl_B_2_1= fitlme(tbl,'GO~Task+(1|SubID)');
mdl_B_2_2= fitlme(tbl,'GO~Task+MS+(1|SubID)');  % winning model
mdl_B_2_3= fitlme(tbl,'GO~Task*MS+(1|SubID)');


mdl_B_3_0= fitlme(tbl,'NOGO~1+(1|SubID)');
mdl_B_3_1= fitlme(tbl,'NOGO~Task+(1|SubID)');
mdl_B_3_2= fitlme(tbl,'NOGO~Task+MS+(1|SubID)');    % winning model
mdl_B_3_3= fitlme(tbl,'NOGO~Task*MS+(1|SubID)');  

tbl.MS=reordercats(tbl.MS,[2 1 3]);
mdl_B_1_2b= fitlme(tbl,'rtGO~Task+MS+(1|SubID)'); % winning model
mdl_B_2_2b= fitlme(tbl,'GO~Task+MS+(1|SubID)');  % winning model
mdl_B_3_2b= fitlme(tbl,'NOGO~Task+MS+(1|SubID)');   % winning model

% mdl_B_3_3= fitlme(tbl,'NOGO~Task+MS+(1|SubID)');

%%
tbl.MS=reordercats(tbl.MS,[2 1 3]);
tbl_norm=tbl;
tbl_norm.SubID=categorical(tbl_norm.SubID);
tbl_norm.Task=categorical(tbl_norm.Task);
tbl_norm.MS=categorical(tbl_norm.MS);
SubID=unique(tbl_norm.SubID);

tbl_norm.GO_ON=nan(size(tbl_norm,1),1);
tbl_norm.NOGO_ON=nan(size(tbl_norm,1),1);
tbl_norm.rtGO_ON=nan(size(tbl_norm,1),1);
for nS=1:length(SubID)
    for nT=1:2
        tbl_norm.rtGO(tbl_norm.SubID==SubID(nS)  & tbl_norm.Task==num2str(nT))=tbl_norm.rtGO(tbl_norm.SubID==SubID(nS)  & tbl_norm.Task==num2str(nT))./mean(tbl_norm.rtGO(tbl_norm.SubID==SubID(nS) & tbl_norm.MS=='1'  & tbl_norm.Task==num2str(nT)));
        tbl_norm.GO(tbl_norm.SubID==SubID(nS)  & tbl_norm.Task==num2str(nT))=tbl_norm.GO(tbl_norm.SubID==SubID(nS)  & tbl_norm.Task==num2str(nT))./mean(tbl_norm.GO(tbl_norm.SubID==SubID(nS) & tbl_norm.MS=='1'  & tbl_norm.Task==num2str(nT)));
        tbl_norm.NOGO(tbl_norm.SubID==SubID(nS)  & tbl_norm.Task==num2str(nT))=tbl_norm.NOGO(tbl_norm.SubID==SubID(nS)  & tbl_norm.Task==num2str(nT))./mean(tbl_norm.NOGO(tbl_norm.SubID==SubID(nS) & tbl_norm.MS=='1'  & tbl_norm.Task==num2str(nT)));
    end
end


mdl2_B_1_0= fitlme(tbl_norm,'rtGO~1+(1|SubID)');
mdl2_B_1_1= fitlme(tbl_norm,'rtGO~Task+(1|SubID)');
mdl2_B_1_2= fitlme(tbl_norm,'rtGO~Task+MS+(1|SubID)'); 
mdl2_B_1_3= fitlme(tbl_norm,'rtGO~Task*MS+(1|SubID)'); % winning model

mdl2_B_2_0= fitlme(tbl_norm,'GO~1+(1|SubID)');
mdl2_B_2_1= fitlme(tbl_norm,'GO~Task+(1|SubID)');
mdl2_B_2_2= fitlme(tbl_norm,'GO~Task+MS+(1|SubID)');  
mdl2_B_2_3= fitlme(tbl_norm,'GO~Task*MS+(1|SubID)');  % winning model


mdl2_B_3_0= fitlme(tbl_norm,'NOGO~1+(1|SubID)');
mdl2_B_3_1= fitlme(tbl_norm,'NOGO~Task+(1|SubID)');
mdl2_B_3_2= fitlme(tbl_norm,'NOGO~Task+MS+(1|SubID)');   
mdl2_B_3_3= fitlme(tbl_norm,'NOGO~Task*MS+(1|SubID)');   % winning model

tbl_norm.MS=reordercats(tbl_norm.MS,[2 1 3]);
mdl2_B_1_2b= fitlme(tbl_norm,'rtGO~Task+MS+(1|SubID)'); % winning model
mdl2_B_2_2b= fitlme(tbl_norm,'GO~Task+MS+(1|SubID)');  % winning model
mdl2_B_3_2b= fitlme(tbl_norm,'NOGO~Task+MS+(1|SubID)');   % winning model
