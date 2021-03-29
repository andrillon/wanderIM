%%
clear all
close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))
addpath(genpath(path_export))
addpath(genpath(path_RainCloudPlot))
% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

% state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
% cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];

% load([pwd filesep 'paper_SubID'])
%%
all_probes_mat=[];
all_probes_mat_fullgo=[];
all_probes_mat_fullnogo=[];
all_probes_mat_fullrtgo=[];
all_probes_mat_fullrtnogo=[];
all_probes_mat2=[];
all_probes_mat3=[];


all_probes_mat_far=[];
all_probes_mat_fullgo_far=[];
all_probes_mat_fullnogo_far=[];
all_probes_mat_fullrtgo_far=[];
all_probes_mat_fullrtnogo_far=[];
all_probes_mat2_far=[];
all_probes_mat3_far=[];

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
            number_STD(nSc,nbl,npr)=length(thesegotrials(end-8:end));
            number_DEV(nSc,nbl,npr)=length(thesenogotrials(end:end));
            
            temp_testres=these_trials([thesegotrials(end-8:end) ; thesenogotrials(end:end)],:);
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
            
            temp_testres=these_trials([thesegotrials(end-8:end) ; thesenogotrials(end:end)],:);
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
            
            
            %%% 20-10 trials
            thesegotrials=find(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx & these_trials(:,5)~=3);
            thesenogotrials=find(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx & these_trials(:,5)==3);
            number_STD_far(nSc,nbl,npr)=length(thesegotrials(end-17:end-9));
            number_DEV_far(nSc,nbl,npr)=length(thesenogotrials(end-1:end-1));
            
            temp_testres=these_trials([thesegotrials(end-17:end-9) ; thesenogotrials(end-1:end-1)],:);
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
            all_probes_mat_far=[all_probes_mat_far ; [nSc nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_go tcorr_nogo num_go num_nogo tdp (these_probes(npr,37)) rt_go rt_nogo]];
            
            temp_testres=these_trials([thesegotrials(end-17:end-9) ; thesenogotrials(end-1:end-1)],:);
            tcorr_go=(temp_testres(~isnan(temp_testres(:,12)),12))';%/corr_go(nSc,these_probes(npr,5));
            rtcorr_go=(temp_testres(~isnan(temp_testres(:,12)),10)-temp_testres(~isnan(temp_testres(:,12)),8))';%/corr_go(nSc,these_probes(npr,5));
            tcorr_nogo=(temp_testres(~isnan(temp_testres(:,11)),11))';%/corr_go(nSc,these_probes(npr,5));
            rtcorr_nogo=(temp_testres(~isnan(temp_testres(:,11)),10)-temp_testres(~isnan(temp_testres(:,11)),8))';%/corr_go(nSc,these_probes(npr,5));
            all_probes_mat_fullgo_far=[all_probes_mat_fullgo_far ; [nSc nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr tcorr_go]];
            all_probes_mat_fullnogo_far=[all_probes_mat_fullnogo_far ; [nSc nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr tcorr_nogo]];
            all_probes_mat_fullrtgo_far=[all_probes_mat_fullrtgo_far ; [nSc nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr rtcorr_go]];
            all_probes_mat_fullrtnogo_far=[all_probes_mat_fullrtnogo_far ; [nSc nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr rtcorr_nogo]];
            
            [dprime, crit]=calc_dprime2(tcorr_go==1,tcorr_nogo==0);
            all_probes_mat2_far=[all_probes_mat2_far ; [nSc nbl these_probes(npr,5) these_probes(npr,32) npr nanmean(tcorr_go) nanmean(tcorr_nogo) num_go num_nogo tdp (these_probes(npr,38)) rt_go rt_nogo dprime crit]];
            all_probes_mat3_far=[all_probes_mat3_far ; [nSc nbl npr these_probes(npr,5) these_probes(npr,32) npr nanmean(tcorr_go) nanmean(tcorr_nogo) num_go num_nogo tdp (these_probes(npr,38)) rt_go rt_nogo dprime crit these_probes(npr,31:38)]];
            
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
% export_fig([path_fig 'Behav_RT_perSubjectAndState.fig'])
% export_fig([path_fig 'Behav_RT_perSubjectAndState.eps'],'-r 300')


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


%%
temp=all_probes_mat2; %(all_probes_mat2(:,5)==nTask,:);
tbl=array2table(temp,'VariableNames',{'SubID','nBl','Task','MS','nPr','GO','NOGO','nGO','nNOGO','tdp','vig','rtGO','rtNOGO','dp','crit'});
tbl.MS=categorical(tbl.MS);
tbl.Task=categorical(tbl.Task);
tbl.MS=reordercats(tbl.MS,[1 2 3]);



temp=all_probes_mat2_far; %(all_probes_mat2(:,5)==nTask,:);
tbl_far=array2table(temp,'VariableNames',{'SubID','nBl','Task','MS','nPr','GO','NOGO','nGO','nNOGO','tdp','vig','rtGO','rtNOGO','dp','crit'});
tbl_far.MS=categorical(tbl_far.MS);
tbl_far.Task=categorical(tbl_far.Task);
tbl_far.MS=reordercats(tbl_far.MS,[1 2 3]);

figure;
hp=[];
set(gcf,'Position',[374   125   923   860]);
for ntask=1:2
    subplot(3,2,ntask); format_fig;
    if ntask==1
        title('Faces')
    else
        title('Digits')
    end
    for nWin=1:2
        for nstate=1:3
            if nWin==1
                this_table=tbl_far;
                this_col=[1 1 1; Colors(nstate,:)];
            elseif nWin==2
                this_table=tbl;
                this_col=[Colors(nstate,:)];
            end
            temp1=this_table.GO(this_table.Task==num2str(ntask) & this_table.MS==num2str(nstate));
            temp2=this_table.SubID(this_table.Task==num2str(ntask) & this_table.MS==num2str(nstate));
            mean1=grpstats(temp1,temp2);
            [hp{nWin}(nstate)]=simpleBarPlot(nWin+(nstate-2)*0.25,100*(1-mean1),this_col,0.21,'k');
        end
    end
    ylim([0 6.5])
    xlim([0.2 2.8])
    set(gca,'XTick',1:2,'XTickLabel',{'[-20, -10]','[-10, 0]'});
%         xlabel('Time Window')
    ylabel('Misses (%)')
    if ntask==2
        legend(hp{2},{'ON','MW','MB'},'Position',[0.9 0.5 0.1 0.1],'Box','off');
    end
    subplot(3,2,2+ntask); format_fig;
    for nWin=1:2
        for nstate=1:3
            if nWin==1
                this_table=tbl_far;
                this_col=[1 1 1; Colors(nstate,:)];
            elseif nWin==2
                this_table=tbl;
                this_col=[Colors(nstate,:)];
            end
            temp1=this_table.NOGO(this_table.Task==num2str(ntask) & this_table.MS==num2str(nstate));
            temp2=this_table.SubID(this_table.Task==num2str(ntask) & this_table.MS==num2str(nstate));
            mean1=grpstats(temp1,temp2);
            simpleBarPlot(nWin+(nstate-2)*0.25,100*(1-mean1),this_col,0.21,'k');
        end
    end
    ylim([0 60])
    xlim([0.2 2.8])
        set(gca,'XTick',1:2,'XTickLabel',{'[-20, -10]','[-10, 0]'});
%     xlabel('Time Window')
    ylabel('FAs (%)')
    
    subplot(3,2,4+ntask); format_fig;
    for nWin=1:2
        for nstate=1:3
            if nWin==1
                this_table=tbl_far;
                this_col=[1 1 1; Colors(nstate,:)];
            elseif nWin==2
                this_table=tbl;
                this_col=[Colors(nstate,:)];
            end
            temp1=this_table.rtGO(this_table.Task==num2str(ntask) & this_table.MS==num2str(nstate));
            temp2=this_table.SubID(this_table.Task==num2str(ntask) & this_table.MS==num2str(nstate));
            mean1=grpstats(temp1,temp2);
            simpleBarPlot(nWin+(nstate-2)*0.25,mean1,this_col,0.21,'k');
        end
    end
    if ntask==1
    ylim([0.55 0.7])
    else
    ylim([0.45 0.6])
    end
    xlim([0.2 2.8])
    set(gca,'XTick',1:2,'XTickLabel',{'[-20, -10]','[-10, 0]'});
    xlabel('Time Window')
    ylabel('RT GO (s)')
end