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

state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];
% load([pwd filesep 'paper_SubID'])

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
nc=0;
all_MS=[];
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    all_MS=[all_MS ; [n*ones(size(probe_res,1),1) probe_res(:,32)]];
    probe_res(probe_res(:,32)==4,32)=3;
    SubID=SubjectInfo.subID;
    if ~ismember(SubID,GoodSubID)
        continue;
    end
    nc=nc+1;
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
        [dprime_test(nc,nbt), crit_test(nc,nbt)]=calc_dprime2((tp_gos==1),(tp_nogos==0));
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
        for nstate=1:3
            pption(nstate)=mean(temp(:,3)==nstate);
            pption2(nstate)=mean(temp(temp(:,3)==2,4)==nstate);
        end
        pption_MS(nc,nbt,:)=pption;
        pption_Ori(nc,nbt,:)=pption2;
    end
    for nb=1:max(probe_res(:,3))
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
            if npr==1
                last_pr_tridx=0;
            else
                last_pr_tridx=these_probes(npr-1,6);
            end
            %                         last_pr_tridx=this_pr_tridx-20;
            thesegotrials=find(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx & these_trials(:,5)~=3);
            thesenogotrials=find(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx & these_trials(:,5)==3);
            number_STD(nc,nbl,npr)=length(thesegotrials(end-17:end));
            number_DEV(nc,nbl,npr)=length(thesenogotrials(end-1:end));
            
            temp_testres=these_trials([thesegotrials(end-17:end) ; thesenogotrials(end-1:end)],:);
            %             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            tcorr_go=nanmean(temp_testres(:,12));%/corr_go(nc,these_probes(npr,5));
            tcorr_nogo=nanmean(temp_testres(:,11));%/corr_nogo(nc,these_probes(npr,5));
            num_go=sum(~isnan(temp_testres(:,12)));
            num_nogo=sum(~isnan(temp_testres(:,11)));
            rt_go=nanmean(temp_testres(~isnan(temp_testres(:,12)),10)-temp_testres(~isnan(temp_testres(:,12)),8));
            rt_nogo=nanmean(temp_testres(~isnan(temp_testres(:,11)),10)-temp_testres(~isnan(temp_testres(:,11)),8));
            
            tdp=calc_dprime2((temp_testres(~isnan(temp_testres(:,12)),12)==1),(temp_testres(~isnan(temp_testres(:,11)),11)==0));
            all_probes_mat=[all_probes_mat ; [nc nbl these_probes(npr,5) these_probes(npr,32) npr tcorr_go tcorr_nogo num_go num_nogo tdp (these_probes(npr,37)) rt_go rt_nogo]];
            
            temp_testres=these_trials([thesegotrials(end-17:end) ; thesenogotrials(end-1:end)],:);
            tcorr_go=(temp_testres(~isnan(temp_testres(:,12)),12))';%/corr_go(nc,these_probes(npr,5));
            rtcorr_go=(temp_testres(~isnan(temp_testres(:,12)),10)-temp_testres(~isnan(temp_testres(:,12)),8))';%/corr_go(nc,these_probes(npr,5));
            tcorr_nogo=(temp_testres(~isnan(temp_testres(:,11)),11))';%/corr_go(nc,these_probes(npr,5));
            rtcorr_nogo=(temp_testres(~isnan(temp_testres(:,11)),10)-temp_testres(~isnan(temp_testres(:,11)),8))';%/corr_go(nc,these_probes(npr,5));
            all_probes_mat_fullgo=[all_probes_mat_fullgo ; [nc nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr tcorr_go]];
            all_probes_mat_fullnogo=[all_probes_mat_fullnogo ; [nc nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr tcorr_nogo]];
            all_probes_mat_fullrtgo=[all_probes_mat_fullrtgo ; [nc nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr rtcorr_go]];
            all_probes_mat_fullrtnogo=[all_probes_mat_fullrtnogo ; [nc nbl these_probes(npr,5) these_probes(npr,32)  these_probes(npr,33) npr rtcorr_nogo]];
            
            [dprime, crit]=calc_dprime2(tcorr_go==1,tcorr_nogo==0);
            all_probes_mat2=[all_probes_mat2 ; [nc nbl these_probes(npr,5) these_probes(npr,32) npr nanmean(tcorr_go) nanmean(tcorr_nogo) num_go num_nogo tdp (these_probes(npr,38)) rt_go rt_nogo dprime crit]];
            all_probes_mat3=[all_probes_mat3 ; [nc nbl npr these_probes(npr,5) these_probes(npr,32) npr nanmean(tcorr_go) nanmean(tcorr_nogo) num_go num_nogo tdp (these_probes(npr,38)) rt_go rt_nogo dprime crit these_probes(npr,31:38)]];
        end
    end
end

%%
allpption_alongtask=[];
allpption_alongtask2=[];
forANOVA=[];
forANOVAgr=[];
forANOVA2=[];
for n=1:nc
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
            if length(this_sub_nogo)>=3
                [dprime_bysubj_state(n,ntask,nSta), crit_bysubj_state(n,ntask,nSta)] = calc_dprime2(this_sub_go,this_sub_nogo==0);
                go_bysubj_state(n,ntask,nSta)=nanmean(this_sub_go);
                nogo_bysubj_state(n,ntask,nSta)=nanmean(this_sub_nogo);
                corr_bysubj_state(n,ntask,nSta)=nanmean([this_sub_go this_sub_nogo]);
                
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
                corr_bysubj_state(n,ntask,nSta)=NaN;
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
        fatigue_alongtask(n,nblock)=mean(this_sub_part2);
        for nstate=1:3
            %                 if nstate~=2
            pption_alongtask(n,nblock,nstate)=mean(this_sub_part==nstate);
            %                 else
            %                 pption_alongtask(n,nblock,nstate)=mean(this_sub_part==nstate & this_sub_part3==nstate);
            %                 end
            allpption_alongtask=[allpption_alongtask ; [n nblock nstate mean(this_sub_part==nstate)]];
            allpption_alongtask2=[allpption_alongtask2 ; [n nblock nstate mean(this_sub_part2)]];
        end
    end
    for npart=1:10 %length(partidx)
        this_sub_part=all_probes_mat3(all_probes_mat3(:,1)==n & ismember(all_probes_mat3(:,3),npart),end-6);
        this_sub_part2=all_probes_mat3(all_probes_mat3(:,1)==n & ismember(all_probes_mat3(:,3),npart),end);
        fatigue_alongblock(n,npart)=mean(this_sub_part2);
        for nstate=1:3
            pption_alongblock(n,npart,nstate)=mean(this_sub_part==nstate);
        end
    end
end

%% Figure 1c
% read into cell array of the appropriate dimensions
for i = 1:3
    for j = 1:2
        data{i, j} = pption_MS(:,j,i)';
    end
end

% make figure
figure;
set(gcf,'position',[11   513   375   445])
h   = rm_raincloud(data, Colors(1:2,:),0, 'ks', [],100);
set(gca, 'XLim', [0 1],'YTick','','Ylim',[-3 13]);
% title(['Figure M9' newline 'Repeated measures raincloud plot']);

% change one subset to new colour and alter dot size
for i=1:3
    % scatter
    h.s{i, 2}.MarkerFaceColor   = Colors(i,:);
    h.s{i, 2}.MarkerEdgeColor   = [1 1 1]*0.5;
    h.s{i, 2}.SizeData          = 100;
    h.s{i, 2}.LineWidth          = 2;
    h.s{i,2}.YData=h.s{i,2}.YData-0.2;
    
    % patch
    h.p{i, 2}.FaceColor     = Colors(i,:);
    h.p{i, 2}.LineWidth          = 2;
    
    % scatter mean
    h.m(i,2).MarkerFaceAlpha=0.7;
            h.m(i,2).MarkerFaceColor=[1 1 1]*0.5;
            h.m(i,2).MarkerEdgeColor=[0 0 0];%[190 90 255]./255;
            h.m(i,2).SizeData=300;
            h.m(i,2).LineWidth=4;
    % line mean
    if i<3
        h.l(i,2).Color=[1 1 1]*0;
    end
end

for i=1:3
    h.s{i, 1}.Marker         = 'd';
    h.s{i, 1}.MarkerEdgeColor   = [1 1 1]*0.5;
    h.s{i, 1}.MarkerFaceColor   = Colors(i,:);
    h.s{i, 1}.SizeData          = 100;
    h.s{i, 1}.LineWidth          = 2;
    h.s{i,1}.YData=h.s{i,1}.YData+0.2;
    
    h.m(i, 1).MarkerEdgeColor   = 'none';
    h.m(i, 1).MarkerFaceColor   = Colors(i,:);
    
    h.p{i, 1}.FaceAlpha   = 0;
    h.p{i, 1}.EdgeColor= Colors(i,:);
    h.p{i, 1}.LineWidth          = 1;
    
    
    h.m(i,1).Marker='d';
    h.m(i,1).MarkerFaceAlpha=0.7;
    h.m(i,1).MarkerFaceColor=[1 1 1]*0.5;
    h.m(i,1).MarkerEdgeColor=[0 0 0];
    h.m(i,1).SizeData=300;
    h.m(i,1).LineWidth=4;
   
    if i<3
        h.l(i,1).Color=[1 1 1]*0;
        h.l(i,1).LineStyle='-';
    end
end
format_fig;

%export_fig([path_fig filesep 'WanderIM_behav_pption_byStateAndTask.eps'],'-r 300')


%%
temp=all_probes_mat2(all_probes_mat2(:,4)<4,:);
tbl=array2table(temp,'VariableNames',{'SubID','nBl','Task','MS','nPr','GO','NOGO','nGO','nNOGO','tdp','vig','rtGO','rtNOGO','dp','crit'});
tbl.MS=categorical(tbl.MS);
tbl.Task=categorical(tbl.Task);
tbl.MS=reordercats(tbl.MS,[1 2 3]);


mdl6_0= fitlme(tbl,'vig~1+(1|SubID)');
mdl6_1= fitlme(tbl,'vig~Task+(1|SubID)');
mdl6_2= fitlme(tbl,'vig~Task+MS+(1|SubID)'); %wining model
mdl6_3= fitlme(tbl,'vig~Task*MS+(1|SubID)'); %wining model

tbl.MS=reordercats(tbl.MS,[2 1 3]);
mdl6_2b= fitlme(tbl,'vig~Task+MS+(1|SubID)'); %wining model

%%
data=[];
for i = 1:3
    for j = 1:2
        temp=all_probes_mat2(all_probes_mat2(:,3)==j & all_probes_mat2(:,4)==i,11);
        tempS=squeeze(all_probes_mat2(all_probes_mat2(:,3)==j & all_probes_mat2(:,4)==i,1));
        tempbyS=[]; myS=unique(tempS);
        for nS=1:length(myS)
            tempbyS(nS)=nanmean(temp(tempS==myS(nS)));
        end
        data{i, j} = tempbyS;
        
        myS=unique(all_probes_mat2(:,1));
        datattest_n{i, j}=nan(1,length(myS));
        datattest{i, j}=nan(1,length(myS));
        for nS=1:length(myS)
            datattest{i, j}(nS)=nanmean(temp(tempS==myS(nS)));
            datattest_n{i, j}(nS)=sum((tempS==myS(nS)));
        end
        datattest{i, j}(datattest_n{i, j}<3)=NaN;
    end
            datattest{i, 3} = nanmean([datattest{i, 1} ; datattest{i, 2}],1);

end

% % % % make figure
% % % figure;
% % % set(gcf,'position',[11   513   375   445])
% % % h   = rm_raincloud(data, Colors(1:2,:));
% % % set(gca, 'XLim', [1 4],'YTick','');
% % % 
% % % % change one subset to new colour and alter dot size
% % % for i=1:3
% % %     % scatter
% % %     h.s{i, 2}.MarkerFaceColor   = Colors(i,:);
% % %     h.s{i, 2}.MarkerEdgeColor   = [1 1 1]*0.5;
% % %     h.s{i, 2}.SizeData          = 100;
% % %     h.s{i, 2}.LineWidth          = 2;
% % %     h.s{i,2}.YData=h.s{i,2}.YData-0.1;
% % %     
% % %     % patch
% % %     h.p{i, 2}.FaceColor     = Colors(i,:);
% % %     h.p{i, 2}.LineWidth          = 2;
% % %     
% % %     % scatter mean
% % %   h.m(i,2).MarkerFaceAlpha=0.7;
% % %             h.m(i,2).MarkerFaceColor=[1 1 1]*0.5;
% % %             h.m(i,2).MarkerEdgeColor=[0 0 0];%[190 90 255]./255;
% % %             h.m(i,2).SizeData=300;
% % %             h.m(i,2).LineWidth=4;
% % %     % line mean
% % %     if i<3
% % %         h.l(i,2).Color=[1 1 1]*0;
% % %     end
% % % end
% % % 
% % % for i=1:3
% % %     h.s{i, 1}.Marker         = 'd';
% % %     h.s{i, 1}.MarkerEdgeColor   = [1 1 1]*0.5;
% % %     h.s{i, 1}.MarkerFaceColor   = Colors(i,:);
% % %     h.s{i, 1}.SizeData          = 100;
% % %     h.s{i, 1}.LineWidth          = 2;
% % %     h.s{i,1}.YData=h.s{i,1}.YData+0.1;
% % %     
% % %     h.m(i, 1).MarkerEdgeColor   = 'none';
% % %     h.m(i, 1).MarkerFaceColor   = Colors(i,:);
% % %     
% % %     h.p{i, 1}.FaceAlpha   = 0;
% % %     h.p{i, 1}.EdgeColor= Colors(i,:);
% % %     h.p{i, 1}.LineWidth          = 1;
% % %     
% % %     
% % %     h.m(i,1).Marker='d';
% % %     h.m(i,1).MarkerFaceAlpha=0.7;
% % %     h.m(i,1).MarkerFaceColor=[1 1 1]*0.5;
% % %     h.m(i,1).MarkerEdgeColor=[0 0 0];
% % %     h.m(i,1).SizeData=300;
% % %     h.m(i,1).LineWidth=4;
% % %     if i<3
% % %         h.l(i,1).Color=[1 1 1]*0;
% % %         h.l(i,1).LineStyle='-';
% % %     end
% % % end
% % % format_fig;
% % % %export_fig([path_fig filesep 'WanderIM_behav_fatigue_byStateAndTask.eps'],'-r 300')

[h, pV, ~, stats]=ttest([datattest{2,3}],[datattest{1,3}]);
fprintf('post-hoc ttests (av task): MW vs ON t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)
[h, pV, ~, stats]=ttest([datattest{3,3}],[datattest{1,3}]);
fprintf('post-hoc ttests (av task): MB vs ON t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)
[h, pV, ~, stats]=ttest([datattest{3,3}],[datattest{2,3}]);
fprintf('post-hoc ttests (av task): MB vs MW t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)
%% FA (Figure 2b)
% read into cell array of the appropriate dimensions
data=[];
data_n=[];
for i = 1:3
    for j = 1:2
%         data{i, j} = 1-nogo_bysubj_state(:,j,i);
%         data{i, j} = data{i, j} (~isnan(data{i, j}));
        
        temp=all_probes_mat2(all_probes_mat2(:,3)==j & all_probes_mat2(:,4)==i,7);
        tempS=squeeze(all_probes_mat2(all_probes_mat2(:,3)==j & all_probes_mat2(:,4)==i,1));
        tempbyS=[]; myS=unique(tempS);
        tempbyS_n=[]; 
        for nS=1:length(myS)
            tempbyS(nS)=nanmean(temp(tempS==myS(nS)));
            tempbyS_n(nS)=sum(~isnan(temp(tempS==myS(nS))));
        end
        data{i, j} = 1-tempbyS;
        data_n{i, j} = 50+100*minmax(tempbyS_n);
    end
end

% make figure
figure;
set(gcf,'position',[11   513   375   445])
h   = rm_raincloud(data, Colors(1:2,:),0, 'ks', [],data_n);
set(gca, 'XLim', [0 1],'YTick','','Ylim',[-3 13]);
% title(['Figure M9' newline 'Repeated measures raincloud plot']);

% change one subset to new colour and alter dot size
for i=1:3
    % scatter
    h.s{i, 2}.MarkerFaceColor   = Colors(i,:);
    h.s{i, 2}.MarkerEdgeColor   = [1 1 1]*0.5;
%     h.s{i, 2}.SizeData          = 100;
    h.s{i, 2}.LineWidth          = 2;
    h.s{i,2}.YData=h.s{i,2}.YData-0.2;
    
    % patch
    h.p{i, 2}.FaceColor     = Colors(i,:);
    h.p{i, 2}.LineWidth          = 2;
    
    % scatter mean
    h.m(i,2).MarkerFaceAlpha=0.7;
            h.m(i,2).MarkerFaceColor=[1 1 1]*0.5;
            h.m(i,2).MarkerEdgeColor=[0 0 0];%[190 90 255]./255;
            h.m(i,2).SizeData=300;
            h.m(i,2).LineWidth=4;
    % line mean
    if i<3
        h.l(i,2).Color=[1 1 1]*0;
    end
end

for i=1:3
    h.s{i, 1}.Marker         = 'd';
    h.s{i, 1}.MarkerEdgeColor   = [1 1 1]*0.5;
    h.s{i, 1}.MarkerFaceColor   = Colors(i,:);
%     h.s{i, 1}.SizeData          = 100;
    h.s{i, 1}.LineWidth          = 2;
    h.s{i,1}.YData=h.s{i,1}.YData+0.2;
    
    h.m(i, 1).MarkerEdgeColor   = 'none';
    h.m(i, 1).MarkerFaceColor   = Colors(i,:);
    
    h.p{i, 1}.FaceAlpha   = 0;
    h.p{i, 1}.EdgeColor= Colors(i,:);
    h.p{i, 1}.LineWidth          = 1;
    
    
    h.m(i,1).Marker='d';
    h.m(i,1).MarkerFaceAlpha=0.7;
    h.m(i,1).MarkerFaceColor=[1 1 1]*0.5;
    h.m(i,1).MarkerEdgeColor=[0 0 0];
    h.m(i,1).SizeData=300;
    h.m(i,1).LineWidth=4;
   
    if i<3
        h.l(i,1).Color=[1 1 1]*0;
        h.l(i,1).LineStyle='-';
    end
end
format_fig;

%export_fig([path_fig filesep 'WanderIM_behav_NOGO_byStateAndTask.eps'],'-r 300')


%% MISSES (Figure 2s)
% read into cell array of the appropriate dimensions
data=[];
data_n=[];
for i = 1:3
    for j = 1:2
%         data{i, j} = 1-nogo_bysubj_state(:,j,i);
%         data{i, j} = data{i, j} (~isnan(data{i, j}));
        
        temp=all_probes_mat2(all_probes_mat2(:,3)==j & all_probes_mat2(:,4)==i,6);
        tempS=squeeze(all_probes_mat2(all_probes_mat2(:,3)==j & all_probes_mat2(:,4)==i,1));
        tempbyS=[]; myS=unique(tempS);
        tempbyS_n=[]; 
        for nS=1:length(myS)
            tempbyS(nS)=nanmean(temp(tempS==myS(nS)));
            tempbyS_n(nS)=sum(~isnan(temp(tempS==myS(nS))));
        end
        data{i, j} = 1-tempbyS;
        data_n{i, j} = 50+100*minmax(tempbyS_n);
    end
end

% make figure
figure;
set(gcf,'position',[11   513   375   445])
h   = rm_raincloud(data, Colors(1:2,:),0, 'ks', [],data_n);
set(gca, 'XLim', [0 0.25],'YTick','','Ylim',[-50 200]);
% title(['Figure M9' newline 'Repeated measures raincloud plot']);

% change one subset to new colour and alter dot size
for i=1:3
    % scatter
    h.s{i, 2}.MarkerFaceColor   = Colors(i,:);
    h.s{i, 2}.MarkerEdgeColor   = [1 1 1]*0.5;
%     h.s{i, 2}.SizeData          = 100;
    h.s{i, 2}.LineWidth          = 2;
    h.s{i,2}.YData=h.s{i,2}.YData-0.2;
    
    % patch
    h.p{i, 2}.FaceColor     = Colors(i,:);
    h.p{i, 2}.LineWidth          = 2;
    
    % scatter mean
    h.m(i,2).MarkerFaceAlpha=0.7;
            h.m(i,2).MarkerFaceColor=[1 1 1]*0.5;
            h.m(i,2).MarkerEdgeColor=[0 0 0];%[190 90 255]./255;
            h.m(i,2).SizeData=300;
            h.m(i,2).LineWidth=4;
    % line mean
    if i<3
        h.l(i,2).Color=[1 1 1]*0;
    end
end

for i=1:3
    h.s{i, 1}.Marker         = 'd';
    h.s{i, 1}.MarkerEdgeColor   = [1 1 1]*0.5;
    h.s{i, 1}.MarkerFaceColor   = Colors(i,:);
%     h.s{i, 1}.SizeData          = 100;
    h.s{i, 1}.LineWidth          = 2;
    h.s{i,1}.YData=h.s{i,1}.YData+0.2;
    
    h.m(i, 1).MarkerEdgeColor   = 'none';
    h.m(i, 1).MarkerFaceColor   = Colors(i,:);
    
    h.p{i, 1}.FaceAlpha   = 0;
    h.p{i, 1}.EdgeColor= Colors(i,:);
    h.p{i, 1}.LineWidth          = 1;
    
    
    h.m(i,1).Marker='d';
    h.m(i,1).MarkerFaceAlpha=0.7;
    h.m(i,1).MarkerFaceColor=[1 1 1]*0.5;
    h.m(i,1).MarkerEdgeColor=[0 0 0];
    h.m(i,1).SizeData=300;
    h.m(i,1).LineWidth=4;
   
    if i<3
        h.l(i,1).Color=[1 1 1]*0;
        h.l(i,1).LineStyle='-';
    end
end
format_fig;


export_fig([path_fig filesep 'WanderIM_behav_GO_byStateAndTask.eps'],'-r 300')

%%% post hoc tests
% [h, pV, ~, stats]=ttest([datattest{2,1} ; datattest{2,2}],[datattest{1,1} ; datattest{1,2}]);
% fprintf('post-hoc ttests (av task): MW vs ON t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)
% [h, pV, ~, stats]=ttest([datattest{3,1} ; datattest{3,2}],[datattest{1,1} ; datattest{1,2}]);
% fprintf('post-hoc ttests (av task): MB vs ON t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)
% [h, pV, ~, stats]=ttest([datattest{3,1} ; datattest{3,2}],[datattest{2,1} ; datattest{2,2}]);
% fprintf('post-hoc ttests (av task): MB vs MW t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)

[h, pV, ~, stats]=ttest([datattest{2,3}],[datattest{1,3}]);
fprintf('post-hoc ttests (av task): MW vs ON t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)
[h, pV, ~, stats]=ttest([datattest{3,3}],[datattest{1,3}]);
fprintf('post-hoc ttests (av task): MB vs ON t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)
[h, pV, ~, stats]=ttest([datattest{3,3}],[datattest{2,3}]);
fprintf('post-hoc ttests (av task): MB vs MW t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)

