%%
clear all
close all

run ../localdef_wanderIM.m
addpath(genpath(lscpTools_path))
addpath(genpath(path_RainCloudPlot))
addpath(genpath(path_export));
%% load data
filename = '/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_localsleep_pup_Dec21_v5.txt';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

res_table = table;
res_table.SubID = dataArray{:, 1};
res_table.BlockN = dataArray{:, 2};
res_table.TrialN = dataArray{:, 3};
res_table.ProbeN = dataArray{:, 4};
res_table.DistProbe = dataArray{:, 5};
res_table.Task = dataArray{:, 6};
res_table.StimCat = dataArray{:, 7};
res_table.Perf = dataArray{:, 8};
res_table.RCode = dataArray{:, 9};
res_table.RT = dataArray{:, 10};
res_table.W_Fp1 = dataArray{:, 11};
res_table.W_Fz = dataArray{:, 12};
res_table.W_F3 = dataArray{:, 13};
res_table.W_F7 = dataArray{:, 14};
res_table.W_FT9 = dataArray{:, 15};
res_table.W_FC5 = dataArray{:, 16};
res_table.W_FC1 = dataArray{:, 17};
res_table.W_C3 = dataArray{:, 18};
res_table.W_T7 = dataArray{:, 19};
res_table.W_TP9 = dataArray{:, 20};
res_table.W_CP5 = dataArray{:, 21};
res_table.W_CP1 = dataArray{:, 22};
res_table.W_Pz = dataArray{:, 23};
res_table.W_P3 = dataArray{:, 24};
res_table.W_P7 = dataArray{:, 25};
res_table.W_O1 = dataArray{:, 26};
res_table.W_Oz = dataArray{:, 27};
res_table.W_O2 = dataArray{:, 28};
res_table.W_P4 = dataArray{:, 29};
res_table.W_P8 = dataArray{:, 30};
res_table.W_TP10 = dataArray{:, 31};
res_table.W_CP6 = dataArray{:, 32};
res_table.W_CP2 = dataArray{:, 33};
res_table.W_Cz = dataArray{:, 34};
res_table.W_C4 = dataArray{:, 35};
res_table.W_T8 = dataArray{:, 36};
res_table.W_FT10 = dataArray{:, 37};
res_table.W_FC6 = dataArray{:, 38};
res_table.W_FC2 = dataArray{:, 39};
res_table.W_F4 = dataArray{:, 40};
res_table.W_F8 = dataArray{:, 41};
res_table.W_Fp2 = dataArray{:, 42};
res_table.W_AF7 = dataArray{:, 43};
res_table.W_AF3 = dataArray{:, 44};
res_table.W_F1 = dataArray{:, 45};
res_table.W_F5 = dataArray{:, 46};
res_table.W_FT7 = dataArray{:, 47};
res_table.W_FC3 = dataArray{:, 48};
res_table.W_C1 = dataArray{:, 49};
res_table.W_C5 = dataArray{:, 50};
res_table.W_TP7 = dataArray{:, 51};
res_table.W_CP3 = dataArray{:, 52};
res_table.W_P1 = dataArray{:, 53};
res_table.W_P5 = dataArray{:, 54};
res_table.W_PO7 = dataArray{:, 55};
res_table.W_PO3 = dataArray{:, 56};
res_table.W_POz = dataArray{:, 57};
res_table.W_PO4 = dataArray{:, 58};
res_table.W_PO8 = dataArray{:, 59};
res_table.W_P6 = dataArray{:, 60};
res_table.W_P2 = dataArray{:, 61};
res_table.W_CPz = dataArray{:, 62};
res_table.W_CP4 = dataArray{:, 63};
res_table.W_TP8 = dataArray{:, 64};
res_table.W_C6 = dataArray{:, 65};
res_table.W_C2 = dataArray{:, 66};
res_table.W_FC4 = dataArray{:, 67};
res_table.W_FT8 = dataArray{:, 68};
res_table.W_F6 = dataArray{:, 69};
res_table.W_AF8 = dataArray{:, 70};
res_table.W_AF4 = dataArray{:, 71};
res_table.W_F2 = dataArray{:, 72};
res_table.W_FCz = dataArray{:, 73};
res_table.State = dataArray{:, 74};
res_table.Vig = dataArray{:, 75};
res_table.Pup = dataArray{:, 76};
res_table.pcPup = dataArray{:, 77};
res_table.pcRT = dataArray{:, 78};

res_mat=table2array(res_table);

%%
% clean RTs
res_table.RT(res_table.RT<0.3)=NaN;

res_table.State=categorical(res_table.State);
res_table.SubID=categorical(res_table.SubID);
res_table.StimCat=categorical(res_table.StimCat);
res_table.Task=categorical(res_table.Task);

% res_table.pcPup=ordinal(res_table.pcPup);
res_table.fastRT=double(res_table.pcRT==1);
res_table.slowRT=double(res_table.pcRT==10);

%% LME models
mdl_0= fitlme(res_table,sprintf('Pup~1+(1|SubID)'));
mdl_1= fitlme(res_table,sprintf('Pup~1+Task+(1|SubID)')); %StimCat has no effect on Pup
mdl_2= fitlme(res_table,sprintf('Pup~1+State+Task+(1|SubID)'));
mdl_3= fitlme(res_table,sprintf('Pup~1+State*Task+(1|SubID)'));

mdl_0b= fitlme(res_table,sprintf('pcPup~1+(1|SubID)'));
mdl_1b= fitlme(res_table,sprintf('pcPup~1+State+(1|SubID)')); %StimCat and Task have no effect on pcPup

mdl_0c= fitlme(res_table,sprintf('RT~1+(1|SubID)'));
mdl_1c= fitlme(res_table,sprintf('RT~1+StimCat+(1|SubID)')); %Task have no effect on RT (go)
mdl_2c= fitlme(res_table,sprintf('RT~1+Task+StimCat+(1|SubID)'));
mdl_3c= fitlme(res_table,sprintf('RT~1+Task*StimCat+(1|SubID)'));
mdl_4c= fitlme(res_table,sprintf('RT~1+State*Task*StimCat+(1|SubID)'));

mdl_0d= fitlme(res_table,sprintf('Vig~1+(1|SubID)'));
mdl_1d= fitlme(res_table,sprintf('Vig~1+Task+(1|SubID)')); %StimCat has no effect on Pup
mdl_2d= fitlme(res_table,sprintf('Vig~1+State+Task+(1|SubID)'));
mdl_3d= fitlme(res_table,sprintf('Vig~1+State*Task+(1|SubID)'));

mdl_0e= fitlme(res_table(res_table.StimCat=='0',:),sprintf('fastRT~1+(1|SubID)'));
mdl_1e= fitlme(res_table(res_table.StimCat=='0',:),sprintf('fastRT~1+State+(1|SubID)')); %StimCat has no effect on Pup
% mdl_2e= fitglme(res_table,sprintf('fastRT~1+StimCat+State+(1|SubID)')); %StimCat has no effect on Pup
% mdl_3e= fitglme(res_table,sprintf('fastRT~1+StimCat*State+(1|SubID)')); %StimCat has no effect on Pup

mdl_0f= fitlme(res_table(res_table.StimCat=='0',:),sprintf('slowRT~1+(1|SubID)'));
mdl_1f= fitlme(res_table(res_table.StimCat=='0',:),sprintf('slowRT~1+State+(1|SubID)')); %StimCat has no effect on Pup

% %%
% figure; set(gcf,'Position',[ 440   315   800   440]);
% for nS=1:length(mySubs)
%     subplot(5,5,nS)
%     h1 = raincloud_plot(res_table.RT(res_table.StimCat=='0' & res_table.Task=='1' & res_table.SubID==mySubs(nS)), 'box_on', 1, 'color', [1 0 0], 'alpha', 0.5,...
%         'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
%         'box_col_match', 0);
%     h1{2}.SizeData=72; h1{2}.MarkerFaceAlpha=0.5;
%     h2 = raincloud_plot(res_table.RT(res_table.StimCat=='0' & res_table.Task=='2' & res_table.SubID==mySubs(nS)), 'box_on', 1, 'color',[0 0 1], 'alpha', 0.5,...
%         'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
%     h2{2}.SizeData=72; h2{2}.MarkerFaceAlpha=0.5;
%
% line([1 1]*0.3,ylim,'LineStyle','--','Color','k');
% format_fig;
% set(gca,'XLim', [0.2 1.2])
% end

%%
% findslowRTs=find(res_table.RT<0.3);
% preceedingRT=[]; preceedingRT2=[];
% for k=1:length(findslowRTs)
%     if res_table.TrialN(findslowRTs(k)-1)==res_table.TrialN(findslowRTs(k))-1
%         preceedingRT=[preceedingRT ; res_table.RT(findslowRTs(k)-1)];
%         preceedingRT2=[preceedingRT2 ; [res_table.StimCat(findslowRTs(k)-1) res_table.StimCat(findslowRTs(k))]];
%     end
% end

%%
Subs=unique(res_table.SubID);
distrib_pcRT_group=[];
distrib_pcRT=[];
maxBins=5;

for nS=1:length(Subs)
    for nPr=1:60
        temp = res_table.pcRT(res_table.SubID==Subs(nS) & res_table.ProbeN==(nPr));
        [nout]=hist(temp,1:1:10); nout=nout./sum(nout)*100;
        if maxBins==5
        nout=nout(1:2:10)+nout(2:2:10);
        end
        distrib_pcRT=[distrib_pcRT ; nout];
        temp2 = double(res_table.State(res_table.SubID==Subs(nS) & res_table.ProbeN==(nPr)));
        temp3 = double(res_table.Task(res_table.SubID==Subs(nS) & res_table.ProbeN==(nPr)));
        distrib_pcRT_group=[distrib_pcRT_group ; [nS nPr unique(temp2) unique(temp3)]];
        
        %         p = polyfit(1:10,nout,2);
    end
end

data=[]; data2=[];
for nTask=1:2
    for nState=1:3
        for k=1:maxBins
            data{nTask}{k,nState}=distrib_pcRT(distrib_pcRT_group(:,3)==nState & distrib_pcRT_group(:,4)==nTask,k);
            for nS=1:length(Subs)
                if ~isempty(distrib_pcRT(distrib_pcRT_group(:,1)==nS & distrib_pcRT_group(:,3)==nState & distrib_pcRT_group(:,4)==nTask,k))
                    data2{nTask}{k,nState}(nS)=nanmean(distrib_pcRT(distrib_pcRT_group(:,1)==nS & distrib_pcRT_group(:,3)==nState & distrib_pcRT_group(:,4)==nTask,k));
                else
                    data2{nTask}{k,nState}(nS)=NaN;
                end
            end
        end
        if nState==2
            [h, pVs{nTask}]=ttest2(distrib_pcRT(distrib_pcRT_group(:,3)==nState & distrib_pcRT_group(:,4)==nTask,:),distrib_pcRT(distrib_pcRT_group(:,3)==3 & distrib_pcRT_group(:,4)==nTask,:));
        end
    end
end

%%
figure; set(gcf,'Position',[ 440   315   800   360]);
Tasks={'Face','Digit'};
hs=[];
for nTask=1:2
    subplot(1,2,nTask); format_fig; hold on;
        line([0.5 maxBins+0.5],[1 1]*(100/maxBins),'LineStyle','--','Color',[1 1 1]*0.5,'LineWidth',3)
for nState=2:3
        temp=[];
        for k=1:maxBins
            line([1 1]*k+nState*0.2-0.4,[-1 1]*sem(data{nTask}{k,nState})+nanmean(data{nTask}{k,nState}),'Color',Colors(nState,:),'LineWidth',4);
            temp(k)=nanmean(data{nTask}{k,nState});
        end
        plot((1:maxBins)+nState*0.2-0.4,temp,'Color',Colors(nState,:),'LineWidth',4);
        hs(nState)=scatter((1:maxBins)+nState*0.2-0.4,temp,'MarkerEdgeColor',Colors(nState,:),'MarkerFaceColor','w','SizeData',288,'LineWidth',4,'MarkerFaceAlpha',0.7);
    end
    set(gca,'XTick',1:maxBins,'XTickLabel',{'20^t^h','40^t^h','60^t^h','80^t^h','100^t^h'})
    xlim([0.5 maxBins+0.5])
    xlabel('RT percentile')
    ylabel('% trials')
    title(Tasks{nTask})
    ylim([15 27])
%     legend(hs(1:3),{'MW','MB'})
end
export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Behav_RT_percentile.fig'])
export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Behav_RT_percentile.eps'],'-r 300')
