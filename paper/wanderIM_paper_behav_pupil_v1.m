%%
clear all
close all

run ../localdef_wanderIM.m
addpath(genpath(lscpTools_path))
addpath(genpath(path_RainCloudPlot))

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


%%
% clean RTs
res_table(res_table.RT<0.3,:)=[];
res_mat=table2array(res_table);

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

%%
res_probe=[];
mySubs=unique(res_table.SubID);
for nS=1:length(mySubs)
    for nPr=1:60
        for nSt=0:1
            temp=res_mat(res_table.SubID==mySubs(nS) & res_table.ProbeN==nPr & res_table.StimCat==num2str(nSt),:);
            res_probe=[res_probe ; nanmean(temp,1)];
        end
    end
end

%%
var_table=res_table.Properties.VariableNames;
var_interest={'Vig','pcPup'};
Tasks={'F','D'};
data=[];
for nvar=1:length(var_interest)
    for i=1:2
        for j=1:3
            if strcmp(var_interest{nvar},'RTgo')
                temp_1=res_probe(res_probe(:,match_str(var_table,'Task'))==i & res_probe(:,match_str(var_table,'State'))==j & res_probe(:,match_str(var_table,'StimCat'))==0,match_str(var_table,'RT'));
                temp_2=res_probe(res_probe(:,match_str(var_table,'Task'))==i & res_probe(:,match_str(var_table,'State'))==j & res_probe(:,match_str(var_table,'StimCat'))==0,match_str(var_table,'SubID'));
            elseif strcmp(var_interest{nvar},'RTnogo')
                temp_1=res_probe(res_probe(:,match_str(var_table,'Task'))==i & res_probe(:,match_str(var_table,'State'))==j & res_probe(:,match_str(var_table,'StimCat'))==1,match_str(var_table,'RT'));
                temp_2=res_probe(res_probe(:,match_str(var_table,'Task'))==i & res_probe(:,match_str(var_table,'State'))==j & res_probe(:,match_str(var_table,'StimCat'))==1,match_str(var_table,'SubID'));
            elseif strcmp(var_interest{nvar},'Vig')
                temp_1=abs(res_probe(res_probe(:,match_str(var_table,'Task'))==i & res_probe(:,match_str(var_table,'State'))==j,match_str(var_table,var_interest{nvar}))-6);
                temp_2=res_probe(res_probe(:,match_str(var_table,'Task'))==i & res_probe(:,match_str(var_table,'State'))==j,match_str(var_table,'SubID'));
            else
                temp_1=res_probe(res_probe(:,match_str(var_table,'Task'))==i & res_probe(:,match_str(var_table,'State'))==j,match_str(var_table,var_interest{nvar}));
                temp_2=res_probe(res_probe(:,match_str(var_table,'Task'))==i & res_probe(:,match_str(var_table,'State'))==j,match_str(var_table,'SubID'));
            end
            temp_3=nan(length(mySubs),1);
            for nS=1:length(mySubs)
                if nvar==1
                temp_3(nS)=nanmean(temp_1(temp_2==double(mySubs(nS))));
                else
                temp_3(nS)=nanmedian(temp_1(temp_2==double(mySubs(nS))));
                end
            end
            data{j,i}=temp_3(~isnan(temp_3)); %
            data_perP{nvar}{j,i}=temp_1(~isnan(temp_1));
            data_perS{nvar}{j,i}=temp_3(~isnan(temp_3));
            if strcmp(var_interest{nvar},'Vig')
                vig{j,i}=temp_3;
            elseif strcmp(var_interest{nvar},'pcPup')
                pup{j,i}=temp_3;
            elseif strcmp(var_interest{nvar},'RT')
                RT{j,i}=temp_3;
            end
        end
        %         if strcmp(var_interest{nvar},'RT')
        %             temp_1=res_probe(res_probe(:,match_str(var_table,'Task'))==i,match_str(var_table,var_interest{nvar}));
        %             temp_2=res_probe(res_probe(:,match_str(var_table,'Task'))==i,match_str(var_table,'SubID'));
        %             temp_3=nan(length(mySubs),1);
        %             for nS=1:length(mySubs)
        %                 temp_3(nS)=nanmean(temp_1(temp_2==double(mySubs(nS))));
        %             end
        %             data2{i}=temp_1;
        %         end
    end
    
     
%     figure; set(gcf,'Position',[ 440   315   800   440]); hold on;
%     for ntask=1:2
%                         subplot(1,2,ntask);
%                 h1 = raincloud_plot(data{1,ntask}, 'box_on', 1, 'color', Colors(1,:), 'alpha', 0.5,...
%                     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
%                     'box_col_match', 0); h1{2}.SizeData=72; h1{2}.MarkerFaceAlpha=0.5;
%                 h2 = raincloud_plot(data{2,ntask}, 'box_on', 1, 'color', Colors(2,:), 'alpha', 0.5,...
%                     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
%                 h2{2}.SizeData=72; h2{2}.MarkerFaceAlpha=0.5;
%                 h3 = raincloud_plot(data{3,ntask}, 'box_on', 1, 'color', Colors(3,:), 'alpha', 0.5,...
%                     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
%                 h3{2}.SizeData=72; h3{2}.MarkerFaceAlpha=0.5;

% make figure
figure;set(gcf,'Position',[ 440   315   850   300]); hold on;
for ntask=1:2
                        subplot(1,2,ntask);
h   = rm_raincloud(data(:,ntask), Colors(1,:));
% set(gca, 'XLim', [0 1],'YTick','');
% title(['Figure M9' newline 'Repeated measures raincloud plot']);
        box off
        format_fig;
    
% change one subset to new colour and alter dot size
for i=1:3
    % scatter
    h.s{i, 1}.MarkerFaceColor   = Colors(i,:);
    h.s{i, 1}.MarkerEdgeColor   = [1 1 1]*0.5;
    h.s{i, 1}.MarkerFaceAlpha   = .5;
    h.s{i, 1}.MarkerEdgeAlpha   = .5;
    h.s{i, 1}.SizeData          = 144;
    h.s{i, 1}.LineWidth          = 1;
    h.s{i,1}.YData=h.s{i,1}.YData;%+0.1
    
    % patch
    h.p{i, 1}.FaceColor     = Colors(i,:);
    h.p{i, 1}.LineWidth          = 2;
    
    % scatter mean
    h.m(i,1).MarkerFaceAlpha=1;
    h.m(i,1).MarkerFaceColor=Colors(i,:);
    h.m(i,1).MarkerEdgeColor=[1 1 1]*0.5;
    h.m(i,1).SizeData=576;
    % line mean
    if i<3
        h.l(i,1).Color=[1 1 1]*0.5;
    end
end
format_fig;
if nvar==1
set(gca,'Xlim',[0.9 5.1],'XTick',1:5); %,'Ylim',[-0.75 2.5]);
else
set(gca,'Xlim',[0.9 5.1],'XTick',1:5);
end
title(sprintf('%s - %s',Tasks{ntask},var_interest{nvar}))
end
export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Behav_' var_interest{nvar} '_perState.fig'])
    export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Behav_' var_interest{nvar} '_perState.eps'],'-r 300')
   

%         for nState=1:3
%             line(nState*[1 1],[-1 1]*sem(data{nState,ntask})+nanmean(data{nState,ntask}),'Color',Colors(nState,:),'LineWidth',4);
%           scatter(nState,nanmean(data{nState,ntask}),'MarkerEdgeColor',Colors(nState,:),'MarkerFaceColor','w','SizeData',288,'LineWidth',4,'MarkerFaceAlpha',0.7);
%       end

        %     if ntask==1
        %         set(gca,'XLim', plotLims(nPlot,:), 'YLim', [-3.5 6],'YTick',''); xlabel(plotNames{nPlot});
        %     else
        %         set(gca,'XLim', plotLims(nPlot,:), 'YLim', [-3.5 6],'YTick',''); xlabel(plotNames{nPlot});
        %     end

end

%%
figure;set(gcf,'Position',[ 440   315   850   300]); hold on;
for ntask=1:2
    subplot(1,2,ntask); hold on; format_fig;
    all=[];
    for nstate=1:3
        scatter(data_perS{1}{nstate,ntask},data_perS{2}{nstate,ntask},'MarkerFaceColor','w','MarkerEdgeColor',Colors(nstate,:),'SizeData',72,'LineWidth',2)
        all=[all ; [data_perS{1}{nstate,ntask} data_perS{2}{nstate,ntask}]];
    end
    for nstate=1:3
                scatter(nanmean(data_perS{1}{nstate,ntask}),nanmean(data_perS{2}{nstate,ntask}),'MarkerFaceAlpha',0.5,'MarkerFaceColor',Colors(nstate,:),'MarkerEdgeColor',Colors(nstate,:),'SizeData',544,'LineWidth',3)
    end
     [b,rstats] = robustfit(all(:,1),all(:,2));
            plot(xlim,b(1)+b(2)*xlim,'Color',[1 1 1]*0.5,'LineStyle','--','LineWidth',3);
    xlim([1 5])
    ylim([1 5])
    
end