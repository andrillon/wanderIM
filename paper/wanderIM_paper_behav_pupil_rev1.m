%%
clear all
close all

run ../localdef_wanderIM.m
addpath(genpath(lscpTools_path))
addpath(genpath(path_RainCloudPlot))
<<<<<<< HEAD
<<<<<<< HEAD
addpath(path_export);
=======
>>>>>>> a6f52afce47f471b424cdbfe1c1da28e9133484d
=======
>>>>>>> a6f52afce47f471b424cdbfe1c1da28e9133484d

%% load data
filename = '/Users/tand0009/Data/WanderIM/hddm/HDDM_WIM_localsleep_amp_pup_thrE90P2P_Dec21_v5.txt';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
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
res_table.A_Fp1 = dataArray{:, 74};
res_table.A_Fz = dataArray{:, 75};
res_table.A_F3 = dataArray{:, 76};
res_table.A_F7 = dataArray{:, 77};
res_table.A_FT9 = dataArray{:, 78};
res_table.A_FC5 = dataArray{:, 79};
res_table.A_FC1 = dataArray{:, 80};
res_table.A_C3 = dataArray{:, 81};
res_table.A_T7 = dataArray{:, 82};
res_table.A_TP9 = dataArray{:, 83};
res_table.A_CP5 = dataArray{:, 84};
res_table.A_CP1 = dataArray{:, 85};
res_table.A_Pz = dataArray{:, 86};
res_table.A_P3 = dataArray{:, 87};
res_table.A_P7 = dataArray{:, 88};
res_table.A_O1 = dataArray{:, 89};
res_table.A_Oz = dataArray{:, 90};
res_table.A_O2 = dataArray{:, 91};
res_table.A_P4 = dataArray{:, 92};
res_table.A_P8 = dataArray{:, 93};
res_table.A_TP10 = dataArray{:, 94};
res_table.A_CP6 = dataArray{:, 95};
res_table.A_CP2 = dataArray{:, 96};
res_table.A_Cz = dataArray{:, 97};
res_table.A_C4 = dataArray{:, 98};
res_table.A_T8 = dataArray{:, 99};
res_table.A_FT10 = dataArray{:, 100};
res_table.A_FC6 = dataArray{:, 101};
res_table.A_FC2 = dataArray{:, 102};
res_table.A_F4 = dataArray{:, 103};
res_table.A_F8 = dataArray{:, 104};
res_table.A_Fp2 = dataArray{:, 105};
res_table.A_AF7 = dataArray{:, 106};
res_table.A_AF3 = dataArray{:, 107};
res_table.A_F1 = dataArray{:, 108};
res_table.A_F5 = dataArray{:, 109};
res_table.A_FT7 = dataArray{:, 110};
res_table.A_FC3 = dataArray{:, 111};
res_table.A_C1 = dataArray{:, 112};
res_table.A_C5 = dataArray{:, 113};
res_table.A_TP7 = dataArray{:, 114};
res_table.A_CP3 = dataArray{:, 115};
res_table.A_P1 = dataArray{:, 116};
res_table.A_P5 = dataArray{:, 117};
res_table.A_PO7 = dataArray{:, 118};
res_table.A_PO3 = dataArray{:, 119};
res_table.A_POz = dataArray{:, 120};
res_table.A_PO4 = dataArray{:, 121};
res_table.A_PO8 = dataArray{:, 122};
res_table.A_P6 = dataArray{:, 123};
res_table.A_P2 = dataArray{:, 124};
res_table.A_CPz = dataArray{:, 125};
res_table.A_CP4 = dataArray{:, 126};
res_table.A_TP8 = dataArray{:, 127};
res_table.A_C6 = dataArray{:, 128};
res_table.A_C2 = dataArray{:, 129};
res_table.A_FC4 = dataArray{:, 130};
res_table.A_FT8 = dataArray{:, 131};
res_table.A_F6 = dataArray{:, 132};
res_table.A_AF8 = dataArray{:, 133};
res_table.A_AF4 = dataArray{:, 134};
res_table.A_F2 = dataArray{:, 135};
res_table.A_FCz = dataArray{:, 136};
res_table.State = dataArray{:, 137};
res_table.Vig = dataArray{:, 138};
res_table.Pup = dataArray{:, 139};
res_table.pcPup = dataArray{:, 140};
res_table.pcRT = dataArray{:, 141};
res_table.stimulus = dataArray{:, 142};
res_table.response = dataArray{:, 143};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;


%%
% res_table(res_table.SubID==328,:)=[];
% clean RTs
res_table(res_table.RT<0.3,:)=[];
warning('State-focused analysis, discarding trials further from 20s of Probe onset');
res_table(res_table.DistProbe<-20,:)=[];

res_table.Vig=abs(res_table.Vig-5); % correct Vig scale
res_mat=table2array(res_table);

res_table.State=categorical(res_table.State);
res_table.SubID=categorical(res_table.SubID);
res_table.StimCat=categorical(res_table.StimCat);
res_table.Task=categorical(res_table.Task);

% res_table.pcPup=ordinal(res_table.pcPup);
res_table.fastRT=double(res_table.pcRT==1);
res_table.slowRT=double(res_table.pcRT==10);

%%
uniqueSubID=unique(res_table.SubID);
res_table2=[];
for nS=1:length(uniqueSubID)
    for nP=1:60
        res_table2=[res_table2 ; mean(res_mat(res_table.SubID==uniqueSubID(nS) & res_table.ProbeN==nP,:),1)];
    end
end
res_table2=array2table(res_table2,'VariableNames',res_table.Properties.VariableNames(1:end-2));


res_table2.State=categorical(res_table2.State);
res_table2.SubID=categorical(res_table2.SubID);
res_table2.StimCat=categorical(res_table2.StimCat);
res_table2.Task=categorical(res_table2.Task);


%% LME models
mdl_0c= fitlme(res_table,sprintf('RT~1+(1|SubID)'));
mdl_1c= fitlme(res_table,sprintf('RT~1+StimCat+(1|SubID)')); %Task have no effect on RT (go)
mdl_2c= fitlme(res_table,sprintf('RT~1+Task+StimCat+(1|SubID)'));
mdl_3c= fitlme(res_table,sprintf('RT~1+Task*StimCat+(1|SubID)'));
mdl_4c= fitlme(res_table,sprintf('RT~1+State*Task*StimCat+(1|SubID)'));

mdl_0d= fitlme(res_table,sprintf('Vig~1+(1|SubID)'));
mdl_1d= fitlme(res_table,sprintf('Vig~1+State+(1|SubID)'));
mdl_2d= fitlme(res_table,sprintf('Vig~1+State+Task+(1|SubID)'));
mdl_3d= fitlme(res_table,sprintf('Vig~1+State*Task+(1|SubID)')); %WINING MODEL

mdl_0e= fitlme(res_table2,sprintf('pcPup~1+(1|SubID)'));
mdl_1e= fitlme(res_table2,sprintf('pcPup~1+Task+(1|SubID)'));
mdl_2e= fitlme(res_table2,sprintf('pcPup~1+State+Task+(1|SubID)'));
mdl_3e= fitlme(res_table2,sprintf('pcPup~1+State*Task+(1|SubID)')); %WINING MODEL

res_table2.State=reordercats(res_table2.State,[2 1 3]);
mdl_2eb= fitlme(res_table2,sprintf('pcPup~1+State+Task+(1|SubID)'));
mdl_2db= fitlme(res_table2,sprintf('Vig~1+State+Task+(1|SubID)'));


%%
res_probe=[];
mySubs=unique(res_table.SubID);
for nS=1:length(mySubs)
    for nPr=1:60
        %         for nSt=0:1
        temp=res_mat(res_table.SubID==mySubs(nS) & res_table.ProbeN==nPr,:);
        res_probe=[res_probe ; nanmean(temp,1)];
        %         end
    end
end

%%
var_table=res_table.Properties.VariableNames;
var_interest={'Vig','pcPup'};
Tasks={'F','D'};
data=[];
for nvar=1:length(var_interest)
    for j=1:3
        for i=1:2
            if strcmp(var_interest{nvar},'RTgo')
                temp_1=res_probe(res_probe(:,match_str(var_table,'Task'))==(i) & res_probe(:,match_str(var_table,'State'))==(j) & res_probe(:,match_str(var_table,'StimCat'))==0,match_str(var_table,'RT'));
                temp_2=res_probe(res_probe(:,match_str(var_table,'Task'))==(i) & res_probe(:,match_str(var_table,'State'))==(j) & res_probe(:,match_str(var_table,'StimCat'))==0,match_str(var_table,'SubID'));
            elseif strcmp(var_interest{nvar},'RTnogo')
                temp_1=res_probe(res_probe(:,match_str(var_table,'Task'))==(i) & res_probe(:,match_str(var_table,'State'))==(j) & res_probe(:,match_str(var_table,'StimCat'))==1,match_str(var_table,'RT'));
                temp_2=res_probe(res_probe(:,match_str(var_table,'Task'))==(i) & res_probe(:,match_str(var_table,'State'))==(j) & res_probe(:,match_str(var_table,'StimCat'))==1,match_str(var_table,'SubID'));
            elseif strcmp(var_interest{nvar},'Vig')
                temp_1=res_probe(res_probe(:,match_str(var_table,'Task'))==(i) & res_probe(:,match_str(var_table,'State'))==(j),match_str(var_table,var_interest{nvar}));
                temp_2=res_probe(res_probe(:,match_str(var_table,'Task'))==(i) & res_probe(:,match_str(var_table,'State'))==(j),match_str(var_table,'SubID'));
            else
                temp_1=res_probe(res_probe(:,match_str(var_table,'Task'))==(i) & res_probe(:,match_str(var_table,'State'))==(j),match_str(var_table,var_interest{nvar}));
                temp_2=res_probe(res_probe(:,match_str(var_table,'Task'))==(i) & res_probe(:,match_str(var_table,'State'))==(j),match_str(var_table,'SubID'));
            end
            temp_3=nan(length(mySubs),1);
            temp_4=nan(length(mySubs),1);
            for nS=1:length(mySubs)
                if nvar==1
                    temp_3(nS)=nanmean(temp_1(temp_2==str2num(char(mySubs(nS)))));
                else
                    temp_3(nS)=nanmedian(temp_1(temp_2==str2num(char(mySubs(nS)))));
                end
                temp_4(nS)=sum(temp_2==str2num(char(mySubs(nS))) & ~isnan(temp_1));
            end
            data{j,i}=temp_3(~isnan(temp_3)); %
            ndata{j,i}=temp_4(~isnan(temp_3));
            ndata2{nvar}{j,i}=temp_4;
            data_perP{nvar}{j,i}=temp_1(~isnan(temp_1));
            datattest{nvar}{j,i}=temp_3; datattest{nvar}{j,i}(ndata2{1}{j,i}<3)=NaN;
            data_perS{nvar}{j,i}=temp_3(~isnan(temp_3));
            ndata_perS{nvar}{j,i}=ndata{j,i};
            ndata{j,i}=50+100*minmax(ndata{j,i});%
            if strcmp(var_interest{nvar},'Vig')
                vig{j,i}=temp_3;
            elseif strcmp(var_interest{nvar},'pcPup')
                pup{j,i}=temp_3;
            elseif strcmp(var_interest{nvar},'RT')
                RT{j,i}=temp_3;
            end
        end
        datattest{nvar}{j,3}=nanmean([datattest{nvar}{j,1} datattest{nvar}{j,2}],2);
        %         datattest{nvar}{j,3}(ndata2{j,1}<3 | ndata2{j,2}<3)=NaN;
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
    % figure;set(gcf,'Position',[ 440   315   850   300]); hold on;
    figure;set(gcf,'Position',[ 440   315   420   420]); hold on;
    h   = rm_raincloud(data, Colors(1:2,:),0, 'ks', [],ndata);
    box off
    format_fig;
    for ntask=1:2
        
        %                         subplot(1,2,ntask);
        % set(gca, 'XLim', [0 1],'YTick','');
        % title(['Figure M9' newline 'Repeated measures raincloud plot']);
        
        % change one subset to new colour and alter dot size
        for i=1:3
            % scatter
            h.s{i, 2}.MarkerFaceColor   = Colors(i,:);
            h.s{i, 2}.MarkerEdgeColor   = [1 1 1]*0.5;
            %h.s{i, 2}.SizeData          = 100;%50+50*minmax(ndata{i,2})';
            h.s{i, 2}.LineWidth          = 2;
            h.s{i,2}.YData=h.s{i,2}.YData-0.05;
            %             h.s{i,2}.ZData=44.*(ndata{i,2})';
            
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
                h.l(i,2).Color=[0 0 0];%[190 90 255]./255;
            end
        end
        for i=1:3
            h.s{i, 1}.Marker         = 'd';
            h.s{i, 1}.MarkerEdgeColor   = [1 1 1]*0.5;
            h.s{i, 1}.MarkerFaceColor   = Colors(i,:);
            %     h.s{i, 1}.SizeData          = 100;
            h.s{i, 1}.LineWidth          = 2;
            h.s{i,1}.YData=h.s{i,1}.YData+0.05;
            
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
                h.l(i,1).Color=[0 0 0];
                h.l(i,1).LineStyle='-';
            end
        end
        format_fig;
        if nvar==1
            set(gca,'Xlim',[0.9 4.1],'XTick',1:4); %,'Ylim',[-0.75 2.5]);
        else
            set(gca,'Xlim',[0.9 5.1],'XTick',1:5);
            ylim([-1.1 3.3])
        end
        title(sprintf('%s - %s','B',var_interest{nvar}))
    end
    export_fig([path_fig filesep 'Behav_' var_interest{nvar} '_bothTaks_perState.fig'])
    export_fig([path_fig filesep 'Behav_' var_interest{nvar} '_bothTaks_perState.eps'],'-r 300')
    
    
    %         for nState=1:3
    %             line(nState*[1 1],[-1 1]*sem(data{nState,ntask})+nanmean(data{nState,ntask}),'Color',Colors(nState,:),'LineWidth',4);
    %           scatter(nState,nanmean(data{nState,ntask}),'MarkerEdgeColor',Colors(nState,:),'MarkerFaceColor','w','SizeData',288,'LineWidth',4,'MarkerFaceAlpha',0.7);
    %       end
    
    %     if ntask==1
    %         set(gca,'XLim', plotLims(nPlot,:), 'YLim', [-3.5 6],'YTick',''); xlabel(plotNames{nPlot});
    %     else
    %         set(gca,'XLim', plotLims(nPlot,:), 'YLim', [-3.5 6],'YTick',''); xlabel(plotNames{nPlot});
    %     end
    fprintf('%s - %s\n','B',var_interest{nvar})
    [h, pV, ~, stats]=ttest([datattest{nvar}{2,3}],[datattest{nvar}{1,3}]);
    fprintf('post-hoc ttests (av task): MW vs ON t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)
    [h, pV, ~, stats]=ttest([datattest{nvar}{3,3}],[datattest{nvar}{1,3}]);
    fprintf('post-hoc ttests (av task): MB vs ON t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)
    [h, pV, ~, stats]=ttest([datattest{nvar}{3,3}],[datattest{nvar}{2,3}]);
    fprintf('post-hoc ttests (av task): MB vs MW t(%g)=%2.4f (p=%1.5f)\n',stats.df,stats.tstat,pV)
    
end

%%
% figure;set(gcf,'Position',[ 440   315   850   300]); hold on;
figure;set(gcf,'Position',[ 440   315   420   420]); hold on;
TMarker={'d','o'};
for ntask=1:2
    %     subplot(1,2,ntask);
    hold on; format_fig;
    all=[];
    for nstate=1:3
        scatter(data_perS{1}{nstate,ntask},data_perS{2}{nstate,ntask},ndata_perS{1}{nstate,ntask},'Marker',TMarker{ntask},'MarkerFaceColor',Colors(nstate,:),'MarkerFaceAlpha',0.3,'MarkerEdgeColor',Colors(nstate,:),'LineWidth',2)
        all=[all ; [data_perS{1}{nstate,ntask} data_perS{2}{nstate,ntask}]];
        
    end
end
for ntask=1:2
    for nstate=1:3
        scatter(sum(data_perS{1}{nstate,ntask}.*ndata_perS{1}{nstate,ntask})/sum(ndata_perS{1}{nstate,ntask}),sum(data_perS{2}{nstate,ntask}.*ndata_perS{2}{nstate,ntask})/sum(ndata_perS{2}{nstate,ntask}),'Marker',TMarker{ntask},'MarkerFaceAlpha',0.8,'MarkerFaceColor',Colors(nstate,:),'MarkerEdgeColor',[1 1 1]*0.5,'SizeData',544,'LineWidth',3)
    end
    [b,rstats] = robustfit(all(:,1),all(:,2));
    %      b=polyfit(all(:,1),all(:,2),1); b=flip(b);
    
    xlim([1 5])
    ylim([1 5])
    %     plot(xlim,b(1)+b(2)*xlim,'Color',[1 1 1]*0.5,'LineStyle','--','LineWidth',3);
    line(1:5,1:5,'Color',[1 1 1]*0.5,'LineStyle',':','LineWidth',2);
    set(gca,'XTick',1:5,'YTick',1:5)
    xlabel(var_interest{1})
    ylabel(var_interest{2})
    axis equal
end
xlim([1 5])
ylim([1 5])
%export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Behav_VigXPup_bothTaks_perState.fig'])
%export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Behav_VigXPupbothTaks_perState.eps'],'-r 300')

[rall, pVall]=corr(all(:,1),all(:,2),'Type','Pearson','Rows','Pairwise');
% %% RT * Corr * pcPup
% rtbins=10;
% figure; hold on;
% cb=colormap('hot'); cb=flipud(cb);
% mySubs=unique(res_table.SubID);
% for ntask=1%:2
%     %     subplot(1,2,ntask); hold on;
%     clear mean_pup sem_pup mean_corr
%     for nS=1:length(mySubs)
%         %         tempRT=res_table.RT(res_table.Task==num2str(ntask) & res_table.SubID==mySubs(nS));
%         %         tempPup=res_table.pcPup(res_table.Task==num2str(ntask) & res_table.SubID==mySubs(nS));
%         %         tempPerf=res_table.Perf(res_table.Task==num2str(ntask) & res_table.SubID==mySubs(nS));
%         tempRT=res_table.RT(res_table.SubID==mySubs(nS));
%         tempPup=res_table.pcPup(res_table.SubID==mySubs(nS));
%         tempPerf=res_table.Perf(res_table.SubID==mySubs(nS));
%
%         boundaries=prctile(tempRT,0:rtbins:100);
%         temppcRT=nan(size(tempRT,1),1);
%         for nPer=1:length(0:rtbins:100)-2
%             temppcRT(tempRT>=boundaries(nPer) & tempRT<boundaries(nPer+1))=nPer;
%         end
%         temppcRT(tempRT>=boundaries(length(0:rtbins:100)-1))=length(0:rtbins:100)-1;
%
%         for nPer=1:rtbins
%             mean_pup(nPer,nS)=nanmean(tempPup(temppcRT==nPer));
%             mean_corr(nPer,nS)=nanmean(tempPerf(temppcRT==nPer));
%         end
%     end
%
%     for k=1:10
%
%         sem_pup(k)=sem(mean_pup(k,:));
%         thisC=cb(round(64*((mean(mean_corr(k,:))-0.9)/(1-0.9)/2+0.5)),:);
%         line([1 1]*k,[-1 1]*sem_pup(k)+mean(mean_pup(k,:)),'Color',thisC,'LineWidth',3);
%         scatter(k,mean(mean_pup(k,:)),'MarkerFaceColor',thisC,'MarkerEdgeColor',thisC,'LineWidth',3,'SizeData',2^8);
%     end
% end
% %  end
%
% %% RT * Corr * pcPup
% rtbins=10;
% figure; hold on;
% cb=colormap('hot'); cb=flipud(cb);
% mySubs=unique(res_table.SubID);
% for nSt=1:3
%     subplot(1,3,nSt); hold on;
%     clear mean_pup sem_pup mean_corr sem_rt mean_rt
%     for ntask=1:2
%         %     subplot(1,2,ntask); hold on;
%         for nS=1:length(mySubs)
%             tempRT=res_table.RT(res_table.Task==num2str(ntask) & res_table.SubID==mySubs(nS) & res_table.State==num2str(nSt));
%             tempPup=res_table.Pup(res_table.Task==num2str(ntask) & res_table.SubID==mySubs(nS) & res_table.State==num2str(nSt));
%             tempPerf=res_table.Perf(res_table.Task==num2str(ntask) & res_table.SubID==mySubs(nS) & res_table.State==num2str(nSt));
%             %         tempRT=res_table.RT(res_table.SubID==mySubs(nS));
%             %         tempPup=res_table.Pup(res_table.SubID==mySubs(nS));
%             %         tempPerf=res_table.Perf(res_table.SubID==mySubs(nS));
%             %
%             boundaries=prctile(tempRT,0:rtbins:100);
%             temppcRT=nan(size(tempRT,1),1);
%             for nPer=1:length(0:rtbins:100)-2
%                 temppcRT(tempRT>=boundaries(nPer) & tempRT<boundaries(nPer+1))=nPer;
%             end
%             temppcRT(tempRT>=boundaries(length(0:rtbins:100)-1))=length(0:rtbins:100)-1;
%
%             boundaries=prctile(tempPup,0:rtbins:100);
%             temppcPup=nan(size(tempPup,1),1);
%             for nPer=1:length(0:rtbins:100)-2
%                 temppcPup(tempPup>=boundaries(nPer) & tempPup<boundaries(nPer+1))=nPer;
%             end
%             temppcPup(tempPup>=boundaries(length(0:rtbins:100)-1))=length(0:rtbins:100)-1;
%
%
%             for nPer=1:rtbins
%                 mean_pup(ntask,nPer,nS)=nanmean(temppcPup(temppcRT==nPer));
%                 mean_corr2(ntask,nPer,nS)=nanmean(tempPerf(temppcRT==nPer));
%                 mean_rt(ntask,nPer,nS)=nanmean(temppcRT(temppcPup==nPer));
%                 mean_corr(ntask,nPer,nS)=nanmean(tempPerf(temppcPup==nPer));
%             end
%         end
%     end
%
%
%     mean_rt=squeeze(nanmean(mean_rt,1));
%     mean_corr=squeeze(nanmean(mean_corr2,1));
%     mean_pup=squeeze(nanmean(mean_pup,1));
%
%     for k=1:rtbins
%
%         sem_pup(k)=sem(mean_pup(k,:));
%         thisC=cb(round(64*((nanmean(mean_corr(k,:))-nanmin(nanmean(mean_corr,2)))/(nanmax(nanmean(mean_corr,2))-nanmin(nanmean(mean_corr,2)))/2+0.5)),:);
%         line([1 1]*k,[-1 1]*sem_pup(k)+nanmean(mean_pup(k,:)),'Color',thisC,'LineWidth',3);
%         scatter(k,nanmean(mean_pup(k,:)),'MarkerFaceColor',thisC,'MarkerEdgeColor',thisC,'LineWidth',3,'SizeData',2^8);
%     end
% end

%%
mdl_3e= fitlme(res_table2,sprintf('All_Waves~1+Task+ProbeWithin+BlockN+(1|SubID)')); %WINING MODEL

%%
figure;
simpleCorPlotsetbin(res_probe(:,138),res_probe(:,140),1:4,{'o','k',[1 1 1]*0.7,244,5});
set(gca,'XTick',1:4);
xlim([0.5 4.5])
% ylim([0.5 5.5])
ylabel('Vigilance Score')
xlabel('Pupil Size')
format_fig;

res_probe_table=array2table(res_probe(:,[1 2 6 138 140]), 'VariableNames',{'SubID','Block','Task','Vig','Pup'});
res_probe_table.SubID=categorical(res_probe_table.SubID);
res_probe_table.Task=categorical(res_probe_table.Task);

mdl_pupvig0=fitlme(res_probe_table,sprintf('Pup~1+Task+(1|SubID)'));
mdl_pupvig1=fitlme(res_probe_table,sprintf('Pup~1+Task+Vig+(1|SubID)'));
[r, pV]=corr(res_probe_table.Vig,res_probe_table.Pup,'type','Spearman','rows','pairwise');


uniqueS=unique(res_probe_table.SubID);
r_pupvig=cell(1,2);
for nT=1:2
    for nS=1:length(uniqueS)
        vig=res_probe_table.Vig(res_probe_table.SubID==uniqueS(nS) & res_probe_table.Task==num2str(nT));
        pup=res_probe_table.Pup(res_probe_table.SubID==uniqueS(nS) & res_probe_table.Task==num2str(nT));
        [r_pupvig{nT}(nS), ~]=corr(vig,pup,'type','Spearman','rows','pairwise');
    end
end

figure;
simpleDotPlot(1,r_pupvig{1},366,[138,43,226]./255,1,'k','o',[],3,1);
simpleDotPlot(2,r_pupvig{2},366,[138,43,226]./255,1,'k','d',[],3,1);
xlim([0.5 2.5])
ylim([-1 1])
line(xlim,[0 0],'Color',[1 1 1]*0.7,'LineStyle','--');
format_fig;
set(gca,'XTick',1:2,'XTickLabel',{'Face','Digit'});
ylabel('Spearman Coefficients');