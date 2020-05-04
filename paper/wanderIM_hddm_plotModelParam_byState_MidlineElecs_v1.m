%%
clear all;
close all;
run ../localdef_wanderIM.m

addpath(genpath(path_export))
addpath(genpath(path_RainCloudPlot))
addpath(genpath(lscpTools_path))

%% Initialize variables.
np=0;
% Elec_Labels=layout.label;
Task_Labels={'Face','Digit'};
Elec_Labels={'Fz','Cz','Pz','Oz'};
Mat_Param=[];
for nTask=1:2
    for nE=1:length(Elec_Labels)
        fprintf('... ... Task %s - Electrode %s\n',Task_Labels{nTask},Elec_Labels{nE})
        if nTask==1
            filename = sprintf('/Users/tand0009/WorkGit/projects/inprogress/MW_HDDM/Models/Faces/model_SW/model_mid_W_%s_stats.csv',Elec_Labels{nE});
        else
            filename = sprintf('/Users/tand0009/WorkGit/projects/inprogress/MW_HDDM/Models/Digits/model_SW/model_mid_W_%s_stats.csv',Elec_Labels{nE});
        end
        delimiter = ',';
        startRow = 2;
        formatSpec = '%s%f%f%f%f%f%f%f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        fclose(fileID);
        
        model = table;
        model.VarName1 = cellstr(dataArray{:, 1});
        model.mean = dataArray{:, 2};
        model.std = dataArray{:, 3};
        model.q = dataArray{:, 4};
        model.q1 = dataArray{:, 5};
        model.q2 = dataArray{:, 6};
        model.q3 = dataArray{:, 7};
        model.q4 = dataArray{:, 8};
        model.mcerr = dataArray{:, 9};
        clearvars filename delimiter startRow formatSpec fileID dataArray ans;
        
        %%
        var_2factor={'v'};
        var_1factor={'a','t','z'};
        
        mat_LME=[];
        for n=1:length(var_2factor)
            rows=find(~cellfun(@isempty,regexp(model.VarName1,sprintf('^%s_subj*',var_2factor{n}))));
            rownames=strvcat(model.VarName1);
            mat_LME{n}(:,1)=model.mean(rows);
            mat_LME{n}(:,2)=str2num(rownames(rows,8));
            mat_LME{n}(:,3)=str2num(rownames(rows,10));
            mat_LME{n}(:,4)=str2num(rownames(rows,end-2:end));
            Mat_Param=[Mat_Param ; [repmat([1 nTask nE],size(mat_LME{n},1),1) mat_LME{n}]];
        end
        tbl_v=array2table(mat_LME{1},'VariableNames',{'v','sw','stimtype','subj'});
        tbl_v.stimtype=categorical(tbl_v.stimtype);
        tbl_v.subj=categorical(tbl_v.subj);
        
        mdl{nTask,nE}.v.M0= fitlme(tbl_v,'v~1+stimtype+(1|subj)');
        mdl{nTask,nE}.v.M1= fitlme(tbl_v,sprintf('v~1+stimtype+sw+(1|subj)'));
        mdl{nTask,nE}.v.M2= fitlme(tbl_v,sprintf('v~1+stimtype*sw+(1|subj)'));
        
        mdl{nTask,nE}.vGO.M0= fitlme(tbl_v(tbl_v.stimtype=='1',:),'v~1+(1|subj)');
        mdl{nTask,nE}.vGO.M1= fitlme(tbl_v(tbl_v.stimtype=='1',:),sprintf('v~1+sw+(1|subj)'));
        
        mdl{nTask,nE}.vNOGO.M0= fitlme(tbl_v(tbl_v.stimtype=='0',:),'v~1+(1|subj)');
        mdl{nTask,nE}.vNOGO.M1= fitlme(tbl_v(tbl_v.stimtype=='0',:),sprintf('v~1+sw+(1|subj)'));
        
        %%
        mat_LME2=[];
        for n=1:length(var_1factor)
            rows=find(~cellfun(@isempty,regexp(model.VarName1,sprintf('^%s_subj*',var_1factor{n}))));
            rownames=strvcat(model.VarName1);
            mat_LME2{n}(:,1)=model.mean(rows);
            mat_LME2{n}(:,2)=str2num(rownames(rows,8));
            mat_LME2{n}(:,3)=str2num(rownames(rows,11:13));
            Mat_Param=[Mat_Param ; [repmat([2+n nTask nE],size(mat_LME2{n},1),1) mat_LME2{n}(:,1:2) nan(size(mat_LME2{n},1),1) mat_LME2{n}(:,3)]];
        end
        tbl_a=array2table(mat_LME2{1},'VariableNames',{'a','sw','subj'});
        tbl_a.subj=categorical(tbl_a.subj);
        
        tbl_t=array2table(mat_LME2{2},'VariableNames',{'t','sw','subj'});
        tbl_t.subj=categorical(tbl_t.subj);
        
        tbl_z=array2table(mat_LME2{3},'VariableNames',{'z','sw','subj'});
        tbl_z.subj=categorical(tbl_z.subj);
        
        mdl{nTask,nE}.a.M0= fitlme(tbl_a,'a~1+(1|subj)');
        mdl{nTask,nE}.a.M1= fitlme(tbl_a,sprintf('a~1+sw+(1|subj)'));
        
        mdl{nTask,nE}.t.M0= fitlme(tbl_t,'t~1+(1|subj)');
        mdl{nTask,nE}.t.M1= fitlme(tbl_t,sprintf('t~1+sw+(1|subj)'));
        
        mdl{nTask,nE}.z.M0= fitlme(tbl_z,'z~1+(1|subj)');
        mdl{nTask,nE}.z.M1= fitlme(tbl_z,sprintf('z~1+sw+(1|subj)'));
        
    end
end

%%
Mat_Res=[];
for nTask=1:2
    for nE=1:length(Elec_Labels)
        Mat_Res=[Mat_Res ; [nTask nE 1 double(mdl{nTask,nE}.vGO.M1.Coefficients(2,2:end))] ; ...
            [nTask nE 2 double(mdl{nTask,nE}.vNOGO.M1.Coefficients(2,2:end))] ; ...
            [nTask nE 3 double(mdl{nTask,nE}.a.M1.Coefficients(2,2:end))] ; ...
            [nTask nE 4 double(mdl{nTask,nE}.t.M1.Coefficients(2,2:end))] ; ...
            [nTask nE 5 double(mdl{nTask,nE}.z.M1.Coefficients(2,2:end))]];
    end
end
% export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Model_Param_perState.fig'])
% export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Model_Param_perState.eps'],'-r 300')

%%
figure; set(gcf,'Position',[91         522        1149         456]);
for nEl=1:4
    temp_GO=Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nEl & Mat_Param(:,6)==1,3:5);
    subplot(2,4,nEl); hold on;
    simpleBarPlot(1,temp_GO(temp_GO(:,3)==1,2),[.7 .7 .7; 0 0 1],0.9,'k',[],3)
    simpleBarPlot(2,temp_GO(temp_GO(:,3)==0,2),[.7 .7 .7; 1 0 0],0.9,'k',[],3)
    set(gca,'XTick',1:2,'XTickLabel',{'SW','noSW'})
    %     histogram(temp_GO(temp_GO(:,3)==1,2),12)
    %     hold on
    %     histogram(temp_GO(temp_GO(:,3)==0,2),12)
    %    legend({'SW','noSW'})
    title([Elec_Labels{nEl} ' - Go']);
    format_fig;
    
    temp_NOGO=Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nEl & Mat_Param(:,6)==0,3:5);
    subplot(2,4,4+nEl); hold on;
    simpleBarPlot(1,temp_NOGO(temp_NOGO(:,3)==1,2),[.7 .7 .7; 0 0 1],0.9,'k',[],3)
    simpleBarPlot(2,temp_NOGO(temp_NOGO(:,3)==0,2),[.7 .7 .7; 1 0 0],0.9,'k',[],3)
    set(gca,'XTick',1:2,'XTickLabel',{'SW','noSW'})
    %     histogram(temp_NOGO(temp_NOGO(:,3)==1,2),12)
    %     hold on
    %     histogram(temp_NOGO(temp_NOGO(:,3)==0,2),12)
    %    legend({'SW','noSW'})
    title([Elec_Labels{nEl} ' - NoGo']);
    format_fig;
end


%%
figure; set(gcf,'Position',[ 91   198   957   780]);
var_1factor={'a','t','z'};
for nTask=1:2
    for nP=1:3
        subplot(2,3,nP+3*(nTask-1)); hold on;
        
        for nEl=1:4
            temp=Mat_Param(Mat_Param(:,1)==2+nP & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nEl,3:5);
            
            simpleBarPlot(nEl-0.2,temp(temp(:,3)==1,2),[.7 .7 .7; 0 0 1],0.35,'k',[],3);
            simpleBarPlot(nEl+0.2,temp(temp(:,3)==0,2),[.7 .7 .7; 1 0 0],0.35,'k',[],3);
            %     set(gca,'XTick',1:2,'XTickLabel',{'SW','noSW'})
            %     histogram(temp_GO(temp_GO(:,3)==1,2),12)
            %     hold on
            %     histogram(temp_GO(temp_GO(:,3)==0,2),12)
            %    legend({'SW','noSW'})
        end
        set(gca,'XTick',1:4,'XTickLabel',Elec_Labels)
        title(var_1factor{nP});
        format_fig;
        if nP==1
            if nTask==1
                ylabel('Face')
            elseif nTask==2
                ylabel('Digit')
            end
        end
    end
end