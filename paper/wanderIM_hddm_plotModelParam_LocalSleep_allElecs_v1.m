%%
clear all;
close all;
run ../localdef_wanderIM.m

addpath(genpath(path_export))
addpath(genpath(path_RainCloudPlot))
addpath(genpath(lscpTools_path))

%% Initialize variables.
np=0;
Elec_Labels=layout.label;
Task_Labels={'Face','Digit'};
Mat_Param=[];
for nE=1:63
    fprintf('... ... Electrode %s\n',Elec_Labels{nE})
    filename=(['/Users/tand0009/Data/WanderIM/hddm/Models/model_W_' Elec_Labels{nE} filesep 'model_stats_W_' Elec_Labels{nE} '.csv']);
    if exist(filename)~=0
        
    else
        fprintf('... ... Electrode %s MISSING\n',Elec_Labels{nE})
        continue;
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
    
    %%%
    var_2factor={'v'};
    var_1factor={'a','t','z'};
    
    mat_LME=[];
    for n=1:length(var_2factor)
        rows=find(~cellfun(@isempty,regexp(model.VarName1,sprintf('^%s_subj*',var_2factor{n}))));
        rownames=strvcat(model.VarName1);
        mat_LME{n}(:,1)=model.mean(rows);
        mat_LME{n}(:,2)=str2num(rownames(rows,8));
        mat_LME{n}(:,3)=str2num(rownames(rows,10));
        mat_LME{n}(:,4)=str2num(rownames(rows,12));
        mat_LME{n}(:,5)=str2num(rownames(rows,end-2:end));
        Mat_Param=[Mat_Param ; [repmat([1 nE],size(mat_LME{n},1),1) mat_LME{n}]];
    end
    tbl_v=array2table(mat_LME{1},'VariableNames',{'v','task','sw','stimtype','subj'});
    tbl_v.stimtype=categorical(tbl_v.stimtype);
    tbl_v.task=categorical(tbl_v.task);
    tbl_v.subj=categorical(tbl_v.subj);
    
    mdl{nE}.v.M0= fitlme(tbl_v,'v~1+stimtype+(1|subj)');
    mdl{nE}.v.M1= fitlme(tbl_v,sprintf('v~1+stimtype+task+(1|subj)'));
    mdl{nE}.v.M2= fitlme(tbl_v,sprintf('v~1+stimtype+task+sw+(1|subj)'));
    mdl{nE}.v.M3= fitlme(tbl_v,sprintf('v~1+stimtype+sw+(1|subj)'));
    
    mdl{nE}.F.vGO.M0= fitlme(tbl_v(tbl_v.stimtype=='1' & tbl_v.task=='1',:),'v~1+(1|subj)');
    mdl{nE}.F.vGO.M1= fitlme(tbl_v(tbl_v.stimtype=='1' & tbl_v.task=='1',:),sprintf('v~1+sw+(1|subj)'));
    
    mdl{nE}.F.vNOGO.M0= fitlme(tbl_v(tbl_v.stimtype=='0' & tbl_v.task=='1',:),'v~1+(1|subj)');
    mdl{nE}.F.vNOGO.M1= fitlme(tbl_v(tbl_v.stimtype=='0' & tbl_v.task=='1',:),sprintf('v~1+sw+(1|subj)'));
    
    mdl{nE}.D.vGO.M0= fitlme(tbl_v(tbl_v.stimtype=='1' & tbl_v.task=='2',:),'v~1+(1|subj)');
    mdl{nE}.D.vGO.M1= fitlme(tbl_v(tbl_v.stimtype=='1' & tbl_v.task=='2',:),sprintf('v~1+sw+(1|subj)'));
    
    mdl{nE}.D.vNOGO.M0= fitlme(tbl_v(tbl_v.stimtype=='0' & tbl_v.task=='2',:),'v~1+(1|subj)');
    mdl{nE}.D.vNOGO.M1= fitlme(tbl_v(tbl_v.stimtype=='0' & tbl_v.task=='2',:),sprintf('v~1+sw+(1|subj)'));
    
    %%%
    mat_LME2=[];
    for n=1:length(var_1factor)
        rows=find(~cellfun(@isempty,regexp(model.VarName1,sprintf('^%s_subj*',var_1factor{n}))));
        rownames=strvcat(model.VarName1);
        mat_LME2{n}(:,1)=model.mean(rows);
        mat_LME2{n}(:,2)=str2num(rownames(rows,8));
        mat_LME2{n}(:,3)=str2num(rownames(rows,10));
        mat_LME2{n}(:,4)=str2num(rownames(rows,13:15));
        Mat_Param=[Mat_Param ; [repmat([2+n nE],size(mat_LME2{n},1),1) mat_LME2{n}(:,1:3) nan(size(mat_LME2{n},1),1) mat_LME2{n}(:,4)]];
    end
    tbl_a=array2table(mat_LME2{1},'VariableNames',{'a','task','sw','subj'});
    tbl_a.subj=categorical(tbl_a.subj);
    tbl_a.task=categorical(tbl_a.task);
    
    tbl_t=array2table(mat_LME2{2},'VariableNames',{'t','task','sw','subj'});
    tbl_t.subj=categorical(tbl_t.subj);
    tbl_t.task=categorical(tbl_t.task);
    
    tbl_z=array2table(mat_LME2{3},'VariableNames',{'z','task','sw','subj'});
    tbl_z.subj=categorical(tbl_z.subj);
    tbl_z.task=categorical(tbl_z.task);
    
    mdl{nE}.a.M0= fitlme(tbl_a,'a~1+(1|subj)');
    mdl{nE}.a.M1= fitlme(tbl_a,sprintf('a~1+task+(1|subj)'));
    mdl{nE}.a.M2= fitlme(tbl_a,sprintf('a~1+task+sw+(1|subj)'));
    mdl{nE}.a.M3= fitlme(tbl_a,sprintf('a~1+sw+(1|subj)'));
    
    mdl{nE}.t.M0= fitlme(tbl_t,'t~1+(1|subj)');
    mdl{nE}.t.M1= fitlme(tbl_t,sprintf('t~1+task+(1|subj)'));
    mdl{nE}.t.M2= fitlme(tbl_t,sprintf('t~1+task+sw+(1|subj)'));
    mdl{nE}.t.M3= fitlme(tbl_t,sprintf('t~1+sw+(1|subj)'));
    
    mdl{nE}.z.M0= fitlme(tbl_z,'z~1+(1|subj)');
    mdl{nE}.z.M1= fitlme(tbl_z,sprintf('z~1+task+(1|subj)'));
    mdl{nE}.z.M2= fitlme(tbl_z,sprintf('z~1+task+sw+(1|subj)'));
    mdl{nE}.z.M3= fitlme(tbl_z,sprintf('z~1+sw+(1|subj)'));
    
    mdl{nE}.F.a.M0= fitlme(tbl_a(tbl_a.task=='1',:),'a~1+(1|subj)');
    mdl{nE}.F.a.M1= fitlme(tbl_a(tbl_a.task=='1',:),sprintf('a~1+sw+(1|subj)'));
    mdl{nE}.D.a.M0= fitlme(tbl_a(tbl_a.task=='2',:),'a~1+(1|subj)');
    mdl{nE}.D.a.M1= fitlme(tbl_a(tbl_a.task=='2',:),sprintf('a~1+sw+(1|subj)'));
    
    mdl{nE}.F.t.M0= fitlme(tbl_t(tbl_t.task=='1',:),'t~1+(1|subj)');
    mdl{nE}.F.t.M1= fitlme(tbl_t(tbl_t.task=='1',:),sprintf('t~1+sw+(1|subj)'));
    mdl{nE}.D.t.M0= fitlme(tbl_t(tbl_t.task=='2',:),'t~1+(1|subj)');
    mdl{nE}.D.t.M1= fitlme(tbl_t(tbl_t.task=='2',:),sprintf('t~1+sw+(1|subj)'));
    
    mdl{nE}.F.z.M0= fitlme(tbl_z(tbl_z.task=='1',:),'z~1+(1|subj)');
    mdl{nE}.F.z.M1= fitlme(tbl_z(tbl_z.task=='1',:),sprintf('z~1+sw+(1|subj)'));
    mdl{nE}.D.z.M0= fitlme(tbl_z(tbl_z.task=='2',:),'z~1+(1|subj)');
    mdl{nE}.D.z.M1= fitlme(tbl_z(tbl_z.task=='2',:),sprintf('z~1+sw+(1|subj)'));
end

%%
Mat_Res=[];
for nTask=1:2
    for nE=1:63
        if isempty(mdl{nE})
            Mat_Res=[Mat_Res ; [nTask nE 1 nan(1,7)] ; ...
                [nTask nE 2 nan(1,7)] ; ...
                [nTask nE 3 nan(1,7)] ; ...
                [nTask nE 4 nan(1,7)] ; ...
                [nTask nE 5 nan(1,7)]];
        else
            if nTask==1
                Mat_Res=[Mat_Res ; [nTask nE 1 double(mdl{nE}.F.vGO.M1.Coefficients(2,2:end))] ; ...
                    [nTask nE 2 double(mdl{nE}.F.vNOGO.M1.Coefficients(2,2:end))] ; ...
                    [nTask nE 3 double(mdl{nE}.F.a.M1.Coefficients(2,2:end))] ; ...
                    [nTask nE 4 double(mdl{nE}.F.t.M1.Coefficients(2,2:end))] ; ...
                    [nTask nE 5 double(mdl{nE}.F.z.M1.Coefficients(2,2:end))]];
            else
                Mat_Res=[Mat_Res ; [nTask nE 1 double(mdl{nE}.D.vGO.M1.Coefficients(2,2:end))] ; ...
                    [nTask nE 2 double(mdl{nE}.D.vNOGO.M1.Coefficients(2,2:end))] ; ...
                    [nTask nE 3 double(mdl{nE}.D.a.M1.Coefficients(2,2:end))] ; ...
                    [nTask nE 4 double(mdl{nE}.D.t.M1.Coefficients(2,2:end))] ; ...
                    [nTask nE 5 double(mdl{nE}.D.z.M1.Coefficients(2,2:end))]];
            end
        end
    end
end
% export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Model_Param_perState.fig'])
% export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/wanderIM/figmaterial/Model_Param_perState.eps'],'-r 300')

%%
addpath(path_fieldtrip);
ft_defaults;
TitlesT={'Faces','Digits'};
TitlesP={'vGO','vNOGO','a','t','z'};
for nTask=1:2
    for nP=1:5
        figure; set(gcf,'Position',[ 562   669   325   316]);
        temp_topo=Mat_Res(Mat_Res(:,1)==nTask & Mat_Res(:,3)==nP,6); %in seconds but for all probes (60) so equivalent of in minutes
        temp_topo([10 21])=NaN;
        temp_pV=Mat_Res(Mat_Res(:,1)==nTask & Mat_Res(:,3)==nP,8);
        temp_pV([10 21])=NaN;
        simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
        colorbar; caxis([-1 1]*6)
        title([TitlesT{nTask} ' - ' TitlesP{nP}])
        format_fig;
        
        load(path_PsychFTlayout);
        if ~isempty(find(temp_pV<fdr(temp_pV,0.05)))
            ft_plot_lay_me(layout, 'chanindx',find(temp_pV<fdr(temp_pV,0.05)),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','no')
        end
    end
end
