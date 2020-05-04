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
for nTask=1:2
    for nE=1:63
        fprintf('... ... Task %s - Electrode %s\n',Task_Labels{nTask},Elec_Labels{nE})
        if nTask==1
            filename = sprintf('/Users/tand0009/WorkGit/projects/inprogress/MW_HDDM/Models/Faces/model_SW/model_W_%s_stats.csv',Elec_Labels{nE});
        else
            filename = sprintf('/Users/tand0009/WorkGit/projects/inprogress/MW_HDDM/Models/Digits/model_SW/model_W_%s_stats.csv',Elec_Labels{nE});
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
    for nE=1:63
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
    colorbar; %caxis([0 9])
    title([TitlesT {nTask} ' - ' TitlesP{nP}])
    format_fig;
    
    load(path_PsychFTlayout);
    if ~isempty(find(temp_pV<fdr(temp_pV,0.05)))
        ft_plot_lay_me(layout, 'chanindx',find(temp_pV<fdr(temp_pV,0.05)),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','no')
    end
end
end
%%
addpath(path_fieldtrip);
ft_defaults;
TitlesP={'vGO','vNOGO','a','t','z'};
TitlesT={'Faces','Digits'};
for nTask=1:2
    for nP=1:5
        figure; set(gcf,'Position',[ 562   669   2*325   316]);
        for nSW=1:2
            subplot(1,2,nSW)
            temp_topo=nan(63,1);
            for nE=1:63
                if nP==1
                    temp_topo(nE)=nanmean(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==nSW-1  & Mat_Param(:,6)==1,4));
                elseif nP==2
                    temp_topo(nE)=nanmean(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==nSW-1  & Mat_Param(:,6)==0,4));
                else
                    temp_topo(nE)=nanmean(Mat_Param(Mat_Param(:,1)==nP & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==nSW-1,4));
                end
            end
            if nP==1
                axesLim=[nanmin(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,6)==1,4)) nanmax(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==1 & Mat_Param(:,6)==1,4))];
            elseif nP==2
                axesLim=[nanmin(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,6)==0,4)) nanmax(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==1 & Mat_Param(:,6)==0,4))];
            else
                axesLim=[nanmin(Mat_Param(Mat_Param(:,1)==nP & Mat_Param(:,2)==nTask,4)) nanmax(Mat_Param(Mat_Param(:,1)==nP & Mat_Param(:,2)==nTask,4))];
            end
            
            temp_topo([10 21])=NaN;
            
            simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
            colorbar; %
            caxis(axesLim)
            if nSW==1
                title([TitlesT{nTask} ' - ' TitlesP{nP} ' - noSW' ])
            else
                title([TitlesT{nTask}  ' - ' TitlesP{nP} ' - SW' ])
            end
            format_fig;
        end
    end
end
%%
addpath(path_fieldtrip);
ft_defaults;
TitlesP={'vGO','vNOGO','a','t','z'};

figure; set(gcf,'Position',[ 13     1   480   977]);
for nTask=1:2
for nP=1:5
    subplot(5,2,nTask+(nP-1)*2)
    temp_topo=nan(63,1);
    for nE=1:63
        if nP==1
            temp_topo(nE)=nanmean(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==1  & Mat_Param(:,6)==1,4))-...
                nanmean(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==0  & Mat_Param(:,6)==1,4));
            fprintf('... np:%s | Elec: %g | %g and %g\n',TitlesP{nP},nE,size(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==1  & Mat_Param(:,6)==1,4),1),size(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==0  & Mat_Param(:,6)==1,4),1))
        elseif nP==2
            temp_topo(nE)=nanmean(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==1  & Mat_Param(:,6)==0,4))-...
                nanmean(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==0  & Mat_Param(:,6)==0,4));
             fprintf('... np:%s | Elec: %g | %g and %g\n',TitlesP{nP},nE,size(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==1  & Mat_Param(:,6)==0,4),1),size(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==0  & Mat_Param(:,6)==0,4),1))
       else
            temp_topo(nE)=nanmean(Mat_Param(Mat_Param(:,1)==nP & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==1,4))-...
                nanmean(Mat_Param(Mat_Param(:,1)==nP & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==0,4));
              fprintf('... np:%s | Elec: %g | %g and %g\n',TitlesP{nP},nE,size(Mat_Param(Mat_Param(:,1)==nP & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==1,4),1),size(Mat_Param(Mat_Param(:,1)==nP & Mat_Param(:,2)==nTask & Mat_Param(:,3)==nE  & Mat_Param(:,5)==0,4),1))
      end
    end
%     if nP==1
%         axesLim=[nanmin(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,6)==1,4)) nanmax(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,6)==1,4))];
%     elseif nP==2
%         axesLim=[nanmin(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,6)==0,4)) nanmax(Mat_Param(Mat_Param(:,1)==1 & Mat_Param(:,2)==nTask & Mat_Param(:,6)==0,4))];
%     else
%         axesLim=[nanmin(Mat_Param(Mat_Param(:,1)==nP & Mat_Param(:,2)==nTask,4)) nanmax(Mat_Param(Mat_Param(:,1)==nP & Mat_Param(:,2)==nTask,4))];
%     end
    temp_topo([10 21])=NaN;
    
    simpleTopoPlot_ft(temp_topo, path_PsychFTlayout,'on',[],0,1);
    colorbar; %
    caxis([-1 1]*nanmax(abs(temp_topo)))
    title([TitlesT{nTask} ' | ' TitlesP{nP} ' | SW - noSW'])
    format_fig;
end
end