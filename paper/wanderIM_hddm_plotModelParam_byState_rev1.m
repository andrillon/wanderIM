%%
clear all;
close all;
run ../localdef_wanderIM.m

addpath(genpath(path_export))
addpath(genpath(path_RainCloudPlot))
addpath(genpath(lscpTools_path))

%% Initialize variables.
np=0;

MarkerT={'d','o'};

filename = '/Users/tand0009/Data/WanderIM/hddm/Models/model_State/model_stats_State.csv';

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
    mat_LME{n}(:,4)=str2num(rownames(rows,12));
    mat_LME{n}(:,5)=str2num(rownames(rows,end-2:end));
end
tbl_v=array2table(mat_LME{1},'VariableNames',{'v','state','task','cond','subj'});
tbl_v.cond=categorical(tbl_v.cond);
tbl_v.state=categorical(tbl_v.state);
tbl_v.task=categorical(tbl_v.task);
tbl_v.subj=categorical(tbl_v.subj);
% tbl_v.state=reordercats(tbl_v.state,[2 1 3]);


mdl.vGO.M0= fitlme(tbl_v(tbl_v.cond=='1',:),'v~1+task+(1|subj)');
mdl.vGO.M1= fitlme(tbl_v(tbl_v.cond=='1',:),sprintf('v~1+task+state+(1|subj)'));

mdl.vNOGO.M0= fitlme(tbl_v(tbl_v.cond=='0',:),'v~1+task+(1|subj)');
mdl.vNOGO.M1= fitlme(tbl_v(tbl_v.cond=='0',:),sprintf('v~1+task+state+(1|subj)'));

mdl.D.vGO.M0= fitlme(tbl_v(tbl_v.cond=='1' & tbl_v.task=='1',:),'v~1+(1|subj)');
mdl.D.vGO.M1= fitlme(tbl_v(tbl_v.cond=='1' & tbl_v.task=='1',:),sprintf('v~1+state+(1|subj)'));

mdl.D.vNOGO.M0= fitlme(tbl_v(tbl_v.cond=='0' & tbl_v.task=='1',:),'v~1+(1|subj)');
mdl.D.vNOGO.M1= fitlme(tbl_v(tbl_v.cond=='0' & tbl_v.task=='1',:),sprintf('v~1+state+(1|subj)'));

mdl.F.vGO.M0= fitlme(tbl_v(tbl_v.cond=='1' & tbl_v.task=='2',:),'v~1+(1|subj)');
mdl.F.vGO.M1= fitlme(tbl_v(tbl_v.cond=='1' & tbl_v.task=='2',:),sprintf('v~1+state+(1|subj)'));

mdl.F.vNOGO.M0= fitlme(tbl_v(tbl_v.cond=='0' & tbl_v.task=='2',:),'v~1+(1|subj)');
mdl.F.vNOGO.M1= fitlme(tbl_v(tbl_v.cond=='0' & tbl_v.task=='2',:),sprintf('v~1+state+(1|subj)'));

uniqueSubs=unique(mat_LME{1}(:,5));
mat_vbias=[];
for nS=1:length(uniqueSubs)
    for nTask=1:2
        for nState=1:3
            tempA=mat_LME{1}(mat_LME{1}(:,5)==uniqueSubs(nS) & mat_LME{1}(:,2)==nState & mat_LME{1}(:,3)==nTask  & mat_LME{1}(:,4)==1,1);
            tempB=mat_LME{1}(mat_LME{1}(:,5)==uniqueSubs(nS) & mat_LME{1}(:,2)==nState & mat_LME{1}(:,3)==nTask  & mat_LME{1}(:,4)==0,1);
            if ~isempty(tempA) && ~isempty(tempB)
                mat_vbias=[mat_vbias ; [tempA+tempB nState nTask uniqueSubs(nS)]];
            end
        end
    end
end
tbl_bv=array2table(mat_vbias,'VariableNames',{'bv','state','task','subj'});
tbl_bv.state=categorical(tbl_bv.state);
tbl_bv.task=categorical(tbl_bv.task);
tbl_bv.subj=categorical(tbl_bv.subj);
% tbl_bv.state=reordercats(tbl_bv.state,[2 1 3]);

mdl.vBias.M0= fitlme(tbl_bv,'bv~1+task+(1|subj)');
mdl.vBias.M1= fitlme(tbl_bv,sprintf('bv~1+task+state+(1|subj)'));

mdl.D.vBias.M0= fitlme(tbl_bv(tbl_bv.task=='1',:),'bv~1+(1|subj)');
mdl.D.vBias.M1= fitlme(tbl_bv(tbl_bv.task=='1',:),sprintf('bv~1+state+(1|subj)'));

mdl.F.vBias.M0= fitlme(tbl_bv(tbl_bv.task=='2',:),'bv~1+(1|subj)');
mdl.F.vBias.M1= fitlme(tbl_bv(tbl_bv.task=='2',:),sprintf('bv~1+state+(1|subj)'));

%%
mat_LME2=[];
for n=1:length(var_1factor)
    rows=find(~cellfun(@isempty,regexp(model.VarName1,sprintf('^%s_subj*',var_1factor{n}))));
    rownames=strvcat(model.VarName1);
    mat_LME2{n}(:,1)=model.mean(rows);
    mat_LME2{n}(:,2)=str2num(rownames(rows,8));
    mat_LME2{n}(:,3)=str2num(rownames(rows,10));
    mat_LME2{n}(:,4)=str2num(rownames(rows,13:15));
end
tbl_a=array2table(mat_LME2{1},'VariableNames',{'a','state','task','subj'});
tbl_a.task=categorical(tbl_a.task);
tbl_a.state=categorical(tbl_a.state);
tbl_a.subj=categorical(tbl_a.subj);
% tbl_a.state=reordercats(tbl_a.state,[2 1 3]);

tbl_t=array2table(mat_LME2{2},'VariableNames',{'t','state','task','subj'});
tbl_t.task=categorical(tbl_t.task);
tbl_t.state=categorical(tbl_t.state);
tbl_t.subj=categorical(tbl_t.subj);
% tbl_t.state=reordercats(tbl_t.state,[2 1 3]);
%         tbl_t.state=reordercats(tbl_t.state,[2 1 3]);

tbl_z=array2table(mat_LME2{3},'VariableNames',{'z','state','task','subj'});
tbl_z.task=categorical(tbl_z.task);
tbl_z.state=categorical(tbl_z.state);
tbl_z.subj=categorical(tbl_z.subj);
% tbl_z.state=reordercats(tbl_z.state,[2 1 3]);

mdl.a.M0= fitlme(tbl_a,'a~1+task+(1|subj)');
mdl.a.M1= fitlme(tbl_a,sprintf('a~1+task+state+(1|subj)'));

mdl.t.M0= fitlme(tbl_t,'t~1+task+(1|subj)');
mdl.t.M1= fitlme(tbl_t,sprintf('t~1+task+state+(1|subj)'));

mdl.z.M0= fitlme(tbl_z,'z~1+task+(1|subj)');
mdl.z.M1= fitlme(tbl_z,sprintf('z~1+task+state+(1|subj)'));

mdl.F.a.M0= fitlme(tbl_a(tbl_a.task=='1',:),sprintf('a~1+(1|subj)'));
mdl.F.t.M0= fitlme(tbl_t(tbl_t.task=='1',:),sprintf('t~1+(1|subj)'));
mdl.F.z.M0= fitlme(tbl_z(tbl_z.task=='1',:),sprintf('z~1+(1|subj)'));

mdl.D.a.M0= fitlme(tbl_a(tbl_a.task=='2',:),sprintf('a~1+(1|subj)'));
mdl.D.t.M0= fitlme(tbl_t(tbl_t.task=='2',:),sprintf('t~1+(1|subj)'));
mdl.D.z.M0= fitlme(tbl_z(tbl_z.task=='2',:),sprintf('z~1+(1|subj)'));

mdl.F.a.M1= fitlme(tbl_a(tbl_a.task=='1',:),sprintf('a~1+state+(1|subj)'));
mdl.F.t.M1= fitlme(tbl_t(tbl_t.task=='1',:),sprintf('t~1+state+(1|subj)'));
mdl.F.z.M1= fitlme(tbl_z(tbl_z.task=='1',:),sprintf('z~1+state+(1|subj)'));

mdl.D.a.M1= fitlme(tbl_a(tbl_a.task=='2',:),sprintf('a~1+state+(1|subj)'));
mdl.D.t.M1= fitlme(tbl_t(tbl_t.task=='2',:),sprintf('t~1+state+(1|subj)'));
mdl.D.z.M1= fitlme(tbl_z(tbl_z.task=='2',:),sprintf('z~1+state+(1|subj)'));


%% rain cloud plot
for nTask=1:2
    np=0;
    figure; set(gcf,'position',[93         182        1040         616])
    state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
    cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];
    for ncond=1:3
        np=np+1;
        subplot(2, 3, np); format_fig;
        for nstate=1:3
            if ncond==1
                temp=mat_LME{1}(mat_LME{1}(:,2)==nstate & mat_LME{1}(:,3)==nTask  & mat_LME{1}(:,4)==0,1);
            elseif ncond==2
                temp=mat_LME{1}(mat_LME{1}(:,2)==nstate & mat_LME{1}(:,3)==nTask  & mat_LME{1}(:,4)==1,1);
            elseif ncond==3
                temp=mat_vbias(mat_vbias(:,2)==nstate & mat_vbias(:,3)==nTask,1);
                
            end
            
<<<<<<< HEAD
<<<<<<< HEAD
            data_toplot{nTask}{ncond,nstate}=temp';
=======
            
>>>>>>> a6f52afce47f471b424cdbfe1c1da28e9133484d
=======
            
>>>>>>> a6f52afce47f471b424cdbfe1c1da28e9133484d
            h = raincloud_plot(temp, 'band_width',0.25, 'box_on', 1, 'color', Colors(nstate,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', 0.4*(nstate-1)+0.1, 'dot_dodge_amount',  0.4*(nstate-1)+0.3,...
                'box_col_match', 1);
            h{2}.Marker=MarkerT{nTask};
            set(h{2},'SizeData',144,'MarkerFaceAlpha',0.7)
        end
        % legend([h1{1} h2{1} h3{1}], {'ON', 'MW', 'MB'});
        % h   = rm_raincloud(rcp_data, state_colours);
        if ncond==1
            if nTask==2
                set(gca, 'Xlim',[-4 0.5],'YLim', [-0.9 1]*1.2);
                title(['D - v - NOGO']);
            else
                set(gca, 'Xlim',[-4 0.5],'YLim', [-0.9 1]*1.5);
                title(['F - v - NOGO']);
            end
        elseif ncond==2
            if nTask==2
                set(gca, 'Xlim',[1 8],'YLim', [-0.9 1]*0.9);
                title(['D - v - GO']);
            else
                set(gca, 'Xlim',[1 8],'YLim', [-0.9 1]*1.05);
                title(['F - v - GO']);
            end
        elseif ncond==3
            if nTask==2
                set(gca, 'Xlim',[1 6],'YLim', [-0.9 1]*0.9);
                title(['D - v bias']);
            else
                set(gca, 'Xlim',[1 6],'YLim', [-0.9 1]*1.05);
                title(['F - v bias']);
            end
        end
        format_fig;
    end
    
    
    
    for ncond=1:3
        np=np+1;
        subplot(2, 3,np); format_fig;
        for nstate=1:3
            
            temp=mat_LME2{ncond}(mat_LME2{ncond}(:,2)==nstate & mat_LME2{ncond}(:,3)==nTask,1);

            
            h = raincloud_plot(temp, 'band_width',[], 'box_on', 1, 'color', Colors(nstate,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', 0.4*(nstate-1)+0.1, 'dot_dodge_amount',  0.4*(nstate-1)+0.3,...
                'box_col_match', 1);
            h{2}.Marker=MarkerT{nTask};
            set(h{2},'SizeData',144,'MarkerFaceAlpha',0.7)
           
<<<<<<< HEAD
<<<<<<< HEAD
                        data_toplot{nTask}{ncond+3,nstate}=temp';

=======
>>>>>>> a6f52afce47f471b424cdbfe1c1da28e9133484d
=======
>>>>>>> a6f52afce47f471b424cdbfe1c1da28e9133484d
        end
        % legend([h1{1} h2{1} h3{1}], {'ON', 'MW', 'MB'});
        % h   = rm_raincloud(rcp_data, state_colours);
        if ncond==1
            if nTask==2
                set(gca, 'XLim', [.9 3], 'YLim', [-0.9 1]*3);
            else
                set(gca, 'XLim', [.9 3], 'YLim', [-0.9 1]*3);
            end
            title(['HDDM - a']);
        elseif ncond==2
            if nTask==2
                set(gca, 'XLim', [0.1 .4], 'YLim', [-0.9 1]*22);
            else
                set(gca, 'XLim', [0.1 .4], 'YLim', [-0.9 1]*19);
            end
            title(['HDDM - t']);
        elseif ncond==3
            if nTask==2
                set(gca, 'XLim', [0.1 0.45], 'YLim', [-0.9 1]*20.5);
                title(['D - z']);
            else
                set(gca, 'XLim', [0.1 0.45], 'YLim', [-0.9 1]*20.5);
                title(['F - z']);
            end
        end
        format_fig;
    end
    
    
    if nTask==1
        export_fig([path_fig filesep  'Model_Param_perState_FACE_v2.fig'])
        export_fig([path_fig filesep  'Model_Param_perState_FACE_v2.eps'],'-r 300')
    else
        export_fig([path_fig filesep  'Model_Param_perState_DGIT_v2.fig'])
        export_fig([path_fig filesep  'Model_Param_perState_DGIT_v2.eps'],'-r 300')
    end
end

%% Model Comparison
tbl_a.task=categorical(tbl_a.task);
tbl_a.subj=categorical(tbl_a.subj);
tbl_a.state=categorical(tbl_a.state);

tbl_t.task=categorical(tbl_t.task);
tbl_t.subj=categorical(tbl_t.subj);
tbl_t.state=categorical(tbl_t.state);

tbl_z.task=categorical(tbl_z.task);
tbl_z.subj=categorical(tbl_z.subj);
tbl_z.state=categorical(tbl_z.state);

tbl_bv.task=categorical(tbl_bv.task);
tbl_bv.subj=categorical(tbl_bv.subj);
tbl_bv.state=categorical(tbl_bv.state);

tbl_v.task=categorical(tbl_v.task);
tbl_v.subj=categorical(tbl_v.subj);
tbl_v.state=categorical(tbl_v.state);
tbl_v.cond=categorical(tbl_v.cond);

mdl.a.M0= fitlme(tbl_a,'a~1+(1|subj)');
mdl.a.M1= fitlme(tbl_a,sprintf('a~1+task+(1|subj)'));
mdl.a.M2= fitlme(tbl_a,sprintf('a~1+task+state+(1|subj)'));
mdl.a.M3= fitlme(tbl_a,sprintf('a~1+task*state+(1|subj)')); % winning model (*)

mdl.t.M0= fitlme(tbl_t,'t~1+(1|subj)');
mdl.t.M1= fitlme(tbl_t,sprintf('t~1+task+(1|subj)'));
mdl.t.M2= fitlme(tbl_t,sprintf('t~1+task+state+(1|subj)')); % winning model (*)
mdl.t.M3= fitlme(tbl_t,sprintf('t~1+task*state+(1|subj)'));

mdl.z.M0= fitlme(tbl_z,'z~1+(1|subj)');
mdl.z.M1= fitlme(tbl_z,sprintf('z~1+task+(1|subj)'));
mdl.z.M2= fitlme(tbl_z,sprintf('z~1+task+state+(1|subj)'));
mdl.z.M3= fitlme(tbl_z,sprintf('z~1+task*state+(1|subj)')); % winning model (*)

mdl.bv.M0= fitlme(tbl_bv,'bv~1+(1|subj)');
mdl.bv.M1= fitlme(tbl_bv,sprintf('bv~1+task+(1|subj)')); % winning model (*)
mdl.bv.M2= fitlme(tbl_bv,sprintf('bv~1+task+state+(1|subj)'));
mdl.bv.M3= fitlme(tbl_bv,sprintf('bv~1+task*state+(1|subj)'));

mdl.vgo.M0= fitlme(tbl_v(tbl_v.cond=='1',:),'v~1+(1|subj)');
mdl.vgo.M1= fitlme(tbl_v(tbl_v.cond=='1',:),sprintf('v~1+task+(1|subj)'));
mdl.vgo.M2= fitlme(tbl_v(tbl_v.cond=='1',:),sprintf('v~1+task+state+(1|subj)')); % winning model (*)
mdl.vgo.M3= fitlme(tbl_v(tbl_v.cond=='1',:),sprintf('v~1+task*state+(1|subj)'));

mdl.vnogo.M0= fitlme(tbl_v(tbl_v.cond=='0',:),'v~1+(1|subj)');
mdl.vnogo.M1= fitlme(tbl_v(tbl_v.cond=='0',:),sprintf('v~1+task+(1|subj)'));
mdl.vnogo.M2= fitlme(tbl_v(tbl_v.cond=='0',:),sprintf('v~1+task+state+(1|subj)')); % winning model (*)
mdl.vnogo.M3= fitlme(tbl_v(tbl_v.cond=='0',:),sprintf('v~1+task*state+(1|subj)'));

tbl_a.state=reordercats(tbl_a.state,[2 1 3]);
tbl_t.state=reordercats(tbl_t.state,[2 1 3]);
tbl_z.state=reordercats(tbl_z.state,[2 1 3]);
tbl_bv.state=reordercats(tbl_bv.state,[2 1 3]);
tbl_v.state=reordercats(tbl_v.state,[2 1 3]);

mdl.a.M3= fitlme(tbl_a,sprintf('a~1+task*state+(1|subj)')); % winning model (*) 
mdl.t.M2= fitlme(tbl_t,sprintf('t~1+task+state+(1|subj)')); % winning model (*) MB vs ON diff
mdl.z.M3= fitlme(tbl_z,sprintf('z~1+task*state+(1|subj)')); % winning model (*) MB vs ON diff
mdl.vgo.M2= fitlme(tbl_v(tbl_v.cond=='1',:),sprintf('v~1+task+state+(1|subj)')); % winning model (*)
mdl.vnogo.M2= fitlme(tbl_v(tbl_v.cond=='0',:),sprintf('v~1+task+state+(1|subj)')); % winning model (*)

%% Stats - Model Comp
comp=compare(mdl.F.vGO.M0,mdl.F.vGO.M1); pcomp(1,1)=double(comp.pValue(2));
comp=compare(mdl.F.vNOGO.M0,mdl.F.vNOGO.M1); pcomp(1,2)=double(comp.pValue(2));
comp=compare(mdl.F.vBias.M0,mdl.F.vBias.M1); pcomp(1,3)=double(comp.pValue(2));
comp=compare(mdl.F.a.M0,mdl.F.a.M1); pcomp(1,4)=double(comp.pValue(2));
comp=compare(mdl.F.t.M0,mdl.F.t.M1); pcomp(1,5)=double(comp.pValue(2));
comp=compare(mdl.F.z.M0,mdl.F.z.M1); pcomp(1,6)=double(comp.pValue(2));

comp=compare(mdl.D.vGO.M0,mdl.D.vGO.M1); pcomp(2,1)=double(comp.pValue(2));
comp=compare(mdl.D.vNOGO.M0,mdl.D.vNOGO.M1); pcomp(2,2)=double(comp.pValue(2));
comp=compare(mdl.D.vBias.M0,mdl.D.vBias.M1); pcomp(2,3)=double(comp.pValue(2));
comp=compare(mdl.D.a.M0,mdl.D.a.M1); pcomp(2,4)=double(comp.pValue(2));
comp=compare(mdl.D.t.M0,mdl.D.t.M1); pcomp(2,5)=double(comp.pValue(2));
comp=compare(mdl.D.z.M0,mdl.D.z.M1); pcomp(2,6)=double(comp.pValue(2));