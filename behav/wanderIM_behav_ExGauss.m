%%
clear all
close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))
addpath(genpath(exgauss_path))

% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);

state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];

% load([data_path filesep 'CARS_quest'])
%%
exgauss_blocks=[];
exgauss_probes=[];
exgauss_probes2=[];
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
    %     CARS_flag(n)=CARS_bool(CARS_bool(:,1)==str2num(SubID),2);
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
    for nbl=1:6
        temp_RT=test_res(test_res(:,1)==nbl,10)-test_res(test_res(:,1)==nbl,8);
        [X,fVal,exitFlag,solverOutput] = exgauss_fit(temp_RT');
        exgauss_blocks=[exgauss_blocks ; [n nbl unique(test_res(test_res(:,1)==nbl,2)) X]];
    end
    
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
    countpr=0;
    these_RTs=cell(1,4);
    for nbl=1:6
        these_probes=probe_res(probe_res(:,4)==nbl,:);
        these_trials=test_res(test_res(:,1)==nbl,:);
        for npr=1:10
            countpr=countpr+1;
            this_pr_tridx=these_probes(npr,6);
            if npr==1
                last_pr_tridx=0;
            else
                last_pr_tridx=these_probes(npr-1,6);
            end
            temp_testres=these_trials(these_trials(:,4)>this_pr_tridx-21 & these_trials(:,4)<this_pr_tridx,:);
            temp_RT=temp_testres(:,10)-temp_testres(:,8);
            these_RTs{these_probes(npr,32)}=[these_RTs{these_probes(npr,32)} ; temp_RT];
            
            [X,fVal,exitFlag,solverOutput] = exgauss_fit(temp_RT');
            exgauss_probes2=[exgauss_probes2 ; [n these_probes(npr,[1 4:6 31:38]) X]];
            
        end
    end
    for nstate=1:3
        if length(these_RTs{nstate})>100
            tempRT=these_RTs{nstate};
            [X,fVal,exitFlag,solverOutput] = exgauss_fit(temp_RT');
            exgauss_probes=[exgauss_probes ; [n nstate X length(these_RTs{nstate})]];
        end
    end
end

%%
tbl_probes=array2table(exgauss_probes,'VariableNames',{'SubID','State','mu','sigma','tau','ntrial'});
tbl_probes2=array2table(exgauss_probes2,'VariableNames',{'SubID','nPr','nbl','task','ntrial','look','State','ori','awa','int','eng','perf','alert','mu','sigma','tau'});
tbl_test=array2table(exgauss_blocks,'VariableNames',{'SubID','nblock','task','mu','sigma','tau'});

tbl_probes(tbl_probes.State==4,:)=[];
tbl_probes.State=categorical(tbl_probes.State);
tbl_probes.SubID=categorical(tbl_probes.SubID);

tbl_probes2(tbl_probes2.State==4,:)=[];
tbl_probes2.task=categorical(tbl_probes2.task);
tbl_probes2.State=categorical(tbl_probes2.State);
tbl_probes2.SubID=categorical(tbl_probes2.SubID);

tbl_test.task=categorical(tbl_test.task);
tbl_test.SubID=categorical(tbl_test.SubID);

mdl_test_mu=fitlme(tbl_test,'mu~nblock*task+(1|SubID)');
mdl_probes_mu=fitlme(tbl_probes2,'mu~State*task+(1|SubID)');

mdl_test_sigma=fitlme(tbl_test,'sigma~nblock*task+(1|SubID)');
mdl_probes_sigma=fitlme(tbl_probes2,'sigma~State*task+(1|SubID)');

mdl_test_tau=fitlme(tbl_test,'tau~nblock*task+(1|SubID)');
mdl_probes_tau=fitlme(tbl_probes2,'tau~State*task+(1|SubID)');

%%
voi={'mu','sigma','tau'}; %,'vigidx','Widx','NREMidx','REidx','dt','th','al','sp','be'};
figure;
thistbl=tbl_probes2;
for nplot=1:length(voi)
    subplot(1,length(voi),nplot);
    format_fig;
    eval(sprintf('temp=thistbl.%s;',voi{nplot}))
    for nstate=1:3
        for ntask=1:2
            tempplot=temp(thistbl.State==num2str(nstate) & thistbl.task==num2str(ntask));
            if ntask==1
                simpleBarPlot(nstate-0.2,tempplot,Colors(nstate,:),0.35,'k');
            else
                simpleBarPlot(nstate+0.2,tempplot,[1 1 1;Colors(nstate,:)],0.35,'k');
            end
        end
    end
    set(gca,'XTick',1:3,'XTickLabel',{'ON','MW','MB'});
    title([voi{nplot}]);
    xlim([0.2 3.8])
end