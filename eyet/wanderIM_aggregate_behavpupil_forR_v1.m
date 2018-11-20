%%
clear all
close all

run ../localdef_wanderIM

addpath(genpath(lscpTools_path))

% data_path='/Users/Thomas/temp_data/WanderIM/behavonly/';
data_path=[root_path filesep 'behav/'];
files=dir([data_path filesep 'wanderIM_behavres_s3*.mat']);
eyet_path=[root_path filesep 'eyetracker'];

state_colours=[0 146 146 ; 182 109 255; 219 209 0]/255;
cond_colours=[0.9 0.55 0.2 ; 0.2 0.55 0.9];

% load([data_path filesep 'CARS_quest'])
%%
all_test_res=[];
all_probes_mat=[];
RTs=[];
for n=1:length(files)
    % load Behaviour
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s\n',SubID)
    
    % load eye-tracker data
    savename=['wanderIM_eyelink_S' SubID];
    load([eyet_path filesep savename])
    fprintf('... %s loaded\n',SubID)
    Fs=EL_headers.Fs;
    
    % Prepare ET data
    myEventsType=EL_events.Events.type;
    myEventsTime=EL_events.Events.time;
    data_time=EL_data.time;
    rawdata_pupil=EL_data.pupilSize;
    data_pupil=get_cleanpupil(rawdata_pupil,data_time,EL_events);
    
    % Test
    RTs=[RTs; test_res(:,10)-test_res(:,8)];
    temp_perf=min(test_res(:,11:12),[],2);
    temp_cat=(test_res(:,5)==test_res(:,6));
    temp_RT=(test_res(:,10)-test_res(:,8));
    fprintf('%3.0f%%\n',0)
    thisPup=nan(size(test_res,1),1);
    for nTr=1:size(test_res,1)
        fprintf('\b\b\b\b\b%3.0f%%\n',nTr/size(test_res,1)*100)
        thisBl=test_res(nTr,1);
        thisTr=test_res(nTr,4);
        thisEv=match_str(myEventsType,sprintf('B%g_T%g',thisBl,thisTr));
        thisEvTime=myEventsTime(thisEv);
        thisEvIdx=find(data_time==thisEvTime);
        thisEvIdx2=find(data_time==thisEvTime+Fs);
        %         [~,thisEvIdx]=findclosest(data_time,thisEvTime);
        %         [~,thisEvIdx2]=findclosest(data_time,thisEvTime+Fs);
        %         if nTr==1
        %             thisPup=nan(size(test_res,1),length(thisEvIdx(1):thisEvIdx2(1)));
        %         end
        thisPup(nTr)=nanmean(data_pupil(thisEvIdx(1):thisEvIdx2(1)));
    end
    all_test_headers={'SubID','nBlock','Task','nTrial','StimID','Corr','RT','TrCat','Pup'};
    all_test_res=[all_test_res ; [str2num(SubID)*ones(size(test_res,1),1) test_res(:,[1 2 4 5]) temp_perf temp_RT temp_cat thisPup]];
    
    % Probes
    countpr=0;
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
            probe_details(1)=these_probes(npr,31); % look
            probe_details(2)=these_probes(npr,32); % look
            probe_details(3)=these_probes(npr,33); % look
            probe_details(4)=these_probes(npr,34); % look
            probe_details(5)=these_probes(npr,35); % look
            probe_details(6)=these_probes(npr,36); % look
            probe_details(7)=these_probes(npr,37); % look
            probe_details(8)=these_probes(npr,38); % look
            
            temp_testres=these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,:);
            %             temp_testres(1:round(size(temp_testres,1)/2),:)=[];
            tcorr_perf=min(temp_testres(:,11:12),[],2);
            tcorr_catTr=(temp_testres(:,5)==temp_testres(:,6)); %0 for go, 1 for nogo
            tcorr_RT=temp_testres(:,10)-temp_testres(:,8);
            tcorr_distPr=nan(1,length(tcorr_catTr));
            tcorr_distPr(tcorr_catTr==1)=-sum(tcorr_catTr==1):1:-1;
            tcorr_distPr(tcorr_catTr==0)=-sum(tcorr_catTr==0):1:-1;
            
            % take Pupil over 20s;
            temp_trIdx=these_trials(these_trials(:,4)>last_pr_tridx & these_trials(:,4)<this_pr_tridx,4);
            thisPup=nan(length(temp_trIdx),1);
            for nTr=1:length(temp_trIdx)
                thisBl=nbl;
                thisTr=temp_trIdx(nTr);
                thisEv=match_str(myEventsType,sprintf('B%g_T%g',thisBl,thisTr));
                thisEvTime=myEventsTime(thisEv);
                thisEvIdx=find(data_time==thisEvTime);
                thisEvIdx2=find(data_time==thisEvTime+Fs);
                thisPup(nTr)=nanmean(data_pupil(thisEvIdx(1):thisEvIdx2(1)));
            end
            
            % Colum order:
            all_probes_headers={'SubID','nBlock','nProbe','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe','Pup'};
            all_probes_mat=[all_probes_mat ; [repmat([str2num(SubID) nbl countpr these_probes(npr,5) this_pr_tridx probe_details],size(temp_testres,1),1) tcorr_perf tcorr_RT tcorr_catTr tcorr_distPr' thisPup]];
        end
    end
end

%% transform into tables and export
tbl_test=array2table(all_test_res,'VariableNames',all_test_headers);
% 'SubID','nBlock','Task','nTrial','StimID','Corr','RT','TrCat'
tbl_test.SubID=categorical(tbl_test.SubID);
tbl_test.Task=categorical(tbl_test.Task);
tbl_test.StimID=categorical(tbl_test.StimID);
tbl_test.TrCat=categorical(tbl_test.TrCat);

writetable(tbl_test,[data_path filesep 'WanderIM_BehavPupil_TestResults.txt']);


tbl_probe=array2table(all_probes_mat,'VariableNames',all_probes_headers);
% 'SubID','nBlock','Task','nTrial','Look','State','Orig','Awa','Int','Eng','Perf','Vig','Corr','RT','TrCat','DistProbe'
tbl_probe.SubID=categorical(tbl_probe.SubID);
tbl_probe.Task=categorical(tbl_probe.Task);
tbl_probe.Look=categorical(tbl_probe.Look);
tbl_probe.State=categorical(tbl_probe.State);
tbl_probe.Orig=categorical(tbl_probe.Orig);
tbl_probe.TrCat=categorical(tbl_probe.TrCat);

writetable(tbl_probe,[data_path filesep 'WanderIM_BehavPupil_ProbeResults.txt']);


%%
Prop={'o','k','k',128,2};
X=tbl_test.RT;
Y=tbl_test.Pup;
Z=tbl_test.Corr;
bins=0.05:0.05:1.2;

disc=isnan(X) | isnan(Y);
X(disc)=[];
Y(disc)=[];
Z(disc)=[];
clear Xbin Ybin Zbin

figure; format_fig; hold on;
Ybin(1)=mean(Y(X<=bins(1)));
Xbin(1)=mean(X(X<=bins(1)));
Zbin(1)=mean(Z(X<=bins(1)));
Ybin_sem(1)=sem(Y(X<=bins(1)));
line([1 1]*Xbin(1),[-1 1]*Ybin_sem(1)+Ybin(1),'Color',Prop{2},'LineWidth',2)
scatter(Xbin(1),Ybin(1),[],Zbin(1),'filled','SizeData',Prop{4})
for nbin=2:length(bins)-2
    Ybin(nbin)=mean(Y(X>bins(nbin-1) & X<=bins(nbin)));
    Ybin_sem(nbin)=real(sem(Y(X>bins(nbin-1) & X<=bins(nbin))));
    Xbin(nbin)=mean(X(X>bins(nbin-1) & X<=bins(nbin)));
    Xbin_sem(nbin)=real(sem(X(X>bins(nbin-1) & X<=bins(nbin))));
    Zbin(nbin)=mean(Z(X>bins(nbin-1) & X<=bins(nbin)));
    
    line([1 1]*Xbin(nbin),[-1 1]*Ybin_sem(nbin)+Ybin(nbin),'Color',Prop{2},'LineWidth',2)
    line([-1 1]*Xbin_sem(nbin)+Xbin(nbin),[1 1]*Ybin(nbin),'Color',Prop{2},'LineWidth',2)
    scatter(Xbin(nbin),Ybin(nbin),[],Zbin(nbin),'filled','SizeData',Prop{4});
end
cmap=colormap('hot');
cmap=flipud(cmap);
colormap(cmap);
colorbar;
caxis([0.5 1])
plot(Xbin,Ybin,'Color','k','LineWidth',2)
xlabel('RT')
ylabel('Pupil diameter')

%% 2D version

Prop={'o','k','k',128,2};
X=tbl_test.RT;
Y=tbl_test.Pup;
Z=tbl_test.Corr;
binsX=0.05:0.05:1.2; binsX=[binsX Inf];
binsY=400:100:2000; binsY=[binsY Inf];

disc=isnan(X) | isnan(Y);
X(disc)=[];
Y(disc)=[];
Z(disc)=[];

Zbin=zeros(length(binsX)-1,length(binsY)-1);
for nx=1:length(binsX)-1
    for ny=1:length(binsY)-1
        %         if nx==length(binsX) && ny==length(binsX)
        %             Zbin(nx,ny)=Z(X>+binsX(nx) & Y>=binsY(ny));
        %         elseif nx==length(binsX) && ny~=length(binsY)
        %             Zbin(nx,ny)=Z(X>=binsX(nx) & (Y>=binsY(ny) & Y<binsY(n+1)));
        %         elseif nx~=length(binsX) && ny==length(binsY)
        %             Zbin(nx,ny)=Z((X>=binsX(nx) & X<binsX(nx+1)) & Y>=binsY(ny));
        %         else
        if ~isempty(find((X>=binsX(nx) & X<binsX(nx+1)) & (Y>=binsY(ny) & Y<binsY(ny+1))))
            Zbin(nx,ny)=nanmean(Z((X>=binsX(nx) & X<binsX(nx+1)) & (Y>=binsY(ny) & Y<binsY(ny+1))));
        end
        %         end
    end
end
figure; format_fig; hold on;
imagesc(binsX(1:end-1),binsY(1:end-1),Zbin);