%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
% close all;
run ../localdef_wanderIM

% adding relevant toolboxes to the path
% spm12 and LSCPtools
addpath(genpath(spm12_path))
addpath(genpath(lscpTools_path))
addpath(path_localsleep)

% select relevant files, here baseline blocks
eeg_path=[root_path filesep 'preproc_eeg'];
eeg_path2=[root_path filesep 'preproc_ica'];
behav_path=[root_path filesep 'behav'];
bsl_files=dir([eeg_path filesep 'nfEEG_S3*.mat']);
bsl_files(5,:)=[];

%% loop across trials for baseline blocks
prticle_Thr=90; % 80 or 90 or 95
LimFrqW=[1 4]; % [1 4] or [4 10]
AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
% fixThr=75;
fixThr=[];
art_ampl=200;
max_posampl=75;
max_Freq=7;
frontalElecs=[1 32 33 60];
for n=1:length(bsl_files)
    % load file with spm
    filename=bsl_files(n).name;
    
    % load behavioural results
    SubID=filename;
    SubID=SubID(findstr(SubID,'_S3')+2:findstr(SubID,'_S3')+4);
    if exist([eeg_path filesep 'wanderIM_cont_twa2_' SubID '.mat'])~=0 && exist([eeg_path2 filesep 'wanderIM_twa4_' SubID '.mat'])~=0
        fprintf('... load local sleep detection for subject %s\n',SubID)
        
        load([eeg_path filesep 'wanderIM_twa3_' SubID]);
%         load([eeg_path2 filesep 'wanderIM_twa4_' SubID]);
    else
        fprintf('... load local sleep detection subject %s DOES NOT EXIST\n',SubID)
        continue;
    end
    if AmpCriterionIdx==9
        all_Waves(:,AmpCriterionIdx)=-all_Waves(:,AmpCriterionIdx);
    end
    hdr=ft_read_header([eeg_path filesep filename]);
%     find_possibleEMs=find((all_Waves(:,4))>art_ampl);
%     leave=zeros(1,hdr.nSamples);
%     leaveidx=all_Waves(find_possibleEMs,5)+repmat((-0.5*hdr.Fs:0.5*hdr.Fs),length(find_possibleEMs),1);
%     leaveidx=reshape(leaveidx,1,numel(leaveidx)); leaveidx(leaveidx<1 | leaveidx>length(leave))=[]; leaveidx=unique(leaveidx);
%     leave(leaveidx)=1;
%     all_Wtodiscard=find(leave(all_Waves(:,5))==1);
% %     for nww=1:length(find_possibleEMs)
% %         all_Wtodiscard=[all_Wtodiscard ; find(all_Waves(:,5)>(all_Waves(find_possibleEMs(nww),5)-1*hdr.Fs) & all_Waves(:,5)<(all_Waves(find_possibleEMs(nww),5)+1*hdr.Fs))];
% %     end; 
% %     all_Wtodiscard=unique(all_Wtodiscard);
%     fprintf('... ... %g %% waves discarded because of vicinity to high amplitude events\n',length(all_Wtodiscard)/size(all_Waves,1)*100)
    %     all_Waves(all_Waves(:,AmpCriterionIdx)>max_amplW,:)=[];
%     all_Waves(all_Wtodiscard,:)=[];
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./hdr.Fs);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,2)>art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>max_posampl)*100)
    all_Waves(all_freq>max_Freq | all_Waves(:,2)>art_ampl | all_Waves(:,11)>max_posampl,:)=[];
    
%             save([eeg_path filesep 'wanderIM_cont_twa2_' SubID '_cleaned'],'all_Waves');

%     count=0; total=size(all_Waves,1);
%     for nE=1:63
%         tempW=find(all_Waves(:,3)==nE);
%         count=count+length(tempW(find(diff(all_Waves(tempW,5))/hdr.Fs<0.5)+1));
%         all_Waves(tempW(find(diff(all_Waves(tempW,5))/hdr.Fs<0.5)+1),:)=[];
%     end
%     fprintf('... ... %g %% waves discarded because of vicinity\n',count/total*100)
    all_Wpos=[];
    slow_Waves=[];
    vicinity_matrix2(n,:)=zeros(1,63);
    for nE=1:63
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_len=abs((thisE_Waves(:,5)-thisE_Waves(:,7)))./hdr.Fs;
        temp_freq=1./temp_len;
        temp_abs=1./temp_len;
        temp_p2p=thisE_Waves(:,AmpCriterionIdx);
        temp_p2p2=thisE_Waves(:,4);
        
        %         all_len=[all_len ; [repmat([n nE],length(temp_len),1) temp_len temp_abs]];
        %         thr_Wave1(n,nE)=prctile(all_Waves(all_Waves(:,3)==nE,4),80);
        %         thr_Wave2(n,nE)=prctile(temp_p2p(temp_freq>5),80);
        if ~isempty(fixThr)
            thr_Wave(n,nE)=fixThr;
        else
            %             thr_Wave2(n,nE)=prctile(all_Waves((temp_freq>=LimFrqW(1) & temp_freq<=LimFrqW(2)),AmpCriterionIdx),prticle_Thr);
            thr_Wave(n,nE)=prctile(all_Waves(:,AmpCriterionIdx),prticle_Thr);
        end
        dens_Waves(n,nE)=sum(temp_p2p>thr_Wave(n,nE))/(hdr.nSamples/hdr.Fs);
%         dens_Waves(n,nE)=sum(temp_p2p>thr_Wave(n,nE) & temp_p2p2>75)/(hdr.nSamples/hdr.Fs);
        %           len_Waves(n,nE)=sum(all_Waves(all_Waves(:,3)==nE,4)>thr_Wave2(n,nE));
        upSlope_Waves(n,nE)=mean(thisE_Waves(temp_p2p>thr_Wave(n,nE),12));
        downSlope_Waves(n,nE)=mean(thisE_Waves(temp_p2p>thr_Wave(n,nE),11));
        
        % concomittant
        chanA_waves=thisE_Waves(temp_p2p>thr_Wave(n,nE),:);
        for nE2=1:63
            if nE2==nE || isempty(chanA_waves)
                vicinity_matrix(n,nE,nE2)=nan;
                continue;
            end
            % chanel B
            thisE2_Waves=all_Waves(all_Waves(:,3)==nE2,:);
            temp2_len=abs((thisE2_Waves(:,5)-thisE2_Waves(:,7)))./hdr.Fs;
            temp2_freq=1./temp2_len;
            temp2_abs=1./temp2_len;
            temp2_p2p=thisE2_Waves(:,AmpCriterionIdx);
            chanB_waves=thisE2_Waves(temp2_p2p>thr_Wave(n,nE),:);
            all_posW_B=chanB_waves(:,5);
            neigh=nan(1,size(chanA_waves,1));
            neigh2=nan(1,size(chanA_waves,1));
            for nW=1:size(chanA_waves,1)
                posW=chanA_waves(nW,5);
                neigh(nW)=~isempty(find(all_posW_B>posW-0.5*hdr.Fs & all_posW_B<posW+0.5*hdr.Fs & chanB_waves(:,2)==chanA_waves(nW,2)));
            end
            vicinity_matrix(n,nE,nE2)=nanmean(neigh);
        end
        % chanel ALL
        thisE2_Waves=all_Waves(all_Waves(:,3)~=nE,:);
        temp2_len=abs((thisE2_Waves(:,5)-thisE2_Waves(:,7)))./hdr.Fs;
        temp2_freq=1./temp2_len;
        temp2_abs=1./temp2_len;
        temp2_p2p=thisE2_Waves(:,AmpCriterionIdx);
        chanB_waves=thisE2_Waves(temp2_p2p>thr_Wave(n,nE),:);
        all_posW_B=chanB_waves(:,5);
        neigh2=nan(1,size(chanA_waves,1));
        for nW=1:size(chanA_waves,1)
            posW=chanA_waves(nW,5);
            neigh2(nW)=length(unique(chanB_waves(find(all_posW_B>posW-0.5*hdr.Fs & all_posW_B<posW+0.5*hdr.Fs),3)));
        end
        vicinity_matrix2(n,nE)=nanmean(neigh2);
            
        temp=chanA_waves(:,5);
        temp(find(diff(temp)/hdr.Fs<0.5)+1)=[];
%         all_Wpos=[all_Wpos ; [temp nE*ones(length(temp),1)]];

        slow_Waves=[slow_Waves ; [chanA_waves neigh2']];

    end
    
%     save([eeg_path filesep 'wanderIM_cont_twa2_' SubID '_selected'],'slow_Waves');

    %     [~,idx] = sort(all_Wpos(:,1));
    %     all_Wpos=all_Wpos(idx,:);
    %     diffPos=[NaN ; diff(all_Wpos(:,1)/hdr.Fs)];
    %     breakSeq=find(diffPos>0.5);
    %     begS=1;
    %     endS=breakSeq(1)-1;
    %     begE=all_Wpos(1,2);
    %     endE=all_Wpos(endS,2);
    %     begP=all_Wpos(1,1);
    %     endP=all_Wpos(endS,1);
    %     temp_traj=[begS endS begE endE begP endP (endP-begP)/hdr.Fs length(unique(all_Wpos(1:endS,2)))];
    %     for nSeq=2:length(breakSeq)
    %         begS=breakSeq(nSeq-1);
    %         endS=breakSeq(nSeq)-1;
    %         begE=all_Wpos(begS,2);
    %         endE=all_Wpos(endS,2);
    %         begP=all_Wpos(begS,1);
    %         endP=all_Wpos(endS,1);
    %         temp_traj=[temp_traj ; [begS endS begE endE begP endP (endP-begP)/hdr.Fs length(unique(all_Wpos(begS:endS,2)))]];
    %     end
    %     temp_traj(temp_traj(:,end-1)>1,:)=[];
end

%%
dens_Waves(:,[10 21])=NaN; %Mastoids
thr_Wave(:,[10 21])=NaN; %Mastoids
vicinity_matrix(:,[10 21],:)=NaN;
vicinity_matrix(:,:,[10 21])=NaN;
vicinity_matrix2(:,[10 21])=NaN;
%%
% filename = '/Users/tand0009/Data/WanderIM/preproc_eeg/wanderIM_table_localsleep_results3.txt';
figure;
% subplot(1,3,1)
addpath(genpath(path_eeglab));
temp_topo=mean(dens_Waves,1)*60;
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','labels');
rmpath(genpath(path_eeglab));
colormap('parula'); colorbar; caxis([0 max(temp_topo)])
title('SW density /min')

%%
figure;
% subplot(1,3,1)
addpath(genpath(path_eeglab));
rmpath(genpath([path_eeglab filesep 'plugins/Biosig3.3.0']))
temp_topo=mean(thr_Wave,1);
topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off');
rmpath(genpath(path_eeglab));
colormap('parula'); colorbar; caxis([0 50])
title('SW threshold uV')


%%
figure;
% subplot(1,3,1)
addpath(genpath(path_eeglab));
temp_topo=squeeze(nanmean(vicinity_matrix2,1));
topoplot(temp_topo/64, layout.chaninfo,'style','map','whitebk','on','electrodes','off');
rmpath(genpath(path_eeglab));
colormap('parula'); colorbar; %caxis([0 0.4])
title('Spatial Expense Index')
%%
figure;
myLabels={'Fz','Cz','Pz','Oz'};
allLabels={layout.chaninfo.labels};
[X, Y]=pol2cart(pi/180*[layout.chaninfo.theta],[layout.chaninfo.radius]);
for nE=1:4
    myE=find(ismember(allLabels,myLabels{nE}));
    subplot(1,4,nE)
    addpath(genpath(path_eeglab));
    temp_topo=squeeze(nanmean(vicinity_matrix(:,myE,:),1));
%     temp_topo(myE)=1;
%     temp_topo(isnan(temp_topo))=0;
    topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off'); %,'emarker2',{myE});
%     topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off','emarker2',{myE-sum(find(isnan(temp_topo))<myE)});
%     scatter(X(myE),Y(myE),'MarkerFaceColor','r','Marker','o');
    rmpath(genpath(path_eeglab));
    colormap('parula'); 
%     colorbar; 
    caxis([0 0.55]);
    title(myLabels{nE}); format_fig
end

%%
figure;
[X, Y]=pol2cart(pi/180*[layout.chaninfo.theta],[layout.chaninfo.radius]);
for nE=[1 4]
    myE=find(ismember(allLabels,myLabels{nE}));
    posE=[X(myE) Y(myE)];
    for nE2=1:63
        posE2=[X(nE2) Y(nE2)];
        distE(nE,nE2)=sqrt((posE2(1)-posE(1)).^2+(posE2(2)-posE(2)).^2);
    end
    [~,sorteddis]=sort(distE(nE,:));
    hold on;
    tempplot=squeeze(nanmean(vicinity_matrix(:,myE,sorteddis),1));
    tempplot(isnan(tempplot))=[];
    plot(1:length(tempplot),tempplot,'LineWidth',2);
end
legend(myLabels([1 4]));
format_fig
xlabel('Dist to seed electrode')
ylabel('Proba of concomittant SW')
% % % %%
% % % MatCC=squeeze(nanmean(vicinity_matrix(:,:,:),1));
% % % MatCC(MatCC<0.6)=NaN;
% % % myLabels={layout.chaninfo.labels};
% % % simpleCCplot2(MatCC,myLabels(1:63),layout,1);
% % % caxis([0 1])
% % % %%
% % % A=1;
% % %
% % %

% %%
% figure;
% ChanLabels={layout.chaninfo.labels}; lay=layout;newP=1;
% ChanLabels=ChanLabels(1:63);
% Coor2D=layout.pos(match_str(ChanLabels,layout.label),:)';
% scatter(Coor2D(1,:), Coor2D(2,:),'ko');
% for i=1:length(layout.outline)
%     if ~isempty(layout.outline{i})
%         X = layout.outline{i}(:,1);
%         Y = layout.outline{i}(:,2);
%         line(X, Y, 'color', 'k', 'linewidth', 1.5, 'linestyle', '-');
%     end
% end
% allx=[];
% ally=[];
% for nt=1:size(temp_traj,1)
%     xbeg = Coor2D(1,temp_traj(nt,3));
%     ybeg = Coor2D(2,temp_traj(nt,3));
%     xend = Coor2D(1,temp_traj(nt,4));
%     yend = Coor2D(2,temp_traj(nt,4));
%
%     x = [xbeg xend]';
%     y = [ybeg yend]';
%     allx=[allx x];
%     ally=[ally y];
%     line(x, y,'Color', 'k');
% end
% line(mean(allx,2), mean(ally,2),'Color', 'r','LineWidth',4);
%
% %%
% for nE=1:63
%     pprtionSeg(nE,1)=mean(temp_traj(:,3)==nE);
%     pprtionSeg(nE,2)=mean(temp_traj(:,4)==nE);
% end
% figure;
% subplot(1,2,1)
% addpath(genpath(path_eeglab));
% temp_topo=pprtionSeg(:,1);
% topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off');
% rmpath(genpath(path_eeglab));
% colormap('parula'); colorbar; caxis([0 1]*0.05);
%
% subplot(1,2,2)
% addpath(genpath(path_eeglab));
% temp_topo=pprtionSeg(:,2);
% topoplot(temp_topo, layout.chaninfo,'style','map','whitebk','on','electrodes','off');
% rmpath(genpath(path_eeglab));
% colormap('parula'); colorbar; caxis([0 1]*0.05);