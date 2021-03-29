%%
clear all
close all

table_all=readtable('wanderIM_LocalSleep_PPtiesPerProbe_per6sWin.csv');
table_all.SubID=categorical(table_all.SubID);
table_all.Task=categorical(table_all.Task);
table_all.State=categorical(table_all.State);

%%
recompute=0;
totperm=1000;

% if recompute==1
nWin=2;

MWvsON_nWa_est=cell(1,2);
MWvsON_P2P_est=cell(1,2);
MWvsON_DwS_est=cell(1,2);
MWvsON_UpS_est=cell(1,2);

MBvsON_nWa_est=cell(1,2);
MBvsON_P2P_est=cell(1,2);
MBvsON_DwS_est=cell(1,2);
MBvsON_UpS_est=cell(1,2);

MWvsMB_nWa_est=cell(1,2);
MWvsMB_P2P_est=cell(1,2);
MWvsMB_DwS_est=cell(1,2);
MWvsMB_UpS_est=cell(1,2);

sub_table_all=table_all(table_all.Bin==nWin,:);

tic;
fprintf('E:%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    %%%%% MW vs ON
    [real_out, perm_out]=lme_perm(sub_table_all(sub_table_all.Elec==nE,:),'nWave','MWvsON~1+pred+Task+(1|SubID)',totperm);
    MWvsON_nWa_est{1}=[MWvsON_nWa_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MWvsON_nWa_est{2}=[MWvsON_nWa_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
end
save(['HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsON_nWa_Win' num2str(nWin) '.mat'],...
    'MWvsON_nWa_est')
clear MWvsON_nWa_est
toc;
tic;
fprintf('E:%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    [real_out, perm_out]=lme_perm(sub_table_all(sub_table_all.Elec==nE,:),'P2P','MWvsON~1+pred+Task+(1|SubID)',totperm);
    MWvsON_P2P_est{1}=[MWvsON_P2P_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MWvsON_P2P_est{2}=[MWvsON_P2P_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
end
save(['HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsON_P2P_Win' num2str(nWin) '.mat'],...
    'MWvsON_P2P_est')
clear MWvsON_P2P_est
toc;
tic;
fprintf('E:%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    [real_out, perm_out]=lme_perm(sub_table_all(sub_table_all.Elec==nE,:),'DownSlope','MWvsON~1+pred+Task+(1|SubID)',totperm);
    MWvsON_DwS_est{1}=[MWvsON_DwS_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MWvsON_DwS_est{2}=[MWvsON_DwS_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
end
save(['HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsON_DwS_Win' num2str(nWin) '.mat'],...
    'MWvsON_DwS_est')
clear MWvsON_DwS_est
toc;
tic;
fprintf('E:%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    [real_out, perm_out]=lme_perm(sub_table_all(sub_table_all.Elec==nE,:),'UpSlope','MWvsON~1+pred+Task+(1|SubID)',totperm);
    MWvsON_UpS_est{1}=[MWvsON_UpS_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MWvsON_UpS_est{2}=[MWvsON_UpS_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
end
save(['HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsON_UpS_Win' num2str(nWin) '.mat'],...
    'MWvsON_UpS_est')
clear MWvsON_UpS_est
toc;
tic;
fprintf('E:%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    %%%%% MB vs ON
    [real_out, perm_out]=lme_perm(sub_table_all(sub_table_all.Elec==nE,:),'nWave','MBvsON~1+pred+Task+(1|SubID)',totperm);
    MBvsON_nWa_est{1}=[MBvsON_nWa_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MBvsON_nWa_est{2}=[MBvsON_nWa_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
end
save(['HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MBvsON_nWa_Win' num2str(nWin) '.mat'],...
    'MBvsON_nWa_est')
clear MBvsON_nWa_est
toc;
tic;
fprintf('E:%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    [real_out, perm_out]=lme_perm(sub_table_all(sub_table_all.Elec==nE,:),'P2P','MBvsON~1+pred+Task+(1|SubID)',totperm);
    MBvsON_P2P_est{1}=[MBvsON_P2P_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MBvsON_P2P_est{2}=[MBvsON_P2P_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
end
save(['HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MBvsON_P2P_Win' num2str(nWin) '.mat'],...
    'MBvsON_P2P_est')
clear MBvsON_P2P_est
toc;
tic;
fprintf('E:%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    [real_out, perm_out]=lme_perm(sub_table_all(sub_table_all.Elec==nE,:),'DownSlope','MBvsON~1+pred+Task+(1|SubID)',totperm);
    MBvsON_DwS_est{1}=[MBvsON_DwS_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MBvsON_DwS_est{2}=[MBvsON_DwS_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
end
save(['HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MBvsON_DwS_Win' num2str(nWin) '.mat'],...
    'MBvsON_DwS_est')
clear MBvsON_DwS_est
toc;
tic;
fprintf('E:%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    [real_out, perm_out]=lme_perm(sub_table_all(sub_table_all.Elec==nE,:),'UpSlope','MBvsON~1+pred+Task+(1|SubID)',totperm);
    MBvsON_UpS_est{1}=[MBvsON_UpS_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MBvsON_UpS_est{2}=[MBvsON_UpS_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
end
save(['HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MBvsON_UpS_Win' num2str(nWin) '.mat'],...
    'MBvsON_UpS_est')
clear MBvsON_UpS_est
toc;
tic;
fprintf('E:%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    %%%%% MW vs MB
    [real_out, perm_out]=lme_perm(sub_table_all(sub_table_all.Elec==nE,:),'nWave','MWvsMB~1+pred+Task+(1|SubID)',totperm);
    MWvsMB_nWa_est{1}=[MWvsMB_nWa_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MWvsMB_nWa_est{2}=[MWvsMB_nWa_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
end
save(['HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsMB_nWa_Win' num2str(nWin) '.mat'],...
    'MWvsMB_nWa_est')
clear MWvsMB_nWa_est
toc;
tic;
fprintf('E:%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    [real_out, perm_out]=lme_perm(sub_table_all(sub_table_all.Elec==nE,:),'P2P','MWvsMB~1+pred+Task+(1|SubID)',totperm);
    MWvsMB_P2P_est{1}=[MWvsMB_P2P_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MWvsMB_P2P_est{2}=[MWvsMB_P2P_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
end
save(['HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsMB_P2P_Win' num2str(nWin) '.mat'],...
    'MWvsMB_P2P_est')
clear MWvsMB_P2P_est
toc;
tic;
fprintf('E:%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    [real_out, perm_out]=lme_perm(sub_table_all(sub_table_all.Elec==nE,:),'DownSlope','MWvsMB~1+pred+Task+(1|SubID)',totperm);
    MWvsMB_DwS_est{1}=[MWvsMB_DwS_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MWvsMB_DwS_est{2}=[MWvsMB_DwS_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
end
save(['HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsMB_DwS_Win' num2str(nWin) '.mat'],...
    'MWvsMB_DwS_est')
clear MWvsMB_DwS_est
toc;
tic;
fprintf('E:%2.0f/63\n',0)
for nE=1:63
    fprintf('\b\b\b\b\b\b\b\bE:%2.0f/63\n',nE)
    [real_out, perm_out]=lme_perm(sub_table_all(sub_table_all.Elec==nE,:),'UpSlope','MWvsMB~1+pred+Task+(1|SubID)',totperm);
    MWvsMB_UpS_est{1}=[MWvsMB_UpS_est{1} ; [nE real_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
    MWvsMB_UpS_est{2}=[MWvsMB_UpS_est{2} ; [nE*ones(totperm,1) perm_out]];% double(rtGO_mdl.Coefficients(4,2)) double(rtGO_mdl.Coefficients(4,4)) double(rtGO_mdl.Coefficients(4,6))]];
end
save(['HDDM_WIM_thrE90P2P_P2P_Amp_State_LMEandPerm_rev1_MWvsMB_UpS_Win' num2str(nWin) '.mat'],...
    'MWvsMB_UpS_est')
clear MWvsMB_UpS_est
toc;

