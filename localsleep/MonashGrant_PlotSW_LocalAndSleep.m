%% load SCP SWs
load('/Users/tand0009/Work/Documents/Grants/Monash/SPG_2019/SCP_SlowWaves.mat')


%% plot
timeP=-1.0:1/100:1.0;
figure; format_fig;
% simpleTplot(timeP,squeeze(allKC(countKC>16,:,65)),0,[1 1 1]*0.5,0,'-',0.5,1,4);
simpleTplot(timeP,squeeze(allSO(countSO>16,:,65)),0,'k',0,'-',0.5,1,4);
simpleTplot(-1:1/D.fsample:2,squeeze((ERP_Waves(:,:))),0,'r',[0],'-',0.5,1,2,[],2);
xlim([-0.5 1]);
xlabel('Time from wave onset (s)')
ylabel('Amplitude (\muV)')

export_fig('/Users/tand0009/Desktop/SW_ERP_SCP_WanderIM.eps')

%%
figure;
set(gcf,'position',[100   100   220   240])
timeP=-1.0:1/100:1.0;
simpleTplot(-1:1/D.fsample:2,squeeze((ERP_Waves(:,:))),0,'k',[0],'-',0.5,1,2,[],2);
xlim([-0.4 0.8]); format_fig
xlabel('Time (s)')
ylabel('Amplitude (\muV)')

figure;
set(gcf,'position',[100   100   220   240])
timeP=-1.0:1/100:1.0;
simpleTplot(timeP,squeeze(allSO(countSO>16,:,65)),0,'r',0,'-',0.5,1,2,[],2);
xlim([-0.4 0.8]); format_fig
xlabel('Time (s)')
ylabel('Amplitude (\muV)')