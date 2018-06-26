tm=D.time;
ecgsig=D(match_str(D.chanlabels,'ECG'),:,1);

ecgsig = ft_preproc_bandstopfilter(ecgsig,D.fsample,[45 55],5,'but','twopass');
ecgsig=bandpass(ecgsig,D.fsample,0.05,40,4);
% ecgsig=(ecgsig-nanmedian(ecgsig))./(nanstd(ecgsig));

figure
plot(tm,ecgsig)
hold on
xlabel('Seconds')
ylabel('Amplitude')

qrsEx = ecgsig(661000+(-125:125)+40);
[mpdict,~,~,longs] = wmpdictionary(numel(qrsEx),'lstcpt',{{'sym4',3}});
figure
plot(qrsEx)
hold on
plot(2*circshift(mpdict(:,11),[-2 0]),'r')
axis tight
legend('QRS Complex','Sym4 Wavelet')
title('Comparison of Sym4 Wavelet and QRS Complex')

wt = modwt(ecgsig,5);
wtrec = zeros(size(wt));
wtrec(4:5,:) = wt(4:5,:);
y = imodwt(wtrec,'sym4');

y = abs(y).^2;
[qrspeaks,locs] = findpeaks(y,tm,'MinPeakHeight',5,...
    'MinPeakDistance',0.150);
figure
plot(tm,y)
hold on
plot(locs,qrspeaks,'ro')
xlabel('Seconds')
title('R Peaks Localized by Wavelet Transform with Automatic Annotations')

figure
plot(tm,ecgsig,'k--')
hold on
plot(tm,y,'r','linewidth',1.5)
plot(tm,abs(ecgsig).^2,'b')
set(gca,'xlim',[10.2 12])
legend('Raw Data','Wavelet Reconstruction','Raw Data Squared', ...
    'Location','SouthEast');
xlabel('Seconds')


[qrspeaks,locs] = findpeaks(ecgsig.^2,tm,'MinPeakHeight',0.35,...
    'MinPeakDistance',0.150);

figure
plot(tm,y)
title('R-Waves Localized by Wavelet Transform')
hold on
hwav = plot(locs,qrspeaks,'ro');
xlabel('Seconds')
legend([hwav hexp],'Automatic','Expert','Location','NorthEast');