    function [twa_results]=twalldetectnew_TA_v2(datainput,orig_fs,thramp)
%% Theta Wave Detection
% Script modified on 2016.03.31 for Visual Deprivation Project (CIRS)
%
% Half-waves with duration 0.1-1.0s are used. A negative amplitude
% threshold may be defined by the user (however, a zero threshold is
% advised for most analyses).

disp(['**** Delta/Theta Wave Detection ****']);

 %% Baseline Correction

allavg=nanmean(datainput,2);
dataref=datainput-repmat(allavg,[1,size(datainput,2)]);
LoadedEEG.data=dataref; clear dataref;  

%% Filter Definition

fs=128; %sampling rate changes for decimated signal
Wp=[1.0 10.0]/(fs/2); % Filtering parameters
Ws=[0.1 15]/(fs/2); % Filtering parameters
Rp=3;
Rs=25;
[n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
[bbp,abp]=cheby2(n,Rs,Wn); % Loses no more than 3 dB in pass band and has at least 10 dB attenuation in stop band
clear pass* stop* Rp Rs W* n;

%% Detection Starts Here

clear datapoints swa_results channels datax signal dataff EEGder EEG difference ;
clear pos_index neg_index troughs peaks poscross negcross wndx b bx c cx nump maxb maxc lastpk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([' Analyzing ', num2str(size(LoadedEEG.data,1)),' channels']);
% h = waitbar(0,'Detection Progress');
fprintf('\n');
for i=1:size(LoadedEEG.data,1) % (i) is the number of the channel
    fprintf('... channel %3.0f/%3.0f - %3.0f%%',i,size(LoadedEEG.data,1),0)
%     waitbar(i/size(LoadedEEG.data,1),h,'Detection Progress')
    % Data Extraction, Resample and Filtering
    datax = squeeze(LoadedEEG.data(i,1:size(LoadedEEG.data,2),1));
	signal = resample(double(datax),fs,orig_fs);
    EEG=filtfilt(bbp, abp, signal);
    datapoints=length(EEG);
    channels(i).datalength=datapoints;
    
    %% Search Negcross and Poscross
    pos_index=zeros(length(EEG),1); % Create an empty list for positive peaks
    pos_index(find(EEG>0))=1; % Index of all positive points for EEG
    difference=diff(pos_index); % Locate the first positive and negative point (in time) for each series of consecutive points
    poscross=find(difference==1); % Return the position of all first positive points
    negcross=find(difference==-1); % Return the position of all first negative points
    EEGder=meanfilt(diff(EEG),5); % Meanfilt is a function that uses a 5 sample moving window to smooth derivative
    pos_index=zeros(length(EEGder),1); % Repeat above procedure on smoothed signal (?)
    pos_index(find(EEGder>0))=1; % Index of all positive points above minimum threshold %%!!!!!!!!!!!!!!!!changed to 0!!!!!!!!!!!!!
    difference=diff(pos_index); % Locate first positive and negative points
    peaks=find(difference==-1)+1; % Find pos ZX and neg ZX of the derivative (peaks)
    troughs=find(difference==1)+1; % Find pos ZX and neg ZX of the derivative (troughs)
    peaks(EEG(peaks)<0)=[]; % Rejects peaks below zero 
    troughs(EEG(troughs)>0)=[]; % Rejects troughs above zero
    
    %% Makes negcross and poscross same size to start
    if isempty(negcross) || isempty(poscross)
        continue;
    end
    if negcross(1)<poscross(1); 
            start=1;
        else
            start=2;
    end; 
    
    if start==2;
            poscross(1)=[];
    end;
        lastpk=NaN; % Way to look at Peak to Peak parameters if needed
        ch=i;
                
%% Wave parameters initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	channels(ch).negzx{1}=[];
	channels(ch).poszx{1}=[];
	channels(ch).wvend{1}=[];
	channels(ch).negpks{1}=[];
	channels(ch).maxnegpk{1}=[];
	channels(ch).negpkamp{1}=[];
	channels(ch).maxnegpkamp{1}=[];
	channels(ch).pospks{1}=[];
	channels(ch).maxpospk{1}=[];
	channels(ch).pospkamp{1}=[];
	channels(ch).maxpospkamp{1}=[];
	channels(ch).mxdnslp{1}=[];
	channels(ch).mxupslp{1}=[];
	channels(ch).maxampwn{1}=[];
	channels(ch).minampwn{1}=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Locate Peaks
	wvi=1;
	for wndx=start:length(negcross)-1
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b... channel %3.0f/%3.0f - %3.0f%%',i,size(LoadedEEG.data,1),round(100*wndx/(length(negcross)-1)))
% 	disp([num2str(wndx)]);
    wavest=negcross(wndx); % Only used for neg/pos peaks
    wavend=negcross(wndx+1); % Only used for neg/pos peaks
        
    mxdn=abs(min(meanfilt(diff(EEG(wavest:poscross(wndx))),5)))*fs; % Matrix (27) determines instantaneous 1st segement slope on smoothed signal
    mxup=max(meanfilt(diff(EEG(wavest:poscross(wndx))),5))*fs; % Matrix (28) determines for 2nd segement
    negpeaks=troughs(troughs>wavest&troughs<wavend);
    wavepk=negpeaks(EEG(negpeaks)==min(EEG(negpeaks))); % Locate Wave Peak
    
    pospeaks=peaks(peaks>wavest&peaks<=wavend);

    if isempty(pospeaks)
        pospeaks=wavend; 
    end; %if pospeaks is empty set pospeak to pos ZX
                        
    %% period=wavend-wavest; %matrix(11) /fs
    poszx=poscross(wndx); %matrix(10)
    b=EEG(negpeaks);  
    bx=negpeaks;
    c=EEG(pospeaks);
    cx=pospeaks;
    nump=length(negpeaks); %matrix(24)

    %% get max and min amplitude over a window
    maxampwn=max(EEG(max(wavest-128,1):min(poszx+128,length(EEG))));
    minampwn=min(EEG(max(wavest-128,1):min(poszx+128,length(EEG))));
    %% Location of the max neg peak amp
    maxb=min(EEG(negpeaks)); % Max neg peak amp
    if maxb>0
        maxb=maxb(1);
    end;            
    maxbx=negpeaks(EEG(negpeaks)==maxb); % Location of the max neg peak amp

    %% Location of the max pos peak amp
    maxc=max(EEG(pospeaks)); % Max pos peak amp               
    if maxc>0
        maxc=maxc(1);
    end;            
    maxcx=pospeaks(EEG(pospeaks)==maxc); % Location of the max pos peak amp

    lastpk=maxcx;

	waveamp=abs(single(maxc))+abs(single(maxb)); %GBtest (peak2peak amplitude)
%     waveamp=abs(single(maxb)); %GBtest (neg peak amplitude) 
    wavelength=abs((single(wavest)-single(poszx))./fs); %GBtest (length)
	if wavelength>0.1 && wavelength<1.0 %GBtest (length)
        if waveamp>thramp %GBtest (p2p amplitude)
    
%% Wave parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    channels(ch).negzx{wvi}=round((single(wavest)./128).*orig_fs);
    channels(ch).poszx{wvi}=round((single(poszx)./128).*orig_fs);
    channels(ch).wvend{wvi}=round((single(wavend)./128).*orig_fs);
    channels(ch).negpks{wvi}={round((single(bx)./128).*orig_fs)};
    channels(ch).maxnegpk{wvi}=round((single(maxbx)./128).*orig_fs);
    channels(ch).negpkamp{wvi}=single(b);
    channels(ch).maxnegpkamp{wvi}=single(maxb);
    channels(ch).pospks{wvi}={round((single(cx)./128).*orig_fs)};
    channels(ch).maxpospk{wvi}=round((single(maxcx)./128).*orig_fs);
    channels(ch).pospkamp{wvi}=single(c);
    channels(ch).maxpospkamp{wvi}=single(maxc);
    channels(ch).mxdnslp{wvi}=single(mxdn);
    channels(ch).mxupslp{wvi}=single(mxup);
    channels(ch).maxampwn{wvi}=single(maxampwn);
    channels(ch).minampwn{wvi}=single(minampwn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     wvi=wvi+1;

        end; % (amplitude)
    end; % (duration)
    
    clear wavest wavend poszx bx maxbx b maxb cx maxcx x maxc mxdn mxup nump negpeaks pospeaks;
    clear wavelength waveamp;
    
     end; %end wndx loop
     clear EEGder dataff signal;
     fprintf('\n');
end; %end channel loop

%% Save Output

twa_results.channels=channels;
% close(h)

disp(['**** Detection Completed ****']);