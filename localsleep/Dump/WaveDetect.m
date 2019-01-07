function twa_results=WaveDetect(EEG,fs,epochl,seg);

close all;

% find zero crossing index (poscross & negcross) on the original EEG signal
pos_index=zeros(length(EEG),1);
pos_index(find(EEG>0))=1; %index of all positive points for EEG
difference=diff(pos_index);
poscross=find(difference==1)+1; % posZx index
negcross=find(difference==-1)+1; %negZX index

% Find Peaks and Troughs (using the transition of "EEGder"
% EEGder=meanfilt(diff(EEG),5); %meanfilt is a function that uses a 5 samples moving window to smooth derivative
EEGder=diff(EEG); 
pos_index=zeros(length(EEGder),1);
pos_index(find(EEGder>0))=1; %index of all positive points above minimum threshold
difference=diff(pos_index);
peaks=find(difference==-1)+1;
troughs=find(difference==1)+1; %find pos ZX and neg ZX of the derivative (the peaks & troughs)

% Double checking peaks and troughs
peaks(EEG(peaks)<0)=[]; % remove positive reflex which valuse < 0
troughs(EEG(troughs)>0)=[]; % remove negative reflex value of peaks > 0

% %%%%%%%%%%%%%%%%%% Plot EEG, Peaks and Troughs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(EEG,'b','LineWidth',2); hold on;
plot(EEGder,'r','LineWidth',2); hold on;
plot(difference,'g','LineWidth',2); hold on;
line([0,size(EEG,2)],[0 0],'LineWidth',0.5,'Color',[.8 .8 .8]); hold on;
ylim([-300 300]);
title('threshold:0');
hold on;
for i=1:size(peaks,1)
    line([peaks(i) peaks(i)],[-300 300],'LineWidth',0.5,'Color',[1 0.8 0]);
    hold on;
end;
for i=1:size(troughs,1)
    line([troughs(i) troughs(i)],[-300 300],'LineWidth',0.9,'Color',[0 0.6 0.5]);
    hold on;
end;
legend('EEG','EEGder','difference');

hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find out the end point of wave detection (always starts with negZx)

%     if negcross(1)<poscross(1);start=1;else start=2;end
%     if start==2;poscross(1)=[];end; % if the 1st negZx occurred later the 1st posZx, ignore the 1st posZx

Zcross=union(negcross, poscross);
Zcross(find(ismember(Zcross(:,1),poscross)==1),2)=1;
Zcross(find(ismember(Zcross(:,1),negcross)==1),2)=2;
Zcross_pos=find(Zcross(:,2)==1);
Zcross_neg=find(Zcross(:,2)==2);
if Zcross_neg(1)<Zcross_pos(1); start=1; else start=2; end;
if start==2;
    poscross(1)=[];
    Zcross=union(negcross, poscross);
    Zcross(find(ismember(Zcross(:,1),poscross)==1),2)=1;
    Zcross(find(ismember(Zcross(:,1),negcross)==1),2)=2;
    Zcross_pos=find(Zcross(:,2)==1);
    Zcross_neg=find(Zcross(:,2)==2);
end;
% if size(Zcross_neg,1)>size(Zcross_pos,1)
%     if Zcross(Zcross_pos(end)+1,2)-Zcross(Zcross_neg(end),2)~=0
%     negcross(end)=[];
%     Zcross_neg(end)=[];
% end;

wvi=1;
for wndx=1:length(Zcross_neg)-1
    if Zcross(Zcross_neg(wndx)+1,2)==1 & Zcross(Zcross_neg(wndx+1)-1,2)==1 & ...
            Zcross(Zcross_neg(wndx)+1,1)-Zcross(Zcross_neg(wndx),1)>1 & ...
            Zcross(Zcross_neg(wndx+1),1)-Zcross(Zcross_neg(wndx+1)-1,1)>1

        %% Wave Start and Wave End Index
        wavest=negcross(wndx); % matrix(7): wave start index: first negZx
        wavend=negcross(wndx+1)-1; %matrix(8): wave end index: following negZx
        wavest_amp=EEG(wavest); % wave start amplitude
        wavend_amp=EEG(wavend); % wave end amplitude

        plot([wavest:wavend],EEG(wavest:wavend),'k','LineWidth',3); hold on;
        waveidx=num2str(wndx);
        line([wavest wavend],[0 0],'LineWidth',0.5,'Color',[1 0.8 0]);
        text(wavest,-18,waveidx);

%         pause;

        %% Postive Zcross Location
        poszx=poscross(wndx); %matrix(10)
        poszx_amp=EEG(poszx);

        %% Location of All Negative Peaks
        negpeaks=troughs(troughs>wavest&troughs<wavend); % location of negpeaks

        %% Amplitudes of All Negative Peask
        negpeakamp=EEG(negpeaks);

        %% Location of the Max negative peak
        wavepk_neg=negpeaks(EEG(negpeaks)==min(EEG(negpeaks))); % matrix (13)

        %% Amplitude of the Max Negative Peak
        maxnegpeakamp=EEG(wavepk_neg); % matrix (12)

        %% Location of All Positive Peaks
        pospeaks=peaks(peaks>wavest&peaks<=wavend); % matrix (22)

        %% Amplitudes of All Positive Peaks
        pospeakamp=EEG(pospeaks); % matrix (21)

        %% Location of the Max Positive Peak
        wavepk_pos=pospeaks(EEG(pospeaks)==max(EEG(pospeaks)));

        %% Amplitudes of the Max Positive Peak
        maxpospeakamp=EEG(wavepk_pos);  % matrix (14)

        %% Wave Period
        period=length([wavest:wavend])/fs; % in sec % matrix(11)

        %% Negative 1/2 wave period (in sec)
        neghalfperiod=(Zcross(Zcross_neg(wndx)+1,1)-wavest)/fs; %matrix (25)

        %% Wave Mid-Point Location
        mdpt=wavest+ceil(neghalfperiod*fs/2); %matrix(9)
        
        %% Down and Up Slope in the Negative half wave (called Maximum Slope in Brady's paper, 2007)
        if wavest-5<0 && wavend-(poszx+5)<0
%             mxdn=(min(meanfilt(diff(EEG(wavest:poszx)),5)))/(1/fs); % original script
%             mxup=(max(meanfilt(diff(EEG(wavest:poszx)),5)))/(1/fs); % original script
            mxdn=(min(diff(EEG(wavest:poszx))))/(1/fs); % original script
            mxup=(max(diff(EEG(wavest:poszx))))/(1/fs); % original script
        elseif wavest-5<0 && wavend-(poszx+5)==0
%             x=meanfilt(diff(EEG(wavest:poszx+5)),5); %% adding samples for 5 sample moving average to prevent egde effects
            x=diff(EEG(wavest:poszx+5)); %% adding samples for 5 sample moving average to prevent egde effects
            mxdn=(min(x(1:end-5)))/(1/fs);
            mxup=(max(x(1:end-5)))/(1/fs);

        elseif wavest-5<0 && wavend-(poszx+5)>0
%             x=meanfilt(diff(EEG(wavest:poszx+5)),5); %% adding samples for 5 sample moving average to prevent egde effects
            x=diff(EEG(wavest:poszx+5)); %% adding samples for 5 sample moving average to prevent egde effects
            mxdn=(min(x(1:end-5)))/(1/fs);
            mxup=(max(x(1:end-5)))/(1/fs);
        
        elseif  wavest-5==0 && wavend-(poszx+5)<0
%             x=meanfilt(diff(EEG(wavest-4:poszx)),5); %% adding samples for 5 sample moving average to prevent egde effects
            x=diff(EEG(wavest-4:poszx)); %% adding samples for 5 sample moving average to prevent egde effects
            mxdn=(min(x(5:end)))/(1/fs);
            mxup=(max(x(5:end)))/(1/fs);

        elseif  wavest-5==0 && wavend-(poszx+5)==0
%             x=meanfilt(diff(EEG(wavest-4:poszx+5)),5); %% adding samples for 5 sample moving average to prevent egde effects
            x=diff(EEG(wavest-4:poszx+5)); %% adding samples for 5 sample moving average to prevent egde effects
            mxdn=(min(x(5:end-5)))/(1/fs);
            mxup=(max(x(5:end-5)))/(1/fs);
            
        elseif  wavest-5==0 && wavend-(poszx+5)>0
%             x=meanfilt(diff(EEG(wavest-4:poszx+5)),5); %% adding samples for 5 sample moving average to prevent egde effects
            x=diff(EEG(wavest-4:poszx+5)); %% adding samples for 5 sample moving average to prevent egde effects
            mxdn=(min(x(5:end-5)))/(1/fs);
            mxup=(max(x(5:end-5)))/(1/fs);

        elseif wavest-5>0 && wavend-(poszx+5)<0
%             x=meanfilt(diff(EEG(wavest-5:poszx)),5); %% adding samples for 5 sample moving average to prevent egde effects
            x=diff(EEG(wavest-5:poszx)); %% adding samples for 5 sample moving average to prevent egde effects
            mxdn=(min(x(6:end)))/(1/fs);
            mxup=(max(x(6:end)))/(1/fs);

        elseif wavest-5>0 && wavend-(poszx+5)==0
%             x=meanfilt(diff(EEG(wavest-5:poszx+5)),5); %% adding samples for 5 sample moving average to prevent egde effects
            x=diff(EEG(wavest-5:poszx+5)); %% adding samples for 5 sample moving average to prevent egde effects
            mxdn=(min(x(6:end-5)))/(1/fs);
            mxup=(max(x(6:end-5)))/(1/fs);

        else % wavest-5>0 && wavend-(poszx+5)>0
%             x=meanfilt(diff(EEG(wavest-5:poszx+5)),5); %% adding samples for 5 sample moving average to prevent egde effects
            x=diff(EEG(wavest-5:poszx+5)); %% adding samples for 5 sample moving average to prevent egde effects
            mxdn=(min(x(6:end-5)))/(1/fs);% matrix (27): instantaneous descending slope on the negative half wave
            mxup=(max(x(6:end-5)))/(1/fs);% matrix (28): instantaneous ascending slope on the negative half wave
            clear x;
        end;

        %% Max Negative Peak - to - Max Positive Peak "amplitude"
        maxp2p_amp=maxpospeakamp-maxnegpeakamp; % %matrix (16)

        %% Max Negative Peak - to - Max Positive Peak "period"
        maxp2p_duration=length([wavepk_neg:wavepk_pos])/fs; %matrix(26)

        %% Max (Negative Peak) - to - Max (Positive Peak) "slope"
        maxp2p_slope=(maxpospeakamp-maxnegpeakamp)/(length([wavepk_neg:wavepk_pos])/fs);

        %% Number of Negative Peaks (this is to determine if this is a multipeak wave: a wave with more than 1 negative peak between zero crossings)
        nump=length(negpeaks); %matrix(24)

        %% Average of Negative Peak Amplitudes
        meanAmp=abs(mean(EEG(negpeaks))); %matrix(23)

%         %% First Negative Peak Amplitude
%         n1=abs(EEG(negpeaks(1,1))); %matrix(17)
% 
%         %% First Negative Peak Location
%         n1x=negpeaks(1,1); %matrix(18)

        %% First-segment average slope (the amplitude of the most negative
        %% peak/time from the previous zero crossing)
        avgdnslp=(maxnegpeakamp-EEG(wavest))/(length(wavest:wavepk_neg)/fs);

        %% Second-segment average slope (the amplitude of the most negative peak/time from the previous zero crossing)
        avgupslp=(EEG(Zcross(Zcross_neg(wndx)+1,1))-maxnegpeakamp)/(length(wavepk_neg:Zcross(Zcross_neg(wndx)+1,1))/fs);

        %% Wave Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        waves(wvi,1).negzx=wavest; % Location of negative Zero crossing (wave start), Matrix(7)
        waves(wvi,1).poszx=poszx; % Location of positive Zero crossing, Matrix(10)
        waves(wvi,1).wvend=wavend; % Location of next negative zero crossing, Matrix(8)
        waves(wvi,1).negzxamp=wavest_amp; % Amplitude of negative Zero crossing (wave start)
        waves(wvi,1).poszxamp=poszx_amp; % Amplitude of positive Zero crossing
        waves(wvi,1).wvendamp=wavend_amp; % Amplitude of next negative zero crossing (wave end)
        waves(wvi,1).period=period; % whole wave period (=length(wavest:waveend),in seconds), Matrix(11)
        waves(wvi,1).neghalfwaveperiod=neghalfperiod; % Neg 1/2 wave period (in seconds), Matrix(25)
        waves(wvi,1).midpoint=mdpt; % Wave midpoint, Matrix(9)
        waves(wvi,1).numNegPeaks=nump; % Number of negative peaks, Matrix(24): this is used to determine multipeak waves
        waves(wvi,1).negpks=negpeaks; % all negative peak location, Matrix(18)
        waves(wvi,1).maxnegpk=wavepk_neg; % max negative peak location, Matrix(13)
        waves(wvi,1).negpkamp=negpeakamp; % all negative peak amplitude, Matrix(17)
        waves(wvi,1).maxnegpkamp=maxnegpeakamp; % max negative peak amplitude, Matrix(12)
        waves(wvi,1).pospks=pospeaks; % all positive peak location, Matrix(22)
        waves(wvi,1).maxpospk=wavepk_pos; % max positive peak location, Matrix(15)
        waves(wvi,1).pospkamp=pospeakamp; % all positive peak amplitude, Matrix(21)
        waves(wvi,1).maxpospkamp=maxpospeakamp; % max positive peak amplitude, Matrix(14)
        waves(wvi,1).maxp2pslp=maxp2p_slope; % slope from max negative peak to max positive peak
        waves(wvi,1).mxdnslp=mxdn; % first-segment maximum slope, Matrix(27)
        waves(wvi,1).mxupslp=mxup; % second-segment maximum slope, Matrix(28)
        waves(wvi,1).avgdnslp=avgdnslp; % first-segment average slope
        waves(wvi,1).avgupslp=avgupslp; % Second-segment average slope
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        wvi=wvi+1;
    end; % end of reject (non-wave) loop
end; % end of wndx loop
close all;


twa_results.waves=waves;