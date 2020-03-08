function data_pupil=get_EyeLink_cleanpupil(rawdata_pupil,SR,data_time,EL_events)

data_pupil=rawdata_pupil;
blinks=EL_events.Blinks;
fprintf('... ... cleaning blinks (NaN replacement)\n')
fprintf('%3.0f%%\n',0)
for nbl=1:length(blinks.start)
    fprintf('\b\b\b\b\b%3.0f%%\n',nbl/length(blinks.start)*100)
    startB=find(data_time==blinks.start(nbl));
    endB=find(data_time==blinks.end(nbl));
    %     data_pupil(startB(1):endB(1))=nan;
    %     temp_interpolate=[mode(data_pupil(startB-0.2*SR:startB-0.1*SR)) ; mode(data_pupil(endB+0.1*SR:endB+0.2*SR))];
    if startB-0.1*SR<1 || endB+0.1*SR>length(data_pupil)
        data_pupil(max([startB(1)-0.1*SR 1]):min([endB(1)+0.1*SR length(data_pupil)]))=nan;
    else
        temp_interpolate=[(data_pupil(startB-0.1*SR)) ; (data_pupil(endB+0.1*SR))];
        %     [p,S,mu] =polyfit([0 1],temp_interpolate',1);
        p=[];
        p(2)=temp_interpolate(1);
        p(1)=diff(temp_interpolate);
        temp_interpolated=p(1)*(1/length(startB(1)-0.1*SR:endB(1)+0.1*SR):1/length(startB(1)-0.1*SR:endB(1)+0.1*SR):1)+p(2);
        data_pupil(startB(1)-0.1*SR:endB(1)+0.1*SR)=temp_interpolated;
    end
    %     figure;
    %     plot(rawdata_pupil(startB(1)-0.2*SR:endB(1)+0.2*SR),'Color','r');
    %     hold on
    %     plot(data_pupil(startB(1)-0.2*SR:endB(1)+0.2*SR),'Color','k');
    %     pause;
end
fprintf('... ... done\n')
data_pupil(data_pupil<prctile(data_pupil(~isnan(data_pupil)),0.1))=NaN;
