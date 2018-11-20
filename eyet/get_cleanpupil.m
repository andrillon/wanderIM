function data_pupil=get_cleanpupil(rawdata_pupil,data_time,EL_events)

data_pupil=rawdata_pupil;
blinks=EL_events.Blinks;
fprintf('... ... cleaning blinks (NaN replacement)\n')
fprintf('%3.0f%%\n',0)
for nbl=1:length(blinks.start)
    fprintf('\b\b\b\b\b%3.0f%%\n',nbl/length(blinks.start)*100)
    startB=find(data_time==blinks.start(nbl));
    endB=find(data_time==blinks.end(nbl));
    data_pupil(startB(1):endB(1))=nan;
    %     plot(rawdata_pupil(startB(1)-10:endB(1)+10));
    %     hold on
    %     plot(data_pupil(startB(1)-10:endB(1)+10));
    %     pause;
end
fprintf('... ... done\n')
