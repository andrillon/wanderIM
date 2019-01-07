function [filtdata]=meanfilt(datatofilt,pts);
if length(datatofilt)>=pts
    filtdata=[];
    ptsaway=floor(pts/2);
    filtdata([1:pts])=datatofilt([1:pts]);
    filtdata([length(datatofilt)-(pts-1):length(datatofilt)])=datatofilt([length(datatofilt)-(pts-1):length(datatofilt)]);
    for wndw=pts-ptsaway:length(datatofilt)-(ptsaway)
    filtdata(wndw)=mean(datatofilt([wndw-(ptsaway):wndw+(ptsaway)]));
    end
else filtdata=datatofilt;
end
