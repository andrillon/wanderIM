f1=6;
f2=7.5;

flick1 = sin(2*pi*f1*(0:30*500)/500);
flick2 = sin(2*pi*f2*(0:30*500)/500);
flickIM=flick1.*flick2;

figure;
plot(0:1/500:30,flick1,'r');
hold on;
plot(0:1/500:30,flick2,'b');
plot(0:1/500:30,flickIM,'k');
xlim([0 2])
 
epoch1=[];
epoch2=[];
epochIM=[];
time=0:1/500:30;

for k=1:45
    t1=(k-1)*333+1;
    t2=(k)*333;
    epoch1=[epoch1 ; flick1(t1:t2)];
    epoch2=[epoch2 ; flick2(t1:t2)];
    epochIM=[epochIM ; flickIM(t1:t2)];
end

%%
baseline_1=[];
baseline_2=[];
baseline_3=[];
baseline_4=[];
OzIdx=17;
for k=1:3
    t1=((k-1)*4995+1)/D.fsample;
    t2=(k)*4995/D.fsample;
    temp=D_baseline(OzIdx,D_baseline.indsample(t1):D_baseline.indsample(t2),1);
    baseline_1=[baseline_1 ; temp];
    temp=D_baseline(OzIdx,D_baseline.indsample(t1):D_baseline.indsample(t2),5);
    baseline_1=[baseline_1 ; temp];
    
    temp=D_baseline(OzIdx,D_baseline.indsample(t1):D_baseline.indsample(t2),2);
    baseline_2=[baseline_2 ; temp];
    temp=D_baseline(OzIdx,D_baseline.indsample(t1):D_baseline.indsample(t2),6);
    baseline_2=[baseline_2 ; temp];
    
    
    temp=D_baseline(OzIdx,D_baseline.indsample(t1):D_baseline.indsample(t2),3);
    baseline_3=[baseline_3 ; temp];
    temp=D_baseline(OzIdx,D_baseline.indsample(t1):D_baseline.indsample(t2),7);
    baseline_3=[baseline_3 ; temp];
    
    
    temp=D_baseline(OzIdx,D_baseline.indsample(t1):D_baseline.indsample(t2),4);
    baseline_4=[baseline_4 ; temp];
    temp=D_baseline(OzIdx,D_baseline.indsample(t1):D_baseline.indsample(t2),8);
    baseline_4=[baseline_4 ; temp];
end

%%
figure; hold on;
[f,pow1]=get_PowerSpec(mean(baseline_1),D.fsample,0,0); 
plot(f,pow1,'r'); 
[f,pow2]=get_PowerSpec(mean(baseline_2),D.fsample,0,0); 
plot(f,pow2,'r--'); 
[f,pow3]=get_PowerSpec(mean(baseline_3),D.fsample,0,0); 
plot(f,pow3,'b'); 
[f,pow4]=get_PowerSpec(mean(baseline_4),D.fsample,0,0); 
plot(f,pow4,'b--'); 
xlim([1 20])

