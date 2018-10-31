%%
mat_sim=[-1+randn(5000,1) ; randn(500,1) ; 1+randn(5000,1) ; randn(500,1)];
mat_sim2=[[ones(5000,1) ; zeros(500,1) ; ones(5000,1) ; zeros(500,1)] [ones(5500,1) ; zeros(5500,1)]];
tbl_sim=array2table([mat_sim mat_sim2],'VariableNames',{'var','p1','p2'});
tbl_sim.p1=categorical(tbl_sim.p1);
tbl_sim.p2=categorical(tbl_sim.p2);
mdl= fitlme(tbl_sim,'var~p1*p2');