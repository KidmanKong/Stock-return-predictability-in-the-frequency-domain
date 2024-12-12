clear;
addpath('jplv7')
input_file='data.xls';
input_sheet='Equity premium';
y=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:b1153');
input_sheet='Macroeconomic variables';
predictor=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:o1153');
input_sheet='DP_DG_EG';
DP_DG_EG=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:d1153');
T=size(y,1);
N=size(predictor,2);
n_wd=2;
R=240;
P=T-R;

%% Factor predictor
% EW
F_EW=mean(zscore(predictor),2);
% PC
[~,F_PC]=pca(zscore(predictor));
F_PC=F_PC(:,1);
% PLS
F_PLS_temp=zeros(T,4);
F_PLS=zeros(T,4);
for t=2:T
    y_t=y(1:t);
    predictor_t=predictor(1:t,:);
    predictor_t(:,[1 2 4])=detrend(predictor_t(:,[1 2 4]),1);
    y_comp=wavelet_decomposing_function(y_t,'haar',n_wd);
    y_comp=[y_t y_comp(:,end:-1:1)];
    predictor_t_s=zscore(predictor_t);
    for i=1:size(y_comp,2)
        pai=nan(N,1);
        for n=1:N
            predictor_t_s(:,n)=winsor(predictor_t_s(:,n),[2 98]);
            x_t=predictor_t_s(:,n);
            beta=regress(x_t(1:end-1),[ones(length(x_t(1:end-1)),1) y_comp(2:end,i)]);
            pai(n)=beta(end);
        end
        beta=regress(predictor_t_s(t,:)',[ones(length(pai),1) pai]);
        F_PLS_temp(t,i)=beta(end);
        F_PLS(t,i)=beta(end);
    end
end

F_EW=zscore(F_EW(R+1:end));
F_PC=zscore(F_PC(R+1:end));
F_PLS=zscore(F_PLS(R+1:end,4));

DP=DP_DG_EG(R+1:end,1);
y=DP_DG_EG(R+1:end,:);

%% predictor regression
[T,M]=size(y);
results_EW=nan(M,7);
results_PC=nan(M,7);
results_PLS=nan(M,7);

for i=1:M
    % EW
    OLS=nwest(y(2:end,i),[ones(T-1,1) F_EW(1:end-1) DP(1:end-1)],12);
    results_EW(i,:)=[100*OLS.beta(2) OLS.tstat(2) 1-normcdf(abs(OLS.tstat(2)))...
                       100*OLS.beta(3) OLS.tstat(3) 1-normcdf(abs(OLS.tstat(3)))...;
                       100*OLS.rsqr];
    % PC
    OLS=nwest(y(2:end,i),[ones(T-1,1) F_PC(1:end-1) DP(1:end-1)],12);
    results_PC(i,:)=[100*OLS.beta(2) OLS.tstat(2) 1-normcdf(abs(OLS.tstat(2)))...
                       100*OLS.beta(3) OLS.tstat(3) 1-normcdf(abs(OLS.tstat(3)))...;
                       100*OLS.rsqr];
    % PLS
    OLS=nwest(y(2:end,i),[ones(T-1,1) F_PLS(1:end-1) DP(1:end-1)],12);
    results_PLS(i,:)=[100*OLS.beta(2) OLS.tstat(2) 1-normcdf(abs(OLS.tstat(2)))...
                       100*OLS.beta(3) OLS.tstat(3) 1-normcdf(abs(OLS.tstat(3)))...;
                       100*OLS.rsqr];
end