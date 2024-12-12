clear;
addpath('jplv7')
input_file='data.xls';
input_sheet='Equity premium';
y=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:b1153');
input_sheet='Macroeconomic variables';
predictor=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:o1153');
input_sheet='sent';
sent=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:g685');
T=size(y,1);
R=240;
N=size(predictor,2);
n_wd=2;

%% Factor predictor
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

%% in-sample
y=y(463:462+length(sent));
F_PLS=F_PLS(463:462+length(sent),end);
x=[sent F_PLS];
N_sent=size(sent,2);

results_U=nan(N_sent,4);
for i=1:N_sent
    OLS=nwest(y(2:end),[ones(length(y(2:end)),1) zscore(x(1:end-1,i))],12);
    results_U(i,:)=[100*OLS.beta(2) OLS.tstat(2) 1-normcdf(abs(OLS.tstat(2))) 100*OLS.rsqr];
end

results_B=nan(N_sent,7);
for i=1:N_sent
    OLS=nwest(y(2:end),[ones(length(y(2:end)),1) zscore(F_PLS(1:end-1)) zscore(sent(1:end-1,i))],12);
    results_B(i,1:3)=[100*OLS.beta(2) OLS.tstat(2) 1-normcdf(abs(OLS.tstat(2)))];
    results_B(i,4:6)=[100*OLS.beta(3) OLS.tstat(3) 1-normcdf(abs(OLS.tstat(3)))];
    results_B(i,7)=100*OLS.rsqr;
end