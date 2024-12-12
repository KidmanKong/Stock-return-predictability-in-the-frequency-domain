clear;
addpath('jplv7')
input_file='data.xls';
input_sheet='Equity premium';
y=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:b1153');
input_sheet='Macroeconomic variables';
predictor=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:o1153');
business_cycles=readmatrix(input_file,'Sheet',input_sheet,'Range','r2:r1153');
T=size(y,1);
N=size(predictor,2);
n_wd=2;
R=240;
P=T-R;
rec_index=find(business_cycles(R+1:R+P));
exp_index=find(~business_cycles(R+1:R+P));

%% Factor predictor
% EW
F_EW=mean(zscore(predictor),2);
% PC
[~,F_PC]=pca(zscore(predictor),'corr');
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
All_predictor=[predictor F_EW F_PC F_PLS];

%% predictor regression
N_All=size(All_predictor,2);
results_all=nan(N_All,6);
for i=1:N_All
    OLS=nwest(y(R+1:end),[ones(P,1) zscore(All_predictor(R:end-1,i))],12);
    results_all(i,1:2)=[100*OLS.beta(2) OLS.tstat(2)];
    results_all(i,4)=100*OLS.rsqr;
    y_temp=y(R+1:end);
    y_temp_exp=y_temp(exp_index);
    y_temp_rec=y_temp(rec_index);
    y_hat_exp=OLS.yhat(exp_index);
    y_hat_rec=OLS.yhat(rec_index);    
    results_all(i,5)=100*(1-mean((y_temp_exp-y_hat_exp).^2)/mean((y_temp_exp-mean(y_temp)).^2));
    results_all(i,6)=100*(1-mean((y_temp_rec-y_hat_rec).^2)/mean((y_temp_rec-mean(y_temp)).^2));
end

%% Computing wild bootstrapped p-values
load('Generate_wild_bootstrapped_pseudo_samples');
B=size(y_star,3); % number of boostrap replications
disp('Computing wild bootstrapped p-values');
tstat_All_star=nan(B,N_All);
for b=1:B
    y_b=y_star(:,1,b);
    predictor_b=predictor_star(:,:,b);
    F_EW_b=F_predictor_star(:,1,b);
    F_PC_b=F_predictor_star(:,2,b);
    F_PLS_b=F_predictor_star(:,3:end,b);
    All_predictor_b=[predictor_b F_EW_b F_PC_b F_PLS_b];
    for i=1:N_All
        OLS_b=nwest(y_b(R+1:end),[ones(P,1) All_predictor_b(R:end-1,i)],12);
        tstat_All_star(b,i)=OLS_b.tstat(2);
    end
    disp(b);
end
for i=1:N_All
    results_all(i,3)=mean(abs(tstat_All_star(:,i))>abs(results_all(i,2)));
    disp(i);
end
