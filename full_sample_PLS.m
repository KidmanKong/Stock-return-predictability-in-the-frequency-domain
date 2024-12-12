clear;
addpath('jplv7')
input_file='data.xls';
input_sheet='Equity premium';
y=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:b1153');
input_sheet='Macroeconomic variables';
predictor=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:o1153');
business_cycles=readmatrix(input_file,'Sheet',input_sheet,'Range','r2:r1153');

T=size(y,1);
R=240;
P=T-R;
N=size(predictor,2);
n_wd=2;
actual=y(R+1:R+P);
FC_HA=nan(P,1);
FC_PLS=nan(P,4);
rec_index=find(business_cycles(R+1:R+P));
exp_index=find(~business_cycles(R+1:R+P));

%% full sample PLS
F_PLS=zeros(T,4);
y_comp=wavelet_decomposing_function(y,'haar',n_wd);
y_comp=[y y_comp(:,end:-1:1)];
predictor(:,[1 2 4])=detrend(predictor(:,[1 2 4]),1);
predictor_s=zscore(predictor);
for i=1:size(y_comp,2)
    pai=nan(N,1);
    for n=1:N
        predictor_s(:,n)=winsor(predictor_s(:,n),[2 98]);
        x=predictor_s(:,n);
        beta=regress(x(1:end-1),[ones(length(x(1:end-1)),1) y_comp(2:end,i)]);
        pai(n)=beta(end);
    end
    for t=1:T
        beta=regress(predictor_s(t,:)',[ones(length(pai),1) pai]);
        F_PLS(t,i)=beta(end);
    end
end

N_All=size(F_PLS,2);
results_in_sample=nan(N_All,6);

%% in-sample performance
for i=1:N_All
    OLS=nwest(y(R+1:end),[ones(P,1) zscore(F_PLS(R:end-1,i))],12);
    results_in_sample(i,1:4)=[100*OLS.beta(2) OLS.tstat(2) 1-normcdf(abs(OLS.tstat(2))) 100*OLS.rsqr];
    y_temp=y(R+1:end);
    y_temp_exp=y_temp(exp_index);
    y_temp_rec=y_temp(rec_index);
    y_hat_exp=OLS.yhat(exp_index);
    y_hat_rec=OLS.yhat(rec_index);    
    results_in_sample(i,5)=100*(1-mean((y_temp_exp-y_hat_exp).^2)/mean((y_temp_exp-mean(y_temp)).^2));
    results_in_sample(i,6)=100*(1-mean((y_temp_rec-y_hat_rec).^2)/mean((y_temp_rec-mean(y_temp)).^2));
end

%% out-of-sample forecast
for t=1:P
    start=1;
    y_t=y(start:R+(t-1));
    predictor_t=F_PLS(start:R+(t-1),:);
    HA=mean(y(1:R+(t-1)));
    FC_HA(t)=HA;
    for i=1:size(F_PLS,2)
        results_predictor_i_t=ols(y_t(2:end),[ones(length(y_t(2:end)),1) predictor_t(1:end-1,i)]);
        FC_PLS(t,i)=[1 predictor_t(end,i)]*results_predictor_i_t.beta;
    end
end

%% R2OS performance
FC_all=FC_PLS;
R2OS=nan(size(FC_all,2),2);
% MSFE criterion, historical average
e_HA=(actual-FC_HA).^2;
MSFE_HA=mean(e_HA);
% MSFE criterion, predictor
for i=1:size(FC_all,2)
     MSFE_i=mean((actual-FC_all(:,i)).^2);
     R2OS_i=100*(1-(MSFE_i/MSFE_HA));
     [MSFE_adjusted_i,p_value_i_CW]=Perform_CW_test(actual,FC_HA,FC_all(:,i));
    R2OS(i,:)=[R2OS_i p_value_i_CW];   
end

actual_rec=actual(rec_index);
actual_exp=actual(exp_index);
FC_HA_rec=FC_HA(rec_index);
FC_HA_exp=FC_HA(exp_index);

FC_all_rec=FC_all(rec_index,:);
R2OS_rec=nan(size(FC_all_rec,2),2);
% MSFE criterion, historical average
MSFE_HA_res=mean((actual_rec-FC_HA_rec).^2);
% MSFE criterion, predictor
for i=1:size(FC_all_rec,2)
     MSFE_i=mean((actual_rec-FC_all_rec(:,i)).^2);
     R2OS_i=100*(1-(MSFE_i/MSFE_HA_res));
     [MSFE_adjusted_i,p_value_i]=Perform_CW_test(actual_rec,FC_HA_rec,FC_all_rec(:,i));
     R2OS_rec(i,:)=[R2OS_i p_value_i];   
end
FC_all_exp=FC_all(exp_index,:);
R2OS_exp=nan(size(FC_all_exp,2),2);
% MSFE criterion, historical average
MSFE_HA_exp=mean((actual_exp-FC_HA_exp).^2);
% MSFE criterion, predictor
for i=1:size(FC_all_exp,2)
     MSFE_i=mean((actual_exp-FC_all_exp(:,i)).^2);
     R2OS_i=100*(1-(MSFE_i/MSFE_HA_exp));
     [MSFE_adjusted_i,p_value_i]=Perform_CW_test(actual_exp,FC_HA_exp,FC_all_exp(:,i));
     R2OS_exp(i,:)=[R2OS_i p_value_i];   
end

%%
results_all=[results_in_sample R2OS R2OS_exp R2OS_rec];