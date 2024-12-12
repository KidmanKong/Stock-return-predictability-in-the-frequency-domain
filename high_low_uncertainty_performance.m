clear;
addpath('jplv7')
input_file='data.xls';
input_sheet='Equity premium';
Y=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:b1153');
input_sheet='Macroeconomic variables';
predictor=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:o1153');
input_sheet='Uncertainty';
EU=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:b751');
N=size(predictor,2);
T=size(Y,1);
P=size(EU,1);
U=1*(EU>median(EU));
n_wd=2;

%% Factor predictor
% PLS
F_PLS_temp=zeros(T,4);
F_PLS=zeros(T,4);
for t=2:T
    y_t=Y(1:t);
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

M=size(F_PLS,2);
R2=nan(M,2);
state_reg_results=nan(M,7);
R2OS=nan(M,4);

%% in-sample
% high low uncertainty
F_predictor=F_PLS(end-P+1:end,:);
y=Y(end-P+1:end);
for i=1:M
    predictor=F_predictor(:,i);
    high_index=find(U(2:end));
    low_index=find(~U(2:end));
    NWEST=nwest(y(2:end),[ones(length(y(2:end)),1) zscore(predictor(1:end-1))],1);
    y_temp=y(2:end);
    y_high=y_temp(high_index);
    y_low=y_temp(low_index);
    y_hat_high=NWEST.yhat(high_index);
    y_hat_low=NWEST.yhat(low_index);
    R2(i,1)=100*(1-mean((y_high-y_hat_high).^2)/mean((y_high-mean(y_temp)).^2));
    R2(i,2)=100*(1-mean((y_low-y_hat_low).^2)/mean((y_low-mean(y_temp)).^2));
    NWEST=nwest(y(2:end),[ones(length(y(2:end)),1) U(1:end-1).*zscore(predictor(1:end-1)) (1-U(1:end-1)).*zscore(predictor(1:end-1))],1);
    state_reg_results(i,:)=[100*NWEST.beta(2) NWEST.tstat(2) 1-normcdf(abs(NWEST.tstat(2))) 100*NWEST.beta(3) NWEST.tstat(3) 1-normcdf(abs(NWEST.tstat(3))) 100*NWEST.rsqr];
end

%% out-of-sample
% high low uncertainty
R=120;
F_predictor=F_PLS(end-P-R+1:end,:);
y=Y(end-P-R+1:end);
actual=y(R+1:end);
for i=1:M
    FC_P=nan(P,1);
    FC_HA=nan(P,1);
    for t=1:P
        start=1;
        y_t=y(start:R+(t-1));
        FC_HA(t)=mean(y(1:R+(t-1)));
        predictor_t=F_predictor(start:R+(t-1),i);
        OLS=ols(y_t(2:end),[ones(length(y_t(2:end)),1) predictor_t(1:end-1)]);
        FC_P(t)=[1 predictor_t(end)]*OLS.beta;
    end
    
    high_index=find(U);
    low_index=find(~U);
    actual_high=actual(high_index);
    actual_low=actual(low_index);
    FC_HA_high=FC_HA(high_index);
    FC_HA_low=FC_HA(low_index);
    FC_P_high=FC_P(high_index);
    FC_P_low=FC_P(low_index);
    
    MSFE_HA_rec=mean((actual_high-FC_HA_high).^2);
    MSFE_i=mean((actual_high-FC_P_high).^2);
    R2OS_i=100*(1-(MSFE_i/MSFE_HA_rec));
    [~,p_value_i]=Perform_CW_test(actual_high,FC_HA_high,FC_P_high);
    R2OS(i,1:2)=[R2OS_i p_value_i];   
    
    MSFE_HA_low=mean((actual_low-FC_HA_low).^2);
    MSFE_i=mean((actual_low-FC_P_low).^2);
    R2OS_i=100*(1-(MSFE_i/MSFE_HA_low));
    [~,p_value_i]=Perform_CW_test(actual_low,FC_HA_low,FC_P_low);
    R2OS(i,3:4)=[R2OS_i p_value_i];
end