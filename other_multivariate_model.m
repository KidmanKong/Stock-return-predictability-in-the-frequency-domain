clear;
addpath('jplv7')
input_file='data.xls';
input_sheet='Equity premium';
y=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:b1153');
input_sheet='Macroeconomic variables';
predictor=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:o1153');
Rfree_lag=readmatrix(input_file,'Sheet',input_sheet,'Range','q2:q1153');
T=size(y,1);
R=240;
P=T-R;
N=size(predictor,2);
n_wd=2;
wname='haar';

actual=y(R+1:end);
FC_HA=nan(P,1);
FC_All=nan(P,N);
FC_EKS=nan(P,1);
FC_LR=nan(P,1);
FC_RR=nan(P,1);
FC_PLS=nan(P,1);
FC_mean=nan(P,1);
FC_median=nan(P,1);
FC_trimmean=nan(P,1);

%% PLS-LC
for t=1:P
    start=1;
    y_t=y(start:R+(t-1));
    predictor_t=predictor(start:R+(t-1),:);
    % PLS
    predictor_t(:,[1 2 4])=detrend(predictor_t(:,[1 2 4]),1);
    predictor_t_s=zscore(predictor_t);
    y_comp=wavelet_decomposing_function(y_t,'haar',n_wd);
    y_comp=y_comp(:,1);

    F_PLS=nan(size(y_comp));
    for i=1:size(y_comp,2)
        pai=nan(N,1);
        for n=1:N
            predictor_t_s(:,n)=winsor(predictor_t_s(:,n),[2 98]);
            x_t=predictor_t_s(:,n);
            beta=regress(x_t(1:end-1),[ones(length(x_t(1:end-1)),1) y_comp(2:end,i)]);
            pai(n,:)=beta(2:end)';
        end
        for tt=1:R+(t-start)
            beta=regress(predictor_t_s(tt,:)',[ones(length(pai),1) pai]);
            F_PLS(tt,i)=beta(end);
        end
        OLS=ols(y_t(2:end),[ones(length(y_t(2:end)),1) F_PLS(1:end-1,i)]);
        FC_PLS(t)=[1 F_PLS(end,i)]*OLS.beta;
    end
    disp(string(t)+'/'+string(P)+' PLS-LC')
end

%% other methods
for t=1:P
    y_t=y(1:R+t-1);
    FC_HA(t)=mean(y_t);
    % All
    for i=1:N
        x_t=predictor(1:R+t-1,i);
        OLS=ols(y_t(2:end),[ones(length(y_t(2:end)),1) x_t(1:end-1)]);
        FC_All(t,i)=[1 x_t(end)]*OLS.beta;
    end
    x_t=predictor(1:R+t-1,:);
    % EKS
    alpha_list=0.1:0.1:0.9;
    coef_list=nan(size(predictor,2),length(alpha_list));
    coef0_list=nan(1,length(alpha_list));
    mse_list=nan(1,length(alpha_list));
    temp=[];
    for m=1:1
        for n=1:length(alpha_list)        
            [B,FitInfo] = lasso(x_t(1:end-1,:),y_t(2:end),'Alpha',alpha_list(n),'Lambda',0.01:0.01:1,'CV',5);
            idxLambdaMinMSE = FitInfo.IndexMinMSE;
            coef = B(:,idxLambdaMinMSE);
            coef0 = FitInfo.Intercept(idxLambdaMinMSE);
            coef_list(:,n) = coef;
            coef0_list(n) = coef0;
            mse_list(n) = FitInfo.MSE(idxLambdaMinMSE);
        end
        [~,index]=min(mse_list);
        coef=coef_list(:,index);
        coef0=coef0_list(index);
        temp(end+1)=x_t(end,:)*coef+coef0;
    end
    FC_EKS(t)=mean(temp);
    % LR
    temp=[];
    for m=1:1
        [B,FitInfo]=lasso(x_t(1:end-1,:),y_t(2:end),'Alpha',1,'Lambda',0.01:0.01:1,'CV',5);
        idxLambdaMinMSE = FitInfo.IndexMinMSE;
        coef = B(:,idxLambdaMinMSE);
        coef0 = FitInfo.Intercept(idxLambdaMinMSE);
        temp(end+1)=x_t(end,:)*coef+coef0;
    end
    FC_LR(t)=mean(temp);
    % RR
    k_list=0.1:0.1:1;
    b_list=[];
    mse_list=[];
    for n=1:length(k_list)
        b=ridge(y_t(2:end),x_t(1:end-1,:),k_list(n),0);
        b_list(:,end+1)=b;
        mse_list(end+1)=mean(([ones(length(y_t(2:end)),1) x_t(1:end-1,:)]*b-y_t(2:end)).^2);
    end
    [~,index]=min(mse_list);
    b=b_list(:,index);
    FC_RR(t)=[1 x_t(end,:)]*b;
        
    disp(string(t)+'/'+string(P)+' other methods')
end
FC_DMSPE1=DMSPE(FC_All,actual,1);
FC_DMSPE2=DMSPE(FC_All,actual,0.9);
for t=1:P
    FC_mean(t)=mean(FC_All(t,:));
    FC_median(t)=median(FC_All(t,:));
    FC_trimmean(t)=Trimmean(FC_All(t,:));
end
input_sheet='BMA_BMS';
FC_BMA=xlsread(input_file,input_sheet,'b2:b913');
FC_BMS=xlsread(input_file,input_sheet,'c2:c913');
input_sheet='DMA_DMS';
FC_DMA=xlsread(input_file,input_sheet,'b2:b913');
FC_DMS=xlsread(input_file,input_sheet,'c2:c913');

FC_MULT=[FC_LR FC_EKS FC_RR FC_BMA FC_BMS FC_DMA FC_DMS FC_mean FC_median FC_trimmean FC_DMSPE1 FC_DMSPE2 FC_PLS];

R2OS_MULT=nan(size(FC_MULT,2),3);
% MSFE criterion, historical average
e_HA=(actual-FC_HA).^2;
MSFE_HA=mean(e_HA);
% MSFE criterion, predictor
for i=1:size(FC_MULT,2)
     MSFE_predictor_i=mean((actual-FC_MULT(:,i)).^2);
     R2OS_predictor_i=100*(1-(MSFE_predictor_i/MSFE_HA));
     [MSFE_adjusted_predictor_i,p_value_predictor_i]=Perform_CW_test(actual,...
        FC_HA,FC_MULT(:,i));
    R2OS_MULT(i,:)=[R2OS_predictor_i MSFE_adjusted_predictor_i p_value_predictor_i];    
end

% Portfolio performance
VOL_window=60;
FC_VOL=nan(P,1);
for t=1:P
    % Volatility forecast
    FC_VOL(t)=mean(y(R+(t-1)-VOL_window+1:R+(t-1)).^2)-...
        (mean(y(R+(t-1)-VOL_window+1:R+(t-1))))^2;
end
portfolio_MULT=nan(size(FC_MULT,2),2);
c_bp=0;
r_f_lag_P=Rfree_lag(R+1:R+P);
gamma_MV=3;
% Computing average utility: historical average
[v_HA,SR_HA,xxx,xxx,TO_HA]=Perform_asset_allocation(actual,r_f_lag_P,FC_HA,FC_VOL,gamma_MV,c_bp);
portfolio_HA=[1200*v_HA SR_HA 100*TO_HA];
% Computing average utility gains
for i=1:size(FC_MULT,2)
    [v_i,SR_i,xxx,xxx,TO_i]=Perform_asset_allocation(actual,r_f_lag_P,FC_MULT(:,i),FC_VOL,gamma_MV,c_bp);
    portfolio_MULT(i,:)=[1200*(v_i-v_HA) SR_i];
end

results_all=[R2OS_MULT portfolio_MULT];