clear;
addpath('jplv7')
input_file='data.xls';
input_sheet='Equity premium';
y=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:b1153');
input_sheet='Macroeconomic variables';
predictor=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:o1153');
T=length(y);
R=240;
N=size(predictor,2);
n_wd=2;

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

%% Forecast
All_predictor_temp=[F_EW F_PC F_PLS];
y_temp=y;
Y=y;
H=[3 6 9 12]; % forecasting horizon
results_all=nan(size(All_predictor_temp,2),6,length(H));

for hi=1:length(H)
    h=H(hi);
    y=nan(length(y_temp)-h+1,1);
    for t=1:length(y_temp)-h+1
        y(t)=mean(y_temp(t:t+h-1));
    end
    All_predictor=All_predictor_temp(1:length(y_temp)-h+1,:);
    
    % in-sample
    for i=1:size(All_predictor,2)
        x_t=zscore(All_predictor(R+1:end,i));
        y_t=y(R+1:end);
        OLS=nwest(y_t(2:end),[ones(length(y_t(2:end)),1) x_t(1:end-1)],12);
        results_all(i,1:4,hi)=[100*OLS.beta(2) OLS.tstat(2) 1-normcdf(abs(OLS.tstat(2))) 100*OLS.rsqr];
    end
    
    % out-of-sample
    T=length(y);
    P=T-R;
    actual=y(R+1:end);
    FC_HA=nan(P,1);
    FC_all=nan(P,size(All_predictor,2));
    
    for t=1:P
        start=1;
        y_t=y(start:R+t-h);
        Y_t=Y(start:R+t-1);
        FC_HA(t)=mean(y_t);
        predictor_t=predictor(start:R+(t-1),:);
        % EW
        F_EW=mean(zscore(predictor_t),2);
        % PC
        [~,F_PC]=pca(zscore(predictor_t),'corr');
        F_PC=F_PC(:,1);
        % PLS
        predictor_t(:,[1 2 4])=detrend(predictor_t(:,[1 2 4]),1);
        y_comp=wavelet_decomposing_function(Y_t,'haar',n_wd);
        y_comp=[Y_t y_comp(:,end:-1:1)];
        predictor_t_s=zscore(predictor_t);
        F_PLS=nan(size(y_comp));
        for i=1:size(y_comp,2)
            pai=nan(N,1);
            for n=1:N
                predictor_t_s(:,n)=winsor(predictor_t_s(:,n),[2 98]);
                x_t=predictor_t_s(:,n);
                beta=regress(x_t(1:end-1),[ones(length(x_t(1:end-1)),1) y_comp(2:end,i)]);
                pai(n,:)=beta(2:end)';
            end
            for tt=start:R+(t-1)
                beta=regress(predictor_t_s(tt,:)',[ones(length(pai),1) pai]);
                F_PLS(tt,i)=beta(end);
            end
        end
        All_predictor_t=[F_EW F_PC F_PLS];
        for i=1:size(All_predictor_t,2)
            x_t=All_predictor_t(:,i);
            OLS=ols(y_t(2:end),[ones(length(y_t(2:end)),1) x_t(1:end-h)]);
            FC_all(t,i)=[1 x_t(end)]*OLS.beta;
        end
        disp([t h])
    end
    
    % R2OS performance
    % MSFE criterion, historical average
    MSFE_HA=mean((actual-FC_HA).^2);
    % MSFE criterion, predictor
    for i=1:size(FC_all,2)
         MSFE_i=mean((actual-FC_all(:,i)).^2);
         R2OS_i=100*(1-(MSFE_i/MSFE_HA));
         [MSFE_adjusted_i,p_value_i_CW]=Perform_CW_test(actual,FC_HA,FC_all(:,i));
         results_all(i,5:6,hi)=[R2OS_i p_value_i_CW];
    end
end