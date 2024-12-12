clear;
addpath('jplv7')
input_file='data.xls';

%% quarterly
input_sheet='Quarterly';
y=readmatrix(input_file,'Sheet',input_sheet,'Range','q2:q385');
predictor=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:o385');
T=size(y,1);
N=size(predictor,2);
n_wd=2;
R=80;
P=T-R;

actual=y(R+1:R+P);
FC_HA=nan(P,1);
FC_original=nan(P,N);
FC_EW=nan(P,1);
FC_PC=nan(P,1);
FC_PLS=nan(P,4);

% Factor predictor
% EW
F_EW=mean(zscore(predictor),2);
% PC
[~,F_PC]=pca(zscore(predictor),'corr');
F_PC=F_PC(:,1);
% PLS
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
            predictor_t_s(:,n)=winsor(predictor_t_s(:,n),[2,98]);
            x_t=predictor_t_s(:,n);
            beta=regress(x_t(1:end-1),[ones(length(x_t(1:end-1)),1) y_comp(2:end,i)]);
            pai(n)=beta(end);
        end
        beta=regress(predictor_t_s(t,:)',[ones(length(pai),1) pai]);
        F_PLS(t,i)=beta(end);
    end
end
All_predictor=[predictor F_EW F_PC F_PLS];

% in-sample performance
N_All=size(All_predictor,2);
results_in_sample=nan(N_All,4);
for i=1:N_All
    OLS=nwest(y(R+1:end),[ones(P,1) zscore(All_predictor(R:end-1,i))],12);
    results_in_sample(i,1:4)=[100*OLS.beta(2) OLS.tstat(2) 1-normcdf(abs(OLS.tstat(2))) 100*OLS.rsqr];
end

% out-of-sample forecast
for t=1:P
    start=1;
    y_t=y(start:R+(t-1));
    predictor_t=predictor(start:R+(t-1),:);
    HA=mean(y(1:R+(t-1)));
    FC_HA(t)=HA;
    % original
    for i=1:N
        results_predictor_i_t=ols(y_t(2:end),[ones(length(y_t(2:end)),1) predictor_t(1:end-1,i)]);
        FC_original(t,i)=[1 predictor_t(end,i)]*results_predictor_i_t.beta;
    end
    % EW
    F_EW=mean(zscore(predictor_t),2);
    OLS=ols(y_t(2:end),[ones(length(y_t(2:end)),1) F_EW(1:end-1)]);
    FC_EW(t)=[1 F_EW(end)]*OLS.beta;
    % PC
    [~,F_PC]=pca(zscore(predictor_t),'corr');
    F_PC=F_PC(:,1);
    OLS=ols(y_t(2:end),[ones(length(y_t(2:end)),1) F_PC(1:end-1)]);
    FC_PC(t)=[1 F_PC(end)]*OLS.beta;
    % PLS
    predictor_t(:,[1 2 4])=detrend(predictor_t(:,[1 2 4]),1);
    y_comp=wavelet_decomposing_function(y_t,'haar',n_wd);
    y_comp=[y_t y_comp(:,end:-1:1)];
    predictor_t_s=zscore(predictor_t);
    for i=1:size(y_comp,2)
        F_PLS=nan(size(y_comp));
        pai=nan(N,1);
        for n=1:N
            predictor_t_s(:,n)=winsor(predictor_t_s(:,n),[2,98]);
            x_t=predictor_t_s(:,n);
            beta=regress(x_t(1:end-1),[ones(length(x_t(1:end-1)),1) y_comp(2:end,i)]);
            pai(n,:)=beta(2:end)';
        end
        for tt=1:R+(t-start)
            beta=regress(predictor_t_s(tt,:)',[ones(length(pai),1) pai]);
            F_PLS(tt,i)=beta(end);
        end
        OLS=ols(y_t(2:end),[ones(length(y_t(2:end)),1) F_PLS(1:end-1,i)]);
        FC_PLS(t,i)=[1 F_PLS(end,i)]*OLS.beta;
    end
    disp(t)
end

FC_F=[FC_EW FC_PC FC_PLS];

% R2OS performance
FC_all=[FC_original FC_F];
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

results_all_Q=[results_in_sample R2OS];

%% annual
input_sheet='annual';
y=readmatrix(input_file,'Sheet',input_sheet,'Range','q2:q97');
predictor=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:o97');
T=size(y,1);
N=size(predictor,2);
n_wd=2;
R=40;
P=T-R;

actual=y(R+1:R+P);
FC_HA=nan(P,1);
FC_original=nan(P,N);
FC_EW=nan(P,1);
FC_PC=nan(P,1);
FC_PLS=nan(P,4);

% Factor predictor
% EW
F_EW=mean(zscore(predictor),2);
% PC
[~,F_PC]=pca(zscore(predictor),'corr');
F_PC=F_PC(:,1);
% PLS
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
            predictor_t_s(:,n)=winsor(predictor_t_s(:,n),[2,98]);
            x_t=predictor_t_s(:,n);
            beta=regress(x_t(1:end-1),[ones(length(x_t(1:end-1)),1) y_comp(2:end,i)]);
            pai(n)=beta(end);
        end
        beta=regress(predictor_t_s(t,:)',[ones(length(pai),1) pai]);
        F_PLS(t,i)=beta(end);
    end
end
All_predictor=[predictor F_EW F_PC F_PLS];

% in-sample performance
N_All=size(All_predictor,2);
results_in_sample=nan(N_All,4);
for i=1:N_All
    OLS=nwest(y(R+1:end),[ones(P,1) zscore(All_predictor(R:end-1,i))],12);
    results_in_sample(i,1:4)=[100*OLS.beta(2) OLS.tstat(2) 1-normcdf(abs(OLS.tstat(2))) 100*OLS.rsqr];
end

% out-of-sample forecast
for t=1:P
    start=1;
    y_t=y(start:R+(t-1));
    predictor_t=predictor(start:R+(t-1),:);
    HA=mean(y(1:R+(t-1)));
    FC_HA(t)=HA;
    % original
    for i=1:N
        results_predictor_i_t=ols(y_t(2:end),[ones(length(y_t(2:end)),1) predictor_t(1:end-1,i)]);
        FC_original(t,i)=[1 predictor_t(end,i)]*results_predictor_i_t.beta;
    end
    % EW
    F_EW=mean(zscore(predictor_t),2);
    OLS=ols(y_t(2:end),[ones(length(y_t(2:end)),1) F_EW(1:end-1)]);
    FC_EW(t)=[1 F_EW(end)]*OLS.beta;
    % PC
    [~,F_PC]=pca(zscore(predictor_t),'corr');
    F_PC=F_PC(:,1);
    OLS=ols(y_t(2:end),[ones(length(y_t(2:end)),1) F_PC(1:end-1)]);
    FC_PC(t)=[1 F_PC(end)]*OLS.beta;
    % PLS
    predictor_t(:,[1 2 4])=detrend(predictor_t(:,[1 2 4]),1);
    y_comp=wavelet_decomposing_function(y_t,'haar',n_wd);
    y_comp=[y_t y_comp(:,end:-1:1)];
    predictor_t_s=zscore(predictor_t);
    for i=1:size(y_comp,2)
        F_PLS=nan(size(y_comp));
        pai=nan(N,1);
        for n=1:N
            predictor_t_s(:,n)=winsor(predictor_t_s(:,n),[2,98]);
            x_t=predictor_t_s(:,n);
            beta=regress(x_t(1:end-1),[ones(length(x_t(1:end-1)),1) y_comp(2:end,i)]);
            pai(n,:)=beta(2:end)';
        end
        for tt=1:R+(t-start)
            beta=regress(predictor_t_s(tt,:)',[ones(length(pai),1) pai]);
            F_PLS(tt,i)=beta(end);
        end
        OLS=ols(y_t(2:end),[ones(length(y_t(2:end)),1) F_PLS(1:end-1,i)]);
        FC_PLS(t,i)=[1 F_PLS(end,i)]*OLS.beta;
    end
    disp(t)
end

FC_F=[FC_EW FC_PC FC_PLS];

% R2OS performance
FC_all=[FC_original FC_F];
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

results_all_A=[results_in_sample R2OS];