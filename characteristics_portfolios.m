clear;
addpath('jplv7')
input_file='data.xls';
input_sheet='Equity premium';
Y=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:b1153');
input_sheet='industry';
Industry=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:k1153');
input_sheet='size';
Size=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:k1153');
input_sheet='BM';
BM=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:k1153');
input_sheet='momentum';
Momentum=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:k1153');
input_sheet='Macroeconomic variables';
predictor=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:o1153');
y=[Industry Size BM Momentum];
[T,M]=size(y);
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
results_EW=nan(M,6);
results_PC=nan(M,6);
results_PLS=nan(M,6);
n_PLS=4;

%% in-sample regression
for i=1:M
    % EW
    OLS=nwest(y(R+1:end,i),[ones(P,1) zscore(F_EW(R:end-1))],12);
    results_EW(i,1:2)=[OLS.beta(2) OLS.tstat(2)];
    results_EW(i,4)=100*OLS.rsqr;
    % PC
    OLS=nwest(y(R+1:end,i),[ones(P,1) zscore(F_PC(R:end-1))],12);
    results_PC(i,1:2)=[OLS.beta(2) OLS.tstat(2)];
    results_PC(i,4)=100*OLS.rsqr;
    % PLS
    OLS=nwest(y(R+1:end,i),[ones(P,1) zscore(F_PLS(R:end-1,n_PLS))],12);
    results_PLS(i,1:2)=[OLS.beta(2) OLS.tstat(2)];
    results_PLS(i,4)=100*OLS.rsqr;
end

% Computing wild bootstrapped p-values
input_file='Generate_wild_bootstrapped_pseudo_samples';
load(input_file);
B=size(y_star,2); % number of boostrap replications
y_star=y_star(:,2:M+1,:);
tstat_EW_star=nan(B,M);
tstat_PC_star=nan(B,M);
tstat_PLS_star=nan(B,M);
for i=1:M
    for b=1:B
        y_b=y_star(:,i,b);
        F_EW_b=F_predictor_star(:,1,b);
        F_PC_b=F_predictor_star(:,2,b);
        F_PLS_b=F_predictor_star(:,2+n_PLS,b);
        % EW
        OLS=nwest(y_b(R+1:end),[ones(P,1) F_EW_b(R:end-1)],12);
        tstat_EW_star(b,i)=OLS.tstat(2);
        % PC
        OLS=nwest(y_b(R+1:end),[ones(P,1) F_PC_b(R:end-1)],12);
        tstat_PC_star(b,i)=OLS.tstat(2);
        % PLS
        OLS=nwest(y_b(R+1:end),[ones(P,1) F_PLS_b(R:end-1)],12);
        tstat_PLS_star(b,i)=OLS.tstat(2);
    end
    disp(i);
end

for i=1:M
    results_EW(i,3)=mean(abs(tstat_EW_star(:,i))>abs(results_EW(i,2)));
    results_PC(i,3)=mean(abs(tstat_PC_star(:,i))>abs(results_PC(i,2)));
    results_PLS(i,3)=mean(abs(tstat_PLS_star(:,i))>abs(results_PLS(i,2)));
end

%% out-of-sample
actual=y(R+1:R+P,:);
FC_HA=nan(P,M);
FC_EW=nan(P,M);
FC_PC=nan(P,M);
FC_PLS=nan(P,M);

for t=1:P
    start=1;
    predictor_t=predictor(start:R+(t-1),:);
    % EW
    F_EW=mean(zscore(predictor_t),2);
    % PC
    [~,F_PC]=pca(zscore(predictor_t));
    F_PC=F_PC(:,1);
    % PLS
    Y_t=Y(start:R+(t-1));
    predictor_t(:,[1 2 4])=detrend(predictor_t(:,[1 2 4]),1);
    y_comp=wavelet_decomposing_function(Y_t,'haar',n_wd);
    y_comp=y_comp(:,1);
    predictor_t_s=zscore(predictor_t);
    F_PLS=nan(size(y_comp));
    pai=nan(N,1);
    for n=1:N
        predictor_t_s(:,n)=winsor(predictor_t_s(:,n),[2 98]);
        x_t=predictor_t_s(:,n);
        beta=regress(x_t(1:end-1),[ones(length(x_t(1:end-1)),1) y_comp(2:end)]);
        pai(n,:)=beta(2:end)';
    end
    for tt=1:R+(t-start)
        beta=regress(predictor_t_s(tt,:)',[ones(length(pai),1) pai]);
        F_PLS(tt)=beta(end);
    end
    % [~,index]=winsor(F_PLS,[0,100]);
    % F_PLS(index)=mean(F_PLS);
    for i=1:M
        y_t=y(start:R+(t-1),i);
        HA=mean(y(1:R+(t-1),i));
        FC_HA(t,i)=HA;
        % EW
        OLS=ols(y_t(2:end),[ones(length(y_t(2:end)),1) F_EW(1:end-1)]);
        FC_EW(t,i)=[1 F_EW(end)]*OLS.beta;
        % PC
        OLS=ols(y_t(2:end),[ones(length(y_t(2:end)),1) F_PC(1:end-1)]);
        FC_PC(t,i)=[1 F_PC(end)]*OLS.beta;
        % PLS
        OLS=ols(y_t(2:end),[ones(length(y_t(2:end)),1) F_PLS(1:end-1)]);
        FC_PLS(t,i)=[1 F_PLS(end)]*OLS.beta;
    end   
    disp(t)
end

for i=1:M
    % MSFE criterion, historical average
    e_HA=(actual(:,i)-FC_HA(:,i)).^2;
    MSFE_HA=mean(e_HA);
    % MSFE EW
    MSFE_EW=mean((actual(:,i)-FC_EW(:,i)).^2);
    R2OS_EW=100*(1-(MSFE_EW/MSFE_HA));
    [~,p_value_CW]=Perform_CW_test(actual(:,i),FC_HA(:,i),FC_EW(:,i));
    results_EW(i,5:6)=[R2OS_EW p_value_CW];   
    % MSFE PC
    MSFE_PC=mean((actual(:,i)-FC_PC(:,i)).^2);
    R2OS_PC=100*(1-(MSFE_PC/MSFE_HA));
    [~,p_value_CW]=Perform_CW_test(actual(:,i),FC_HA(:,i),FC_PC(:,i));
    results_PC(i,5:6)=[R2OS_PC p_value_CW];
    % MSFE PLS
    MSFE_PLS=mean((actual(:,i)-FC_PLS(:,i)).^2);
    R2OS_PLS=100*(1-(MSFE_PLS/MSFE_HA));
    [~,p_value_CW]=Perform_CW_test(actual(:,i),FC_HA(:,i),FC_PLS(:,i));
    results_PLS(i,5:6)=[R2OS_PLS p_value_CW];
end