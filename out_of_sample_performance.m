clear;
addpath('jplv7')
input_file='data.xls';
input_sheet='Equity premium';
y=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:b1153');
input_sheet='Macroeconomic variables';
predictor=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:o1153');
business_cycles=readmatrix(input_file,'Sheet',input_sheet,'Range','r2:r1153');
Rfree_lag=readmatrix(input_file,'Sheet',input_sheet,'Range','q2:q1153');
T=size(y,1);
R=240;
P=T-R;
N=size(predictor,2);
n_wd=2;
actual=y(R+1:R+P);
FC_HA=nan(P,1);
FC_original=nan(P,N);
FC_EW=nan(P,1);
FC_PC=nan(P,1);
FC_SPC=nan(P,1);
FC_PLS=nan(P,n_wd+2);
rec_index=find(business_cycles(R+1:R+P));
exp_index=find(~business_cycles(R+1:R+P));

%% out-of-sample forecast
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
    results_EW_t=ols(y_t(2:end),[ones(length(y_t(2:end)),1) F_EW(1:end-1)]);
    FC_EW(t)=[1 F_EW(end)]*results_EW_t.beta;
    % PC
    [~,F_PC]=pca(zscore(predictor_t),'corr');
    F_PC=F_PC(:,1);
    results_PC_t=ols(y_t(2:end),[ones(length(y_t(2:end)),1) F_PC(1:end-1)]);
    FC_PC(t)=[1 F_PC(end)]*results_PC_t.beta;
    % PLS
    predictor_t(:,[1 2 4])=detrend(predictor_t(:,[1 2 4]),1);
    predictor_t_s=zscore(predictor_t);
    y_comp=wavelet_decomposing_function(y_t,'haar',n_wd);
    y_comp=[y_t y_comp(:,end:-1:1)];

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
        FC_PLS(t,i)=[1 F_PLS(end,i)]*OLS.beta;
    end
    disp(t)
end

FC_F=[FC_EW FC_PC FC_PLS];
DSPE_F=nan(size(FC_F));
for t=1:size(FC_F,1)
    for i=1:size(FC_F,2)
        DSPE_F(t,i)=(FC_HA(t)-actual(t))^2-(FC_F(t,i)-actual(t))^2;
    end
end
CDSPE_F=cumsum(DSPE_F);

%% plot
figure(4)
for t=1:length(rec_index)
    plot([rec_index(t) rec_index(t)],[-100,100], 'color', [0.9,0.9,0.9], 'LineWidth', 1)
    hold on
end
a=plot([actual FC_EW FC_PC FC_PLS(:,end)],'linewidth',0.5);
xlim([0 P+1])
ylim([-0.1,0.1])
legend(a,'Stock returns','EW forecast','PC forecast','PLS-LC forecast')
set(gca,'XTick',[36 156 276 396 516 636 756 876]);
set(gca,'XTickLabel',{'1950','1960','1970','1980','1990','2000','2010','2020'});
set(gcf,'color','w');

figure(5)
for t=1:length(rec_index)
    plot([rec_index(t) rec_index(t)],[-100,100], 'color', [0.9,0.9,0.9], 'LineWidth', 1)
    hold on
end
a=plot(CDSPE_F,'linewidth',2);
hold on
plot(zeros(P,1),'k--','linewidth',2)
xlim([0 P+1])
ylim([-0.06,0.06])
legend(a,'EW','PC','PLS-R','PLS-SC','PLS-MC','PLS-LC')
set(gca,'XTick',[36 156 276 396 516 636 756 876]);
set(gca,'XTickLabel',{'1950','1960','1970','1980','1990','2000','2010','2020'});
set(gcf,'color','w');

%% forecast encompassing test
FC_all=[FC_original FC_F];
encompassing_test=nan(size(FC_all,2),size(FC_F,2));
for i=1:size(FC_F,2)
    for n=1:size(FC_all,2)
    [lambda,~,MHLN_pval]=Perform_HLN_test(actual,FC_all(:,n),FC_F(:,i));
    encompassing_test(n,i)=MHLN_pval;
    end
end

%% R2OS performance
R2OS=nan(size(FC_all,2),5);
% MSFE criterion, historical average
e_HA=(actual-FC_HA).^2;
MSFE_HA=mean(e_HA);
% MSFE criterion, predictor
for i=1:size(FC_all,2)
     MSFE_i=mean((actual-FC_all(:,i)).^2);
     R2OS_i=100*(1-(MSFE_i/MSFE_HA));
     [MSFE_adjusted_i,p_value_i_CW]=Perform_CW_test(actual,FC_HA,FC_all(:,i));
     e_i=(actual-FC_all(:,i)).^2;
     [DM_i,p_value_i_DM]=dmtest(e_HA.^0.5,e_i.^0.5,1);
    R2OS(i,:)=[R2OS_i MSFE_adjusted_i p_value_i_CW DM_i,p_value_i_DM];   
end

actual_rec=actual(rec_index);
actual_exp=actual(exp_index);
FC_HA_rec=FC_HA(rec_index);
FC_HA_exp=FC_HA(exp_index);

FC_all_rec=FC_all(rec_index,:);
R2OS_rec=nan(size(FC_all_rec,2),3);
% MSFE criterion, historical average
MSFE_HA_res=mean((actual_rec-FC_HA_rec).^2);
% MSFE criterion, predictor
for i=1:size(FC_all_rec,2)
     MSFE_i=mean((actual_rec-FC_all_rec(:,i)).^2);
     R2OS_i=100*(1-(MSFE_i/MSFE_HA_res));
     [MSFE_adjusted_i,p_value_i]=Perform_CW_test(actual_rec,FC_HA_rec,FC_all_rec(:,i));
     R2OS_rec(i,:)=[R2OS_i p_value_i MSFE_i];   
end
FC_all_exp=FC_all(exp_index,:);
R2OS_exp=nan(size(FC_all_exp,2),3);
% MSFE criterion, historical average
MSFE_HA_exp=mean((actual_exp-FC_HA_exp).^2);
% MSFE criterion, predictor
for i=1:size(FC_all_exp,2)
     MSFE_i=mean((actual_exp-FC_all_exp(:,i)).^2);
     R2OS_i=100*(1-(MSFE_i/MSFE_HA_exp));
     [MSFE_adjusted_i,p_value_i]=Perform_CW_test(actual_exp,FC_HA_exp,FC_all_exp(:,i));
     R2OS_exp(i,:)=[R2OS_i p_value_i MSFE_i];   
end
R2OS_all=[R2OS R2OS_exp(:,1) R2OS_rec(:,1)];

%% portfolio performance
VOL_window=60;
FC_VOL=nan(P,1);
for t=1:P
    % Volatility forecast
    FC_VOL(t)=mean(y(R+(t-1)-VOL_window+1:R+(t-1)).^2)-...
        (mean(y(R+(t-1)-VOL_window+1:R+(t-1))))^2;
end
r_f_lag_P=Rfree_lag(R+1:R+P);

C=[0 50]; % transaction cost
G=[3 5]; % risk aversion coefficient
portfolio_all=nan(size(FC_all,2),2*(length(C)+length(G)));

for ci=1:length(C)
    for gi=1:length(G)
        c_bp=C(ci);
        gamma_MV=G(gi);     
        % Computing average utility: historical average
        [v_HA,SR_HA,xxx,xxx,TO_HA]=Perform_asset_allocation(actual,r_f_lag_P,FC_HA,FC_VOL,gamma_MV,c_bp);
        portfolio_HA=[1200*v_HA SR_HA 100*TO_HA];
        % Computing average utility gains
        for i=1:size(FC_all,2)
            [v_i,SR_i,xxx,xxx,TO_i]=Perform_asset_allocation(actual,r_f_lag_P,FC_all(:,i),FC_VOL,gamma_MV,c_bp);
            portfolio_all(i,4*(ci-1)+2*(gi-1)+1:4*(ci-1)+2*(gi-1)+2)=[1200*(v_i-v_HA) SR_i];
        end
    end
end
