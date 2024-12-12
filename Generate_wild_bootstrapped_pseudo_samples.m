clear;
warning('off','all');
addpath('jplv7')
input_file='data.xls';
input_sheet='Equity premium';
SP500=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:b1153');
input_sheet='Macroeconomic variables';
predictor=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:o1153');
input_sheet='sent';
sent=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:g685');
input_sheet='industry';
Industry=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:k1153');
input_sheet='size';
Size=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:k1153');
input_sheet='BM';
BM=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:k1153');
input_sheet='momentum';
Momentum=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:k1153');
input_sheet='DP_DG_EG';
DP_DG_EG=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:d1153');
y=[SP500 Industry Size BM Momentum];
[T,M]=size(y);
N=size(predictor,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating wild bootstrapped pseudo samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preliminaries
disp('Wild bootstrap preliminaries');
B=2000;
y_star=nan(T,M,B);
predictor_star=nan(T,N,B);

% Estimating bias-corrected AR processes for economic variables
disp('Computing bias-adjusted AR parameters');
AR_coefficients=nan(N,2);
v_hat_c=nan(T-1,N);
for i=1:N
    results_AR_i=ols(predictor(2:T,i),[ones(T-1,1) predictor(1:T-1,i)]);
    rho_hat_i=results_AR_i.beta(2);
    rho_hat_c_i=rho_hat_i+((1+3*rho_hat_i)/T)+(3*(1+3*rho_hat_i)/T^2);
    theta_hat_c_i=mean(predictor(2:T,i)-rho_hat_c_i*predictor(1:T-1,i));
    AR_coefficients(i,:)=[theta_hat_c_i rho_hat_c_i];
    v_hat_c(:,i)=predictor(2:T,i)-(theta_hat_c_i*ones(T-1,1)+...
        rho_hat_c_i*predictor(1:T-1,i));
end

% Generating wild bootstrapped draws for economic variables
disp('Generating wild bootstrapped draws');
rng('default');
w=randn(T-1,1,B);
for b=1:B
    w_b=w(:,:,b);
    v_star_b=kron(ones(1,N),w_b).*v_hat_c;
    v_star_b=[zeros(1,N) ; ...
        v_star_b-kron(mean(v_star_b),ones(T-1,1))];
    predictor_star_b=zeros(T,N);
    predictor_star_b(1,:)=predictor(1,:);
    for t=2:T
        predictor_star_b(t,:)=AR_coefficients(:,1)'+...
            AR_coefficients(:,2)'.*predictor_star_b(t-1,:)+v_star_b(t,:);
    end
    predictor_star(:,:,b)=predictor_star_b;
    disp(b);
end

% Generating wild bootstrapped draws for returns
for i=1:M
    disp('Computing OLS residuals');
    results_kitchen_sink=ols(y(2:T,i),...
        [ones(T-1,1) predictor(1:T-1,:)]);
    u_hat=results_kitchen_sink.resid;
    % Generating wild bootstrapped draws
    disp('Generating wild bootstrapped draws');
    alpha_hat_restrict=mean(y(2:T,i));
    rng('default');
    w=randn(T-1,1,B);
    for b=1:B
        w_b=w(:,:,b);
        u_star_b=w_b.*u_hat;
        u_star_b=[0 ; u_star_b-mean(u_star_b)];
        y_star_b=zeros(T,1);
        y_star_b(1)=y(1);
        for t=2:T
            y_star_b(t)=alpha_hat_restrict+u_star_b(t);
        end
        y_star(:,i,b)=y_star_b;
        disp(b);
    end
end

% Generating wild bootstrapped draws for factors
n_wd=2;
F_predictor_star=nan(T,6,B);
for b=1:B
    % EW
    F_EW=mean(zscore(predictor_star(:,:,b)),2);
    % PC
    [~,F_PC]=pca(zscore(predictor_star(:,:,b)));
    F_PC=F_PC(:,1);
    % PLS
    F_PLS=zeros(T,4);
    for t=2:T
        y_t=y_star(1:t,1,b);
        predictor_t=predictor_star(1:t,:,b);
        y_comp=wavelet_decomposing_function(y_t,'haar',n_wd);
        y_comp=[y_t y_comp(:,end:-1:1)];
        predictor_t_s=zscore(predictor_t);
        for i=1:size(y_comp,2)
            pai=nan(N,1);
            for n=1:N
                x_t=predictor_t_s(:,n);
                beta=regress(x_t(1:end-1),[ones(length(x_t(1:end-1)),1) y_comp(2:end,i)]);
                pai(n)=beta(end);
            end
            beta=regress(predictor_t_s(t,:)',[ones(length(pai),1) pai]);
            F_PLS(t,i)=beta(end);
        end
    end
    F_predictor_star(:,:,b)=[F_EW F_PC F_PLS];
    disp(b);
end

% Estimating bias-corrected AR processes for sent variables
[T,N]=size(sent);
sent_star=nan(T,N,B);
disp('Computing bias-adjusted AR parameters');
AR_coefficients=nan(N,2);
v_hat_c=nan(T-1,N);
for i=1:N
    results_AR_i=ols(sent(2:T,i),[ones(T-1,1) sent(1:T-1,i)]);
    rho_hat_i=results_AR_i.beta(2);
    rho_hat_c_i=rho_hat_i+((1+3*rho_hat_i)/T)+(3*(1+3*rho_hat_i)/T^2);
    theta_hat_c_i=mean(sent(2:T,i)-rho_hat_c_i*sent(1:T-1,i));
    AR_coefficients(i,:)=[theta_hat_c_i rho_hat_c_i];
    v_hat_c(:,i)=sent(2:T,i)-(theta_hat_c_i*ones(T-1,1)+...
        rho_hat_c_i*sent(1:T-1,i));
end

% Generating wild bootstrapped draws for sent variables
disp('Generating wild bootstrapped draws');
rng('default');
w=randn(T-1,1,B);
for b=1:B
    w_b=w(:,:,b);
    v_star_b=kron(ones(1,N),w_b).*v_hat_c;
    v_star_b=[zeros(1,N) ; ...
        v_star_b-kron(mean(v_star_b),ones(T-1,1))];
    sent_star_b=zeros(T,N);
    sent_star_b(1,:)=sent(1,:);
    for t=2:T
        sent_star_b(t,:)=AR_coefficients(:,1)'+...
            AR_coefficients(:,2)'.*sent_star_b(t-1,:)+v_star_b(t,:);
    end
    sent_star(:,:,b)=sent_star_b;
    disp(b);
end

% Estimating bias-corrected AR processes for DP_DG_EG variables
[T,N]=size(DP_DG_EG);
DP_DG_EG_star=nan(T,N,B);
disp('Computing bias-adjusted AR parameters');
AR_coefficients=nan(N,2);
v_hat_c=nan(T-1,N);
for i=1:N
    results_AR_i=ols(DP_DG_EG(2:T,i),[ones(T-1,1) DP_DG_EG(1:T-1,i)]);
    rho_hat_i=results_AR_i.beta(2);
    rho_hat_c_i=rho_hat_i+((1+3*rho_hat_i)/T)+(3*(1+3*rho_hat_i)/T^2);
    theta_hat_c_i=mean(DP_DG_EG(2:T,i)-rho_hat_c_i*DP_DG_EG(1:T-1,i));
    AR_coefficients(i,:)=[theta_hat_c_i rho_hat_c_i];
    v_hat_c(:,i)=DP_DG_EG(2:T,i)-(theta_hat_c_i*ones(T-1,1)+...
        rho_hat_c_i*DP_DG_EG(1:T-1,i));
end

% Generating wild bootstrapped draws for DP_DG_EG variables
disp('Generating wild bootstrapped draws');
rng('default');
w=randn(T-1,1,B);
for b=1:B
    w_b=w(:,:,b);
    v_star_b=kron(ones(1,N),w_b).*v_hat_c;
    v_star_b=[zeros(1,N) ; ...
        v_star_b-kron(mean(v_star_b),ones(T-1,1))];
    DP_DG_EG_star_b=zeros(T,N);
    DP_DG_EG_star_b(1,:)=DP_DG_EG(1,:);
    for t=2:T
        DP_DG_EG_star_b(t,:)=AR_coefficients(:,1)'+...
            AR_coefficients(:,2)'.*DP_DG_EG_star_b(t-1,:)+v_star_b(t,:);
    end
    DP_DG_EG_star(:,:,b)=DP_DG_EG_star_b;
    disp(b);
end

output_file='Generate_wild_bootstrapped_pseudo_samples';
save(output_file,'y_star','predictor_star','F_predictor_star','sent_star','DP_DG_EG_star');
