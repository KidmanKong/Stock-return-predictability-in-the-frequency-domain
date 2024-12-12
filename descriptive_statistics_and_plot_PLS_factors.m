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
rec_index=find(business_cycles(R+1:R+P));

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
All_variable=[y*100 predictor F_EW F_PC F_PLS];

% descriptive stat
summary=nan(size(All_variable,2),7);
for i=1:size(All_variable,2)
    summary(i,1)=mean(All_variable(:,i));
    summary(i,2)=std(All_variable(:,i));
    summary(i,3)=min(All_variable(:,i));
    summary(i,4)=max(All_variable(:,i));
    summary(i,5)=skewness(All_variable(:,i));
    summary(i,6)=kurtosis(All_variable(:,i));
    a=autocorr(All_variable(:,i));
    summary(i,7)=a(2);
end    

% variance contribution
y_comp=wavelet_decomposing_function(y,'haar',6);
y_comp=y_comp(:,end:-1:1);
var_ratio=100*var(y_comp)./var(y);

% plot PLS factors
figure(3)
labels={'PLS-R','PLS-SC','PLS-MC','PLS-LC'};
for i=1:size(F_PLS,2)
    subplot(2,2,i)
    for t=1:length(rec_index)
        plot([rec_index(t) rec_index(t)],[-100,100], 'color', [0.9,0.9,0.9], 'LineWidth', 1)
        hold on
    end
    a=plot(zscore(F_PLS(R+1:end,i)),'LineWidth', 1);
    axis([0 P+1 -4 4])
    legend(a,labels{i},Location="northwest")
    set(gca,'XTick',[36 156 276 396 516 636 756 876]);
    set(gca,'XTickLabel',{'1950','1960','1970','1980','1990','2000','2010','2020'});
    set(gcf,'color','w');
end