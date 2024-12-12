clear;
addpath('jplv7')
input_file='data.xls';
input_sheet='Equity premium';
y=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:b1153');
input_sheet='Macroeconomic variables';
business_cycles=readmatrix(input_file,'Sheet',input_sheet,'Range','r2:r1153');
input_sheet='indpro';
indpro=readmatrix(input_file,'Sheet',input_sheet,'Range','b2:b1153');
n_wd=2;
y_comp=wavelet_decomposing_function(y,'haar',n_wd);
y_all=[y y_comp(:,end:-1:1)];
label={'Stock returns','Short-term component','Middle-term component','Long-term component'};
rec_index=find(business_cycles);

figure(1)
for i=1:size(y_all,2)
    subplot(size(y_all,2),1,i)
    for t=1:length(rec_index)
        plot([rec_index(t) rec_index(t)],[-100,100], 'color', [0.9,0.9,0.9], 'LineWidth', 1)
        hold on
    end
    plot(y_all(:,i));
    % ylabel(label(i));
    axis([0 size(y_all,1)+1 -0.25 0.25])
    set(gca,'XTick',[36 156 276 396 516 636 756 876 996 1116]);
    set(gca,'XTickLabel',{'1930','1940','1950','1960','1970','1980','1990','2000','2010','2020'});
    set(gcf,'color','w');
end

figure(2)
for i=1:size(y_all,2)
    subplot(2,2,i)
    x = indpro;
    y = y_all(:,i); 
    h = binscatter(x, y,"NumBins",[30 30]);
    colormap(gca, "summer");
    ylabel(label(i));
    xlabel('Industrial production growth');
    hold on; 
    p = polyfit(x,y,1); 
    yfit = polyval(p,x); 
    a = plot(x, yfit, 'r-', LineWidth=2);
    C = corrcoef(x, y);
    legend(a,['Correlation coefficient = ' sprintf('%.2f',C(2,1))],Location="southeast")
    set(gcf,'color','w');
end