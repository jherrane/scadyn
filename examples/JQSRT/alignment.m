clearvars -except pth; clc; close all;
if(exist('pth') == 0)
   pth = input('gib path: ');
   pth = setpath(pth);
end
files = {'all','sph','ob','pro','ell'};
xs = {'0.3', '1', '3'};
titles = {'Whole ensemble',...
    'Spheres',...
    'Oblate spheroids',...
    'Prolate spheroids',...
    'Ellipsoids'};
fig = figure;
pos = [0,0,1800,900];
fs = 14;

nbins = 15;

set(fig,'Position',pos);
i = 1;
for x = 1:length(xs)
    for file = 1:length(files)
        suffix = files{file};

        fn = [pth,'alignment-',xs{x},'-',files{file}];
        loaded = importdata(fn);
        data = loaded.data(:,1);
        b = loaded.data(:,2);
        proper = find(b>2);
        pseudo = find(b<2);   
        acwq = real(acos(data));
        acwq3 = real(acos(data(proper)));
        acwq1 = real(acos(data(pseudo)));

        subplot(3,5,i)
        countsprop = histcounts(acwq3,nbins);
        countsshit = histcounts(acwq1,nbins);
        counts = countsprop + countsshit;
%         histogram(180*acwq/pi,nbins)
        bar(90/nbins:90/nbins:90,counts,'FaceColor', [0    0.4470    0.7410]);
        hold on;
        axis([0 90 0 inf])
        bar(90/nbins:90/nbins:90,countsprop,'FaceColor', [0.8500    0.3250    0.0980]);
        if x == 1 && file == 1
            h = legend('$\hat{Q}_1$','$\hat{Q}_3$');
            set(h,'Interpreter','latex')
            h.FontSize = fs;
        end
        NumTicks = 4;
        L = get(gca,'XLim');
        set(gca,'XTick',linspace(L(1),L(2),NumTicks))
        
        if x == 1
            title(titles{file},'FontSize',fs);
        end
        if file == 1
           ylabel(['ka = ',xs{x}],'FontSize',fs+2);
        end
        xlabel('$\angle(\hat{\omega}_{av},\hat{Q}_{1/3,av})$ (deg)','Interpreter','latex',...
            'FontSize',fs)
        hold on;
        i = i + 1;
    end
end
print('alignment','-dpng')