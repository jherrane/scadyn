%% COMB
clear; close all; clc;

Qgrs = importdata('Qgrs.dat');
xgrs = Qgrs.data(:,1);
Qfgrs = Qgrs.data(:,2);

Qgrsfine = importdata('Qgrsfine.dat');
xgrsfine = Qgrsfine.data(:,1);
Qfgrsfine = Qgrsfine.data(:,2);

ps = [0,0,1200,1200];
fig1  = figure(1);
lw = 4;
fsize = 24;
ms = 48;
set(fig1, 'Position', ps)
set(gca,'fontsize',fsize);

hold on;

scatter(xgrs,Qfgrs,3*ms,'MarkerFaceColor','m','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
scatter(xgrsfine,Qfgrsfine,1.5*ms,'MarkerFaceColor','r','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
dl1.Color = [0,0,0,0.5];
dl1.LineStyle = '-.';

axis([min(xgrs)*(1-abs(min(xgrs))*0.2) max(xgrs)*(1+abs(min(xgrs))*0.1) 0 1.5])
xlabel('x')
ylabel('Q_{pr}')
legend('Normal tetras','Finer tetras','Location','northwest');
hold off;

%% TORQUES
QT1 = Qgrs.data(:,3);
QT2 = Qgrs.data(:,4);
QT3 = Qgrs.data(:,5);

QT1fine = Qgrsfine.data(:,3);
QT2fine = Qgrsfine.data(:,4);
QT3fine = Qgrsfine.data(:,5);
xfine = Qgrsfine.data(:,1);

ps = [0,0,1200,1200];
fig4  = figure(4);
lw = 2;
fsize = 24;
ms = 48;
set(fig4, 'Position', ps)
set(gca,'fontsize',fsize);

hold on;

scatter(xgrs,QT1,3*ms,'MarkerFaceColor',[1,0.4,0.5],'MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
scatter(xgrsfine,QT1fine,1.5*ms,'MarkerFaceColor',[1,0.6,0.7],'MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
scatter(xgrs,QT2,3*ms,'MarkerFaceColor',[0.3,0.4,1],'MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
scatter(xgrsfine,QT2fine,1.5*ms,'MarkerFaceColor',[0.6,0.7,1],'MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
scatter(xgrs,QT3,3*ms,'MarkerFaceColor',[0.1,1,0.3],'MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
scatter(xgrsfine,QT3fine,1.5*ms,'MarkerFaceColor',[0.3,1,0.5],'MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)




axis([min(xgrs) max(xgrs) -inf inf])
xlabel('x')
ylabel('Q_{\Gamma}')
legend('Q_{\Gamma, 1}','Q_{\Gamma, 1 (Finer)}','Q_{\Gamma, 2}',...
    'Q_{\Gamma, 2 (Finer)}','Q_{\Gamma, 3}', 'Q_{\Gamma, 3 (Finer)}','Location','northwest')

hold off;
