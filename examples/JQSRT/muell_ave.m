clear; clc; close all;
refr={'1.31', '1.70', '2.0'};
refi={'0.0','0.0','0.2'};
x={'3','1','0.3'};
shapes={'sph','ob','pro','ell'};

ds = {'0',''};

% Parameters for plotting
fs = 12;
lw = [2.5,1.5];
ps = [0,0,1000,1000];
xtick = [0,30,60,90,120,150,180];
xlab = 'Scattering angle';
ylabs = {'S_{11}', '-S_{12}/S_{11}', 'S_{13}/S_{11}', 'S_{14}/S_{11}'};

base = 'mueller';
path = 'mueller/';

s12lim = [-1.1,1.1];

for shape = 1:4
    for irf = 1:3
        for ix = 1:3
            close;
            fig = figure;
            set(fig, 'Position', ps)

            for ids = 1:2
                suffix = [refr{irf},'-',refi{irf},'-x-',x{ix}];
                figname = ['ave','-',shapes{shape},base,'-',suffix,'.png'];
%                 dat = importdata([path,base,'_all',ds{ids},'-',suffix]);

                dat = importdata([path,base,'-',shapes{shape},ds{ids},'-',suffix]);
                if ids == 1
                    theta = dat.data(:,2);
                else
                    theta = dat.data(:,2);
                end
                S11 = dat.data(:,3); S12 = dat.data(:,4);
                S13 = dat.data(:,5); S14 = dat.data(:,6);
                datum = [S11/S11(30), -S12./S11, S13./S11, S14./S11];
                s13mx = max(abs(max(datum(:,3))),abs(min(datum(:,3))));
                s14mx = max(abs(max(datum(:,4))),abs(min(datum(:,4))));

                if s13mx>0.05
                    s13lim = 1.1*[-s13mx,s13mx];
                else
                    s13lim = s12lim;
                end
                if s14mx>0.05
                    s14lim = 1.1*[-s14mx,s14mx];
                else
                    s14lim = s12lim;
                end
                yl = [s12lim;s13lim;s14lim];

                for i = 1:4
                    subplot(2,2,i)

                    if i == 1
                        if ids == 1
                            semilogy(theta*180/pi,datum(:,i), 'linewidth', lw(ids), 'Color', [0    0.4470    0.7410])
                        else
                            semilogy(theta*180/pi,datum(:,i), 'linewidth', lw(ids), 'Color', [0.8500    0.3250    0.0980])
                        end
                    else
                        if ids == 1
                            plot(theta*180/pi,datum(:,i), 'linewidth', lw(ids), 'Color', [0    0.4470    0.7410])
                        else
                            plot(theta*180/pi,datum(:,i), 'linewidth', lw(ids), 'Color', [0.8500    0.3250    0.0980])
                        end

                        ylim(yl(i-1,:))
                    end
                    hold on;
                    set(gca,'fontsize',fs);
                    xlabel(xlab)
                    ylabel(ylabs{i})
                    xlim([0,180])
                    xticks(xtick)
                    if ids==2 && i==1
                        legend('Smooth RO','Gaussian shape RO');
                    end
                end
            end

            print(figname,'-dpng')
        end
    end
end

close;