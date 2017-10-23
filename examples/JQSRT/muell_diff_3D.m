clearvars -except pth; clc; close all;
if(exist('pth') == 0)
   pth = input('gib path: ');
   pth = setpath(pth);
end
refr={'1.31', '1.70', '2.0'};
refi={'0.0','0.0','0.2'};
x={'3','1','0.3'};
shapes={'sph','ob','pro','ell'};

ds = {'',''};

% Parameters for plotting
fs = 12;
lw = 2;
ps = [0,0,1000,1000];
xtick = [0,60,120,180,240,300,360];
ytick = [0,30,60,90,120,150,180];
ylab = '\theta';
xlab = '\phi';
zlabs = {'S_{11}', '-S_{12}/S_{11}', 'S_{13}/S_{11}', 'S_{14}/S_{11}'};

bas = 'mueller';
for ish = 1:4
    base = [bas,'-',shapes{ish}];

    s12lim = [-1.1,1.1];

    for irf = 1:3
        for ix = 1:3
            close;
            fig = figure;
            set(fig, 'Position', ps)

            for ids = 1:1
                suffix = [refr{irf},'-',refi{irf},'-x-',x{ix}];
                figname = [base,'-',suffix,'.png'];

                dat = importdata([pth,base,'-',ds{ids},suffix]);
                phi = dat.data(:,2);
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
                s112d = datum(:,1); s112d = reshape(s112d,[180,90]);
                s112d = s112d(:,1:90);

                s122d = datum(:,2); s122d = reshape(s122d,[180,90]);
                s122d = s122d(:,1:90);
                s132d = datum(:,3); s132d = reshape(s132d,[180,90]);
                s132d = s132d(:,1:90);
                s142d = datum(:,4); s142d = reshape(s142d,[180,90]);
                s142d = s142d(:,1:90);
                datum2 = zeros(180,90,3);
                datum2(:,:,1) = s122d; datum2(:,:,2) = s132d;
                datum2(:,:,3) = s142d;

                [X,Y] = meshgrid(0:4:359,0:1:179);

                if s13mx>0.01
                    s13lim = 1.1*[-s13mx,s13mx];
                else
                    s13lim = s12lim;
                end
                if s14mx>0.01
                    s14lim = 1.1*[-s14mx,s14mx];
                else
                    s14lim = s12lim;
                end
                yl = [s12lim;s13lim;s14lim];
                clim = [0, 1000];

                for i = 1:4
                    if i == 1
                        ax1 = subplot(2,2,1);
                        s = surf(X,Y,s112d,'FaceLighting', ...
                            'gouraud','FaceColor','interp');             
                        set(gca,'ylim',[0 180],'xlim',[0 360],'zlim',[0 inf])
                        set(gca,'zscale','log','xscale','linear','yscale','linear');
                        light
                        colormap(ax1,cool)
    %                     caxis(log10(clim))
                        view(120,30)
                    else
                        ax2 = subplot(2,2,i);
                        s = surf(X,Y,datum2(:,:,i-1),'FaceLighting','gouraud','FaceColor','interp');
                        set(gca,'ylim',[0 180],'xlim',[0 360],'zlim',yl(i-1,:))
                        light
                        colormap(ax2,parula)
    %                     set(gca,'Ydir','reverse')
                        view(120,30)
                    end
                    hold on;
                    material dull
                    s.EdgeColor = 'none';
                    p = get(gca,'position');
                    p(4) = p(4)*1.2; 
                    p(3) = p(3)*1.08; 
                    set(gca, 'position', p);

                    set(gca,'fontsize',fs);
                    xlabel(xlab)
                    ylabel(ylab)
                    zlabel(zlabs{i})
                    xticks(xtick)
                    yticks(ytick)
                    if ids==2 && i==1
                        legend('Random orientation','Scattering plane 1');
                    end
                end
            end

            print(figname,'-dpng')
        end
    end
end

close;