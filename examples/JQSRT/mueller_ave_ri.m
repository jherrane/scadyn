clear; clc; close all;
refr={'1.31', '1.70', '2.0'};
refi={'0.0','0.0','0.2'};
x={'3','1','0.3'};
shapes={'sph','ob','pro','ell'};

ds = {'0',''};

% Parameters for plotting
fs = 20;
lw = [2.5,1.5];
ps = [0,0,1500,650];
xtick = [0,30,60,90,120,150,180];
xlab = 'Scattering angle';
ylabs = {'log_{10} S_{11}', '-S_{12}/S_{11}', 'S_{13}/S_{11}', 'S_{14}/S_{11}'};
lspec = {'-','-';'--','--';':',':'};

base = 'mueller';
path = 'mueller/';

s12lim = [-1.1,1.1];

close;
fig = figure;
set(fig, 'Position', ps)
figname = ['ave',base,'.png'];

a = {};
A = zeros(3);
ii = 1;
iii = 0;

val = [-inf,-inf,-inf];
ind = [1,1,1];
for ix = 1:3
    for irf = 1:3
        for ids = 1:2
            suffix = ['x-',x{ix}];
            dat = importdata([path,base,'_all',ds{ids},'-',refr{irf},'-',refi{irf},'-',suffix]);
            theta = dat.data(:,2);
            
            S11 = dat.data(:,3); S12 = dat.data(:,4);
            S13 = dat.data(:,5); S14 = dat.data(:,6);
            
            if ids==2
                iii = iii + 1;
                a{iii,1} = [refr{irf},'-',refi{irf},'-',suffix];
                a{iii,2} = log10(get_Csca(S11,180)/(15*6));
                A(irf,4-ix) = a{iii,2};
            end

            datum = [log10(S11/S11(30)), -S12./S11, S13./S11, S14./S11];
            s13lim = 1.1*[-0.1,0.5];
            s14lim = 1.1*[-0.1,0.5];
            yl = [s12lim;s13lim;s14lim]; 

            for i = 1:2
                subplot(1,2,i)
                if i == 1
                    if ids == 1
                        p(ii) = plot(theta*180/pi,datum(:,i)+(ix-1)*2, lspec{irf,ids},'linewidth', lw(ids), 'Color', [0    0.4470    0.7410]);
                    else
                        p(ii) = plot(theta*180/pi,datum(:,i)+(ix-1)*2, lspec{irf,ids},'linewidth', lw(ids), 'Color', [0.8500    0.3250    0.0980]);
                    end
                    txt = ['ka = ', x{ix}];
                    if max(datum(:,i))>val(ix)
                        [val(ix),ind(ix)] = max(datum(:,i));
                    end
                    if ids == 2 && irf == 3
                        text(theta(ind(ix))*180/pi,val(ix)+(ix-1)*2,txt,'FontSize',fs,'VerticalAlignment','bottom');
                    end
                else
                    if ids == 1
                        p(ii) = plot(theta*180/pi,datum(:,i)+(ix-1)*2, lspec{irf,ids},'linewidth', lw(ids), 'Color', [0    0.4470    0.7410]);
                    else
                        p(ii) = plot(theta*180/pi,datum(:,i)+(ix-1)*2, lspec{irf,ids},'linewidth', lw(ids), 'Color', [0.8500    0.3250    0.0980]);
                    end
                end
                ii = ii + 1;
                hold on;
                grid on;
                grid minor
                set(gca,'fontsize',fs);
                xlabel(xlab)
                ylabel(ylabs{i})
                xlim([0,180])
                xticks(xtick)
                if ids==2 && i==1 && irf == 3 && ix == 3
                    lgnd = legend([p(1),p(3),p(5),p(7),p(9),p(11)],...
                        ['S','-',refr{1},'-',refi{1}], ...
                        ['GE','-',refr{1},'-',refi{1}],...
                        ['S','-',refr{2},'-',refi{2}], ...
                        ['GE','-',refr{2},'-',refi{2}],...
                        ['S','-',refr{3},'-',refi{3}], ...
                        ['GE','-',refr{3},'-',refi{3}], ...
                    'Location','northeast');
                    set(lgnd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.75])); 
                end
            end
        end
    end
end
print(figname,'-dpng')
A2 = [0, 0.3, 1, 3; 1.31, A(1,:); 1.70, A(2,:); 2.0, A(3,:)];
d1 = digits(4);
latex(sym(A2,'d'))
digits(d1);
close;