clear; clc; close all;

fig = figure;
figname = ['mueller.png'];

% dat = importdata('mueller-perf_ori');
% 
% S_places = [3,4,5,6];
% N_phi = 180;
% N_theta = 90;
% 
% phi = dat.data(:,2);
% theta = dat.data(:,1);

points = importdata('points');
dat = importdata('mueller');

S_places = [3,4,5,6]+1;
N_phi = 18;
N_theta = 9;

range = [1:N_theta*N_phi];

theta = points.data(range,1);
phi = points.data(range,2);

range = range+0*N_theta*N_phi;

S11 = dat.data(range,S_places(1)); 
S12 = dat.data(range,S_places(2));
S13 = dat.data(range,S_places(3)); 
S14 = dat.data(range,S_places(4));

datum = [S11/S11(1), -S12./S11, S13./S11, S14./S11];
s13mx = max(abs(max(datum(:,3))),abs(min(datum(:,3))));
s14mx = max(abs(max(datum(:,4))),abs(min(datum(:,4))));
s112d = datum(:,1); s112d = reshape(s112d,[N_phi,N_theta]);
s112d = s112d(:,1:N_theta);

s122d = datum(:,2); s122d = reshape(s122d,[N_phi,N_theta]);
s122d = s122d(:,1:N_theta);
s132d = datum(:,3); s132d = reshape(s132d,[N_phi,N_theta]);
s132d = s132d(:,1:N_theta);
s142d = datum(:,4); s142d = reshape(s142d,[N_phi,N_theta]);
s142d = s142d(:,1:N_theta);
datum2 = zeros(N_phi,N_theta,3);
datum2(:,:,1) = s122d; datum2(:,:,2) = s132d;
datum2(:,:,3) = s142d;

[X,Y] = meshgrid(0:360/N_theta:359,0:180/N_phi:179);

s12lim = [-1.1,1.1];
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
        view(30,30)
    else
        ax2 = subplot(2,2,i);
        s = surf(X,Y,datum2(:,:,i-1),'FaceLighting','gouraud','FaceColor','interp');
        set(gca,'ylim',[0 180],'xlim',[0 360],'zlim',yl(i-1,:))
        light
        colormap(ax2,parula)
%                     set(gca,'Ydir','reverse')
        view(30,30)
    end
    hold on;
    material dull
    s.EdgeColor = 'none';
    p = get(gca,'position');
    p(4) = p(4)*1.2; 
    p(3) = p(3)*1.08; 
    set(gca, 'position', p);

end

print(figname,'-dpng')
