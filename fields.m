fs = 18; FS = 'FontSize';
fw = 'bold'; FW = 'FontWeight';

E_r = hdf5read('E_field.h5','A_r');
Es_r = hdf5read('E_scat.h5','A_r');
grid = hdf5read('grid_xyz.h5','A_r');

a = 2e-7;

Ex = E_r(1,:);
Esx = Es_r(1,:);
n = sqrt(size(Ex,2));
Ex = reshape(Ex,n,n);
Esx = reshape(Esx,n,n);
norm = max(max(abs(Ex)));
Esx(abs(Esx)>3)=NaN;
Ex = Ex/norm;

Y = reshape(grid(2,:),n,n)/a;
Z = reshape(grid(3,:),n,n)/a;
fig = figure(1);
pos = [0 0 500 1200];
set(gcf,'Position',pos)

subplot(2,1,1)
h = pcolor(Z,Y,Ex);
shading interp;
set(h,'EdgeColor','none');
ylabel('ky',FS,fs,FW,fw)
xlabel('kz',FS,fs,FW,fw)
title('Incident field',FS,fs,FW,fw)
xt = get(gca, 'XTick');
set(gca, FS,fs)

subplot(2,1,2);
h2 = pcolor(Z,Y,Esx);
shading interp;
set(h2,'EdgeColor','none');
ylabel('ky',FS,fs,FW,fw)
xlabel('kz',FS,fs,FW,fw)
title('Scattered field',FS,fs,FW,fw)
xt = get(gca, 'XTick');
set(gca, FS,fs)
caxis([-1,1])
