E_r = hdf5read('E_field.h5','A_r');
Es_r = hdf5read('E_scat.h5','A_r');
grid = hdf5read('grid_xyz.h5','A_r');

Ex = E_r(1,:);
Esx = Es_r(1,:);
n = sqrt(size(Ex,2));
Ex = reshape(Ex,n,n);
Esx = reshape(Esx,n,n);
Y = reshape(grid(2,:),n,n);
Z = reshape(grid(3,:),n,n);
figure(1)
subplot(2,1,1)
surf(Z,Y,Ex)
subplot(2,1,2)
surf(Z,Y,Esx)
