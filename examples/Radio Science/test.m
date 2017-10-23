close all; clc; clear;

fs = 14; FS = 'fontsize';
fw = 'bold'; FW = 'fontweight';
lw1 = 3; lw2 = 5; LW = 'linewidth';
% pos = [0 0 1000 1000];
% set(gcf,'Position',pos)

file = 'Ngrs.out';

fid = fopen(file);
C = textscan(fid, '%s','delimiter', '\n');

Ntheta = 180;
for i = 1:Ntheta
    XYZC = strsplit(C{1}{1+i},' ');
    costheta(i) = str2double(XYZC{1});
    Nk(i) = str2double(XYZC{2});
    N0(i) = str2double(XYZC{3});
    N90(i) = str2double(XYZC{4});
    theta(i) = acos(costheta(i));
end
% costheta = flip(costheta); theta = flip(theta); Nk = flip(Nk);
% N0 = flip(N0); N90 = flip(N90);

N1 = Nk/abs(max(Nk));
N2 = N0/abs(max(N0));
N3 = N90/abs(max(N90));
% plot(costheta,N)
Ly = 0:10:1000;

[GX,GY] = meshgrid(costheta,Ly);
U = zeros(size(GX,1),size(GX,2));
V = zeros(size(GY,1),size(GY,2));
for i = 1:length(Ly)
   for j = 1:Ntheta
      Li = Ly(i)*[cos(theta(j)), sin(theta(j)),0];
      Lf = Li + [N1(j),N2(j),N3(j)];
      if(norm(Lf)*norm(Li)~=0)
         U(i,j) = Lf(1)/norm(Lf)-Li(1)/norm(Li);
      else
         U(i,j) = 0;
      end
      V(i,j) = norm(Lf)-norm(Li);
   end
end
U = U/100; V = V/100;
sx = costheta(1:10:Ntheta);
n = 3;
sy = Ly(1+n)*ones(1,2*size(sx,2)); sy(size(sx,2)+1:end) = Ly(floor(length(Ly)/2));
sx = [sx,sx];

% quiver(GX,GY,U,V)
streamline(GX,GY,U,V,sx,sy)

% figure
% plot(costheta,N1,costheta,N2,costheta,N3)
% legend('n1','n2','n3')

% w0 = str2double(strsplit(C{1}{5},' '));
