close all; clc; clear;

file = 'grsfine';

fid = fopen([file,'.dat']);
C = textscan(fid, '%s','delimiter', '\n');
RL=C{1}{1};
parsed = strsplit(RL,'|');
NT = strsplit(parsed{1},' = ');
Ntheta = str2num(NT{2});
NP = strsplit(parsed{2},' = ');
Nphi = str2num(NP{2});
clear NT NP parsed RL 

fs = 18; FS = 'FontSize';
fw = 'bold'; FW = 'FontWeight';

f = NaN([Ntheta,Nphi]);
costheta = linspace(1,-1,Ntheta);
phi = linspace(0,2*pi,Nphi)-pi;
theta = acos(costheta)-pi/2;

ii = 1;
for j = 1:Nphi
    for i = 1:Ntheta
        XYZC = strsplit(C{1}{2+ii},' ');
        x(ii) = str2double(XYZC{1});
        y(ii) = str2double(XYZC{2});
        z(ii) = str2double(XYZC{3});
        c(ii) = str2double(XYZC{4});
        f(i,j) = c(ii);
        ii = ii + 1;
    end
end

c_max = max(abs(c));
x = x'; y = y'; z = z'; c = c';
[T,P] = meshgrid(theta,phi);
[HX,HY] = sph2hammer(P,T);

fig = figure(1);
pos = [0 0 600 1200];
set(gcf,'Position',pos)

% subplot(2,1,1)
pcolor(fliplr(HY)*180/pi+90,HX*180/pi,f');
hold on;
xxx = [-180:180]; yyy = ones(size(xxx));
% plot((73)*yyy,xxx,'--')
% xlim([-pi*1.02, pi*1.02]);
% ylim([-pi/2*1.02,pi/2*1.02]);
% title('Maximum value at \theta = 73.7 deg',FS,fs,FW,fw)
ylabel('\phi (deg)',FS,fs,FW,fw)
yticks([-180 -135 -90 -45 0 45 90 135 180])
yticklabels({'-180','-135','-90','-45','0','45','90','135','180'})
xlabel('\theta (deg)',FS,fs,FW,fw)
xticks([0 15 30 45 60 75 90 105 120 135 150 165 180])
xticklabels({'0','15','30','45','60','75','90','105','120','135','150','165','180'})
colorbar, shading interp, daspect([1 1 1]);
% 
% subplot(2,1,2)
% title('Torque efficiency on sphere',FS,fs,FW,fw)
% tri = delaunaySph([x,y,z]);
% trisurf(tri,x,y,z,c);
% % surf(X,Y,Z,f);
% hold on;
% % plot3(x(i),y(i),z(i),'r.')
% [i,~] = find(c == max(c));
% colorbar; shading interp; daspect([1 1 1]); axis tight;
% view(45,5)
% xlabel('x',FS,fs,FW,fw); ylabel('y',FS,fs,FW,fw); zlabel('z',FS,fs,FW,fw);
% xticks([-1 -0.5 0 0.5 1]); yticks([-1 -0.5 0 0.5 1]); zticks([-1 -0.5 0 0.5 1]);
% grid off
% axis equal

print(file,'-dpng')