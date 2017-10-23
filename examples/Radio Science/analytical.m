close all; clc; clear;

fs = 18; FS = 'fontsize';
fw = 'bold'; FW = 'fontweight';
lw1 = 3; lw2 = 5; LW = 'linewidth';
pos = [0 0 1000 1000];
set(gcf,'Position',pos)

files = {'log','log2','log3'};
fmax = 3;
ns = [2000,20000,200000];
iskip = [1,10,100];

diff1 = zeros(ns(1),1);
diff2 = zeros(ns(2),1);
diff3 = zeros(ns(3),1);

w1 = zeros(3, ns(1));
w2 = zeros(3, ns(2));
w3 = zeros(3, ns(3));
for f = 1:fmax
    file = files{f};
    n = ns(f);

    fid = fopen(file);
    C = textscan(fid, '%s','delimiter', '\n');

    w0 = str2double(strsplit(C{1}{5},' '));
    w0 = w0(3:5);

    Ip = str2double(strsplit(C{1}{16},' '));
    Omega = (Ip(3)-Ip(1))/Ip(1)*w0(3);
    
    Q1 = str2double(strsplit(C{1}{18},' '));
    Q2 = str2double(strsplit(C{1}{19},' '));
    Q3 = str2double(strsplit(C{1}{20},' '));
    Q = [Q1;Q2;Q3];

    format long

    tt1 = strsplit(C{1}{23},'|');
    tt1 = str2double(tt1{8});
    tt2 = strsplit(C{1}{24},'|');
    tt2 = str2double(tt2{8});
    dt = double(tt2-tt1);

    for i = 1:n
        line = strsplit(C{1}{22+i},'|');
        Nt(i) = str2num(line{1});
        wstr = str2double(strsplit(line{4},' '));
        Rstr = str2double(strsplit(line{9},' '));
        tt = str2double(line{8});
        R = reshape(Rstr(2:10),[3,3]);
        ww = R'*Q'*wstr(2:4)';
        wa = [w0(1)*cos(Omega*tt);w0(1)*sin(Omega*tt);w0(3)];
        if(f==1)
            diff1(i) = wa(1)-ww(1);
            w1(:,i) = ww;
            t1(i) = tt;
        elseif (f==2)
            diff2(i) = wa(1)-ww(1);
            w2(:,i) = ww;
            t2(i) = tt;
        else
            diff3(i) = wa(1)-ww(1);
            w3(:,i) = ww;
            t3(i) = tt;
        end
        wx_a = w0(1)*cos(Omega*tt);
    end
    fclose(fid);
    if(f==1)
        fig = figure(1);
        t = [1:n]*2*pi/(n*Omega);
        wx_a = w0(1)*cos(Omega*t);
        wy_a = w0(1)*sin(Omega*t);
        wz_a = w0(3)*ones(size(t));

        wx = squeeze(w1(1,:)); wy = squeeze(w1(2,:)); wz = squeeze(w1(3,:));
        plot(wx,wy,'--',LW,lw2);
        hold on;
        axis tight
        axis equal
        xlabel('\omega_{1,b}',FW,fw,FS,fs)
        ylabel('\omega_{2,b}',FW,fw,FS,fs)
        plot(wx_a,wy_a,LW,lw1)
        h = legend('Numerical result','Analytical result');
        h.FontSize = fs;
        set(gca,FS,fs) 
        print('analytical1','-dpng')
    end
end

Omega = (Ip(3)-Ip(1))/Ip(1)*w0(3);
fig = figure(2);
pos = [0 0 1000 1000];
set(gcf,'Position',pos)
for f = 1:fmax
    if(f==1)
        semilogy(t1(:),abs(diff1(:)))
    elseif (f==2)
        semilogy(t2(:),abs(diff2(:)))
    else
        semilogy(t3(:),abs(diff3(:)))
    end
    hold on
end
xlabel('t (s)',FW,fw)
ylabel('|error of \omega_{1,b}|',FW,fw)
xlim([t1(1),t1(end)])
xticks([t1(1), t1(ns(1)/2), t1(end)])
xticklabels({floor(num2str(t1(1)-1)),floor(num2str(t1(ns(1)/2))),floor(num2str(t1(end)))})
legend('\theta_{max} = 0.1','\theta_{max} = 0.01','\theta_{max} = 0.001')
% ylim([-1,1]*0.01)
set(gca,FS,fs)
print('analytical2','-dpng')