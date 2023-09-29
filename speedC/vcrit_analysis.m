clearvars

Ly = 2; Lx = 8*Ly; h = 0;
N = Lx*Ly;

loc = Lx/2;

filename = sprintf('Ly_%d_Lx_%d_h_%0.2f_c2dHeis.dat',Ly,Lx,h);

%% get 2nd order method
data = importdata(filename,' ',1);
[tval, enPsi, maxBondDim, enPhi, svn] = collectData(data);

deltaSVN = svn; deltaSVN(1,:) = zeros(size(deltaSVN(1,:)));

%% plot energy and bond dimension
h1 = figure(1); h1.set('Position',[50 50 800 400]), clf, box on
subplot(1,2,1, 'Position', [0.1 0.15 0.4 0.8]);
plot(tval, enPhi, 'LineWidth', 2.0)
xlabel('time t'), ylabel('Energy density')
set(gca,'FontSize',15,'FontName','Times')

subplot(1,2,2, 'Position',[0.58 0.15 0.4 0.8]);
plot(tval, maxBondDim, 'LineWidth', 2.0)
xlabel('time Jt'), ylabel('Max bond dimension')
set(gca,'FontSize',15,'FontName','Times')

%% plot svn
h2 = figure(2); h2.set('Position',[50 50 1000 400]), clf, box on

x = (1:Lx-1) - loc + 0.5;
s1 = subplot(1,2,1, 'Position', [0.07 0.14 0.4 0.85]);
imagesc('XData',x,'YData',tval,'CData',abs(deltaSVN))
xlabel('x'), ylabel('time Jt'), colormap('jet'), axis('tight'), colorbar
xlim([min(x) max(x)]), ylim([0 max(tval)])
set(gca,'FontSize',15,'FontName','Times')

dsvnFFT = fft2(deltaSVN);
k = (x-1)/Lx; omega = tval(1:floor(length(tval))/2)/tval(end);
subplot(1,2,2, 'Position', [0.55 0.14 0.4 0.85]), title('\Delta SvN(k,\omega)')
imagesc('XData',k,'YData',omega,'CData',abs(dsvnFFT(1:floor(length(tval))/2,:)))
xlabel('k/pi'), ylabel('\omega/pi'), colormap('jet'), axis('tight'), colorbar
set(gca,'FontSize',15,'FontName','Times')

%% calculate velocity
h3 = figure(3); clf, box on
hold on
for i=1:10
    r = 1 + floor(max(tval))*i;
    plot((x - pi/2*tval(r))/tval(r), deltaSVN(r,:), '-', 'LineWidth',1, 'DisplayName',sprintf('t = %0.1f', tval(r)));
end
hold off
xlabel('x'), ylabel('\Delta SvN(x,t)')
legend('Location','best')
xlim([min(x) max(x)])
set(gca,'FontSize',15,'FontName','Times')

cutoff = 0.01;

T = zeros(size(x));
for j = 1:Lx-1
    temp = find( deltaSVN(:,j) < cutoff, 1, 'last');
    T(j) = tval(temp);
end
figure(h2),subplot(s1)
hold on
plot(x(2:end-1), T(2:end-1), '-', 'LineWidth',2, 'Color','white')
hold off

k = find( T < max(tval), 1, 'last');
[ft, gof] = fit( T(loc:k)', x(loc:k)', 'poly1');
v = ft.p1; vErr = gof.rmse/sqrt(gof.dfe+2);
fprintf("v = %0.4f +/- %0.4f\n", v, round(vErr,1,'significant'))

k = find( T < max(tval), 1, 'first');
[ft, gof] = fit( T(2:loc-1)', x(2:loc-1)', 'poly1');
v = -ft.p1; vErr = gof.rmse/sqrt(gof.dfe+2);
fprintf("v = %0.4f +/- %0.4f\n", v, round(vErr,1,'significant'))


%% function to get data
function [tval, en_psi, maxBondDim, en_phi, svn] = collectData(A)
    tval = A.data(:,1);
    en_psi = A.data(:,2);
    maxBondDim = A.data(:,3);
    en_phi = A.data(:,4);
    svn = A.data(:,5:end);
end
