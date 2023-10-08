clearvars

Ly = 4; Lx = 9; h = 0; maxDim = 800; LyPlot = 5;
gse = 1 ;

v = zeros(size(Ly)); vErr = v; deltaE = v;

for i=1:length(Ly)
    fprintf("Ly = %d\n", Ly(i))
    loc = Lx(i)/2;

    filename = sprintf('Ly_%d_Lx_%d_h_%0.2f_maxDim_%d_gse_%d_c2dHeis.dat',Ly(i),Lx(i),h,maxDim,gse);
    % filename = sprintf('Ly_%d_Lx_%d_h_%0.2f_maxDim_%d_c2dHeis.dat',Ly(i),Lx(i),h,maxDim);

    %% get 2nd order method
    data = importdata(filename,' ',1);
    [tval, enPsi, maxBondDim, enPhi, svn] = collectData(data);
    dt = tval(2) - tval(1);
    
    deltaSVN = svn; deltaSVN(1,:) = zeros(size(deltaSVN(1,:)));
    x = (1:Lx(i)-1) - loc + 0.5;

    %% plot energy and bond dimension
    % h1 = figure(1); h1.set('Position',[50 50 800 400]), clf, box on
    % subplot(1,2,1, 'Position', [0.1 0.15 0.4 0.8]);
    % plot(tval, enPhi, 'LineWidth', 2.0)
    % xlabel('time t'), ylabel('Energy density')
    % set(gca,'FontSize',15,'FontName','Times')
    % 
    % subplot(1,2,2, 'Position',[0.58 0.15 0.4 0.8]);
    % plot(tval, maxBondDim, 'LineWidth', 2.0)
    % xlabel('time Jt'), ylabel('Max bond dimension')
    % set(gca,'FontSize',15,'FontName','Times')

    %% plot svn
    if Ly(i) == LyPlot
        h2 = figure(2); h2.set('Position',[50 50 1000 400]), clf, box on
        
        s1 = subplot(1,2,1, 'Position', [0.07 0.14 0.4 0.85]);
        imagesc('XData',x,'YData',tval,'CData',abs(deltaSVN))
        xlabel('x'), ylabel('time Jt'), colormap('jet'), axis('tight'), colorbar
        xlim([min(x) max(x)]), ylim([0 max(tval)])
        set(gca,'FontSize',15,'FontName','Times')
    
        dsvnFFT = fft2(deltaSVN);
        k = (x-1)/Lx(i); omega = tval(1:floor(length(tval))/2)/tval(end);
        s2 = subplot(1,2,2, 'Position', [0.55 0.14 0.4 0.85]); title('\Delta SvN(k,\omega)')
        imagesc('XData',k,'YData',omega,'CData',abs(dsvnFFT(1:floor(length(tval))/2,:)))
        xlabel('k/pi'), ylabel('\omega/pi'), colormap('jet'), axis('tight'), colorbar
        set(gca,'FontSize',15,'FontName','Times')
    end

    %% calculate velocity
    % h3 = figure(3); clf, box on
    % hold on
    % for i=1:10
    %     r = 1 + floor(max(tval))*i;
    %     plot(x/tval(r), deltaSVN(r,:), '-', 'LineWidth',1, 'DisplayName',sprintf('t = %0.1f', tval(r)));
    % end
    % hold off
    % xlabel('x/t'), ylabel('\Delta SvN(x)')
    % legend('Location','best')
    % xlim([min(x) max(x)])
    % set(gca,'FontSize',15,'FontName','Times')
    
    cutoff = 0.01;

    T = zeros(size(x));
    for j = 1:Lx(i)-1
        temp = find( deltaSVN(:,j) < cutoff, 1, 'last');
        T(j) = tval(temp);
    end
    if Ly(i) == LyPlot
        figure(h2),subplot(s1)
        hold on
        plot(x(1:end), T(1:end), '-', 'LineWidth',2, 'Color','white')
        hold off
    end

    k = min( find( T < max(tval), 1, 'last'), length(T)-1);
    [ft, gof] = fit( T(loc+1:k)', x(loc+1:k)', 'poly1');
    vR = ft.p1; vErrR = gof.rmse/sqrt(gof.dfe+2);
    fprintf("\tvR = %0.4f +/- %0.4f\n", vR, round(vErrR,1,'significant'))

    if Ly(i) == LyPlot
        h4 = figure(4); clf, box on
        plot(ft, T(loc+1:k), x(loc+1:k),'s')
    end

    k = max( find( T < max(tval), 1, 'first'), 2);
    [ft, gof] = fit( T(k:loc-2)', x(k:loc-2)', 'poly1');
    vL = -ft.p1; vErrL = gof.rmse/sqrt(gof.dfe+2);
    fprintf("\tvL = %0.4f +/- %0.4f\n", vL, round(vErrL,1,'significant'))

    if Ly(i) == LyPlot
        hold on, plot(ft, T(k:loc-2), x(k:loc-2), 'o'), hold off
        xlabel('Time'), ylabel('x')
        legend([{'x>0'},{'fit'},{'x<0'},{'fit'}],'Location','best')
        set(gca,'LineWidth',1,'FontName','Times','FontSize',15)
    end

    fprintf("\tv = %0.4f +/- %0.4f\n", (vR+vL)/2, round(sqrt( (vErrR)^2+(vErrL)^2 )/2 ,1,'significant'))
    fprintf("\tEnergy Diff = %0.4f%% \n", (enPhi( max(T)/dt + 1)-enPhi(1))/enPhi(1)*100)

    v(i) = (vR+vL)/2; vErr(i) = sqrt( (vErrR)^2+(vErrL)^2 )/2;
    deltaE(i) = (enPhi( max(T)/dt+1)-enPhi(1))/enPhi(1)*100;

end
figure(5);
hold on
p = plot(Ly,v,'.','MarkerSize',15);
errorbar(Ly,v,vErr,'Color',p.Color)
plot([2 5],pi/2*[1 1], '--r'), hold off
xlabel('Ly'), ylabel('v')
xlim([1.9 5.1])
set(gca,'LineWidth',1,'FontName','Times','FontSize',15)

assymetry = svn(:,(Lx-1)/2:-1:1) - svn(:, (Lx-1)/2+1:Lx-1);
figure(6)
imagesc('XData',x((Lx-1)/2+1:Lx-1), 'YData',tval,'CData',abs(assymetry))
colorbar, xlabel('x'), ylabel('t')

figure(7)
hold on
totalAssym = sum(abs(assymetry),2);
plot(tval, totalAssym, 'DisplayName', sprintf("GSE %d, dt = %0.3f", gse, dt))
hold off
legend()

%% function to get data
function [tval, en_psi, maxBondDim, en_phi, svn] = collectData(A)
    tval = A.data(:,1);
    en_psi = A.data(:,2);
    maxBondDim = A.data(:,3);
    en_phi = A.data(:,4);
    svn = A.data(:,5:end);
end
