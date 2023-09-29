clearvars
Ly = 2:5; Lx = 8*Ly; N = Ly.*Lx;
h = 0:0.1:9;

en0 = zeros(length(Ly),length(h)); var0 = en0; maxBondDim0 = en0;
gap = en0; var1 = en0; maxBondDim1 = en0;

cd data1E-8\
for i=1:length(Ly)

    filename = sprintf('Ly_%d_Lx_%d_2dHeis_gap_checker.dat',Ly(i),Lx(i));
    A = importdata(filename,' ',1);
    [~, en0(i,:), var0(i,:), maxBondDim0(i,:),...
     gap(i,:), var1(i,:), maxBondDim1(i,:)] = accessData(A);   

end
cd ..\

%% Fit the minimum gap
% minGap = gap(:,1);
% gapError = sqrt(var0)+sqrt(var1);
% ftz = fit(Ly, minGap, 'poly1');
% z = 1;
% 
% f3 = figure(3);
% p = plot(ftz, Ly, minGap);
% p.set('MarkerSize',15,'LineWidth',1.5);
% xlabel('log N'), ylabel('log \Delta')
% set(gca, 'FontName','Times','FontSize',15)

nu = 2; z = 1;

f1 = figure(1); f1.set('Position',[50 50 1000 400])
s1 = subplot(1,2,1, 'Position',[0.1 0.15 0.4 0.8]); box on
s2 = subplot(1,2,2, 'Position',[0.58 0.15 0.4 0.8]);
for i=1:length(Ly)

    subplot(s1)
    hold on
    p = plot(h,gap(i,:), 'DisplayName',sprintf('Ly=%d',Ly(i)),'LineWidth',1.0);
    hold off

end
subplot(s1)
xlabel('h'), ylabel('Gap \Delta')
legend('Location','best'), xlim([0 9])
set(gca, 'FontName','Times','FontSize',15)

subplot(s2), box on
minGap = gap(:,1); gapError = (sqrt(var0)+sqrt(var1))./sqrt(N');
% errorbar(Ly, minGap, gapError(:,1),'.-','MarkerSize',10, 'Color', p.Color)
plot(Ly, minGap.*(Ly').^2, '.-','MarkerSize',10, 'Color', p.Color)
xlabel('Ly'), ylabel('\Delta(Ly,h=0)')
xlim([1 6]), set(gca, 'FontName','Times','FontSize',15)

% exportgraphics(f1, sprintf('scaledGap_Lx_%d.png', Lx), 'Resolution',300);

function [h, en0, var0, maxBondDim0, gap, var1, maxBondDim1] = accessData(S)
    h = S.data(:,1);
    en0 = S.data(:,2); var0 = S.data(:,3); maxBondDim0 = S.data(:,4);
    gap = S.data(:,5); var1 = S.data(:,6); maxBondDim1 = S.data(:,7);
end
