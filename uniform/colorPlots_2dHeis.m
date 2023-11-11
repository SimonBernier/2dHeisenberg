clearvars

Ly = 4; Lx = 45; h=6.0; tau = 2; maxDim = 800;
dt = 0.1; tanhshift = 4; gse = 1;

R = Ly/2; center = (Lx+1)/2 + (-R:R);

cd data_1E-8\
filename = sprintf("Ly_%d_Lx_%d_h_%0.2f_tau_%0.1f_maxDim_%d_gse_%d_2dHeis_uni.dat",Ly,Lx,h,tau,maxDim,gse);
data = importdata(filename,' ',1);
cd ..\

[tval, en, enf, enf_en0, svn, localEn, localEn_En0, corrZ, corrPerp] = collectData(data, Lx);
Nt = length(tval); x = (1:Lx-1) - Lx/2;

tval = 0:dt:dt*(length(tval)-1);

minHeatMap = min(0, min(localEn_En0(:)));
maxHeatMap = max(localEn_En0(:));

f1 = figure(1); f1.set('Position',[50 100 1350 400]); clf, box on
s1 = subplot(1,3,1,'Position',[0.05 0.2 0.28 0.7]);
ylabel('time')
ylim([tval(1) tval(end)]), xlim([x(1) x(end)])
set(gca,'FontName','Times','FontSize',15)

imagesc('XData',x,'YData',tval,'CData', localEn_En0)
colormap('jet'), axis('tight'), clim([minHeatMap maxHeatMap]), colorbar

s2 = subplot(1,3,2,'Position',[0.38 0.2 0.28 0.7]);
xlabel('x-x0')
ylim([tval(1) tval(end)]), xlim([x(1) x(end)])
set(gca,'FontName','Times','FontSize',15)

imagesc('XData',x,'YData',tval,'CData', svn)
colormap('jet'), axis('tight'), colorbar

s3 = subplot(1,3,3,'Position',[0.71 0.2 0.28 0.7]);
ylim([tval(1) tval(end)]), xlim([x(1) x(end)])
set(gca,'FontName','Times','FontSize',15)

x = (1:Lx) - (Lx+1)/2;
imagesc('XData',x,'YData',tval,'CData', corrPerp)
colormap('jet'), axis('tight'), colorbar
clim( [min(corrPerp,[],'all'), max(corrPerp(:,[1:(Lx-1)/2 (Lx+1)/2+1:end]),[],'all') ] )

enCenter = mean(localEn_En0(:,center),2);

enAsym = localEn_En0(:,(Lx-1)/2:-1:1) - localEn_En0(:, (Lx-1)/2+1:Lx-1);
enAsymTot = sum(abs(enAsym),2);

svnAsym = svn(:,(Lx-1)/2:-1:1) - svn(:, (Lx-1)/2+1:Lx-1);
svnAsymTot = sum(abs(svnAsym),2);

corrZAsym = corrZ(:,(Lx-1)/2:-1:1) - corrZ(:, (Lx+1)/2+1:Lx);
corrZAsymTot = sum(abs(corrZAsym),2);

corrPAsym = corrPerp(:,(Lx-1)/2:-1:1) - corrPerp(:, (Lx+1)/2+1:Lx);
corrPAsymTot = sum(abs(corrPAsym),2);

figure(2), clf, box on
hold on
plot(tval, enAsymTot, 'DisplayName','en')
plot(tval, svnAsymTot, 'DisplayName','svn')
plot(tval, corrZAsymTot, 'DisplayName','corrZ')
plot(tval, corrPAsymTot, 'DisplayName','corrP')
hold off
ylabel('\Delta L/R'), xlabel('Time'), legend('Location','best')
set(gca,'FontName','Times','FontSize',15)

%% function to get data
function [tval, en, enf, enf_en0, svn, localEn0, localEn, corrZ, corrPerp] = collectData(A,Lx)
    tval = A.data(1:end,1);
    en = A.data(1:end,2);
    enf = A.data(1:end,3);
    enf_en0 = A.data(1:end,4);
    svn = A.data(1:end,5:4+(Lx-1));
    localEn0 = A.data(1:end, 5 + (Lx-1) : 4 + 2*(Lx-1));
    localEn = A.data(1:end,  5 + 2*(Lx-1) : 4 + 3*(Lx-1));
    corrZ = A.data(1:end,    5 + 3*(Lx-1) : 4 + 3*(Lx-1) + Lx);
    corrPerp = A.data(1:end, 5 + 3*(Lx-1) + Lx: end);
end