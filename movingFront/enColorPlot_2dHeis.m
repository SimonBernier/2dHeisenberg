clearvars

Ly = 2; Lx = 8*Ly; h=5.0; v = 2; tau = 0.2; maxDim = 512;
dt = 0.1; tanhshift = 3;

R = Ly; center = Lx/2-R:Lx/2+R; N = Lx*(Ly-1);

% cd data_1E-8\
filename = sprintf("Ly_%d_Lx_%d_h_%0.2f_v_%0.2f_tau_%0.1f_maxDim_%d_2dHeis_mf.dat",Ly,Lx,h,v,tau,maxDim);
data = importdata(filename,' ',1);
% cd ..\

[tval, en, enf, enf_en0, svn, localEn, localEn_En0, corrZ, corrPerp] = collectData(data, Lx, Ly);
Nt = length(tval); x = (1:Lx-1)-Lx/2; y = 1:Ly;

minHeatMap = min(0, min(localEn_En0(:)));
maxHeatMap = max(localEn_En0(:));

f1 = figure(1); f1.set('Position',[100 100 1200 400]); clf, box on
s1 = subplot(1,1,1,'Position',[0.05 0.2 0.92 0.7]);
xlabel('x-x0'), ylabel('time')
ylim([tval(1) tval(end)]), xlim([x(1) x(end)])
set(gca,'FontName','Times','FontSize',15)

imagesc('XData',x,'YData',tval,'CData', localEn_En0)
colormap('jet'), axis('tight'), clim([minHeatMap maxHeatMap]), colorbar

enCenter = mean(localEn_En0(:,center),2);

%% function to get data
function [tval, en, enf, enf_en0, svn, localEn0, localEn, corrZ, corrPerp] = collectData(A,Nx,Ny)
    tval = A.data(1:end,1);
    en = A.data(1:end,2);
    enf = A.data(1:end,3);
    enf_en0 = A.data(1:end,4);
    svn = A.data(1:end,5);
    N = (Nx-1)*Ny;
    localEn0 = A.data(1:end,   6 : Nx+4);
    localEn = A.data(1:end,  Nx+5 : 2*Nx+3);
    corrZ = A.data(1:end,  2*Nx+4 : 2*Nx+N+5);
    corrPerp = A.data(1:end,   2*Nx+N+6:end);
end