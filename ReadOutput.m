fileID = fopen('data.bin');
N = fread(fileID,1,'int');
x = (fread(fileID,N+1,'double'))';
M = fread(fileID,1,'int');
y = (fread(fileID,M+1,'double'))';
S = (fread(fileID,[N+1,M+1],'double'))';
Sconc = (fread(fileID,[N+1,M+1],'double'))';
Sinterface = (fread(fileID,[N+1,M+1],'double'))';

ChanWidth = (fread(fileID,[N+1,M+1],'double'))';
P = (fread(fileID,[N+1,M+1],'double'))';
U = (fread(fileID,[N,M+1],'double'))';
V = (fread(fileID,[N+1,M],'double'))';
Fluxout = (fread(fileID,[N+1,M+1],'double'))';
Growth = (fread(fileID,[N+1,M+1],'double'))';

SconcfromMaxandMin = (2*S+Sinterface)/3;
ICadjust = Sconc(1,1);
Sconcfromnewton = ICadjust*sqrt(S-ICadjust)./sqrt(S-Sinterface) + 1/3*(Sinterface + 2*S - (ICadjust+2*S).*sqrt(S-ICadjust)./sqrt(S-Sinterface));

fclose(fileID);

figure;
imagesc(x,y,(.5-max(ChanWidth,0))/2);
h = colorbar;
ylabel(h, 'Height [mm]')
% title({'Biofilm Height with Outward Growth',... 
%     'after 4 Days'},'fontweight','Normal');
title({'Merged Biofilm after 4 Days'},'fontweight','Normal');

xlabel('x [mm]');
ylabel('y [mm]');
set(gca,'YDir','normal')
ax = gca;
ax.FontSize = 20;


figure
[X,Y] = meshgrid(x,y);
mesh(X,Y,(.5-max(ChanWidth,0))/2);
zlim([0,1])
%title({'Surface Plot of Biofilm with Outward', ...
%    'Growth after 4 Days'},'fontweight','Normal');
%title('Surface Plot of Initial Biofilm','fontweight','Normal');
title({'Initial Separated Biofilm'},'fontweight','Normal');

xlabel('x [mm]');
ylabel('y [mm]');
zlabel('height [mm]');
ax = gca;
ax.FontSize = 20;
