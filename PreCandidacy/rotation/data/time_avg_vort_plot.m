close all; clear; clc;

format shortE

nd = 766^2;

% this reads all data tim from multiple timestamps into 1 array (it ignores the lineskips in the .dat file)
vort = readmatrix("Om1B30.dat");
invRo = 1;
invFr = sqrt(30);

s1 = size(vort)
sd = s1(1)/nd

% in order to plot this by individual timestep, you have to restrict the plot by column 3 which is the timestep. 
num_timesteps = 766^2;

time_col = 3;
[row,col] = find(vort(:,time_col) == vort(s1(1),time_col));
s2 = size(vort(row,time_col));

np = 64;
colors = zeros(np);

xcol_num = 7;
ycol_num = 8;

min1 = min(vort(:,xcol_num));
max1 = max(vort(:,xcol_num));
min2 = min(vort(:,ycol_num));
max2 = max(vort(:,ycol_num));
d1 = (max1 - min1)/(np-1);
d2 = (max2 - min2)/(np-1);

for i = 1:s1(1)
    x = floor((vort(i,xcol_num)-min1)/d1)+1;
    y = np-floor((vort(i,ycol_num)-min2)/d2);
    % the plotting with imagesc is weird (x-y are switched)
    colors(y,x) = colors(y,x) + 1;
end
colors = colors/sd;

figure('color','white')

plot(vort(:,xcol_num)+invRo,vort(:,ycol_num), 'o')
xlim([min1, max1]+2);

figure('color','white')

imagesc(log(colors))
xtickformat('%.5f');
ytickformat('%.5f');
num_ticks = 10;
xticks(linspace(0,64,num_ticks+1));
yticks(linspace(0,64,num_ticks+1));
xticklabels(round(linspace(min1,max1,num_ticks+1)+invRo,2));
yticklabels(round(flip(linspace(min2,max2,num_ticks+1)),2));
colormap sky
xlabel("$\overline{\omega_z} + 2\Omega$", 'Interpreter','latex','Fontsize',22)
ylabel("$\overline{|\nabla T|^2}$", 'Interpreter','latex', 'Rotation', 0, 'FontSize',22)

figure('color','white')

x1 = 1:766;
y1 = 1:766;

z1 = reshape(vort(row,7),[766, 766]);
z2 = abs(invFr./(z1+invRo));
z3 = reshape(vort(row,8),[766, 766])/2;
z4 = reshape(vort(row,9),[766, 766]);
z5 = reshape(vort(row,12),[766, 766])/(4*600);


surf(x1, y1, z1+invRo)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap(redwhiteblue(min(z1+invRo,[],'all'),max(z1+invRo,[],'all')))
colorbar
title("$\widehat{\omega_z} + 2\Omega$", 'Interpreter','latex','FontSize',22)
%clim([-35,35])
saveas(gcf,'Om1B30_vortz_bar.pdf')

figure('color','white')

surf(x1, y1, z2)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
clim([0 1])
colorbar
title("$\frac{Ro_L}{Fr}$",'Interpreter','latex','FontSize',22)
saveas(gcf,'Om1B30_RC.pdf')

figure('color','white')

surf(x1, y1, z3)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
colorbar
%clim([0,3.5])
title("$\frac{\widehat{|\nabla T|^2}}{Fr^2Pe}$",'Interpreter','latex','FontSize',22)
saveas(gcf,'Om1B30_tdisp_bar.pdf')

figure('color','white')

surf(x1, y1, z4)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
colorbar
%clim([0,1.2])
title("$\widehat{w^2}$",'Interpreter','latex','FontSize',22)
saveas(gcf,'Om1B30_uzrms_bar.pdf')

figure('color','white')

surf(x1, y1, z5)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
colorbar
%clim([0,1.2])
title("$\frac{\widehat{|\nabla \mathbf{u}|^2}}{Re}$",'Interpreter','latex','FontSize',22)
saveas(gcf,'Om1B30_mdisp_bar.pdf')

figure('color','white') 

nd = 766^2;

% this reads all data tim from multiple timestamps into 1 array (it ignores the lineskips in the .dat file)
vort = readmatrix("Om3B30.dat");
invRo = 3;

s1 = size(vort);
sd = s1(1)/nd

% in order to plot this by individual timestep, you have to restrict the plot by column 3 which is the timestep. 

time_col = 3;
[row,col] = find(vort(:,time_col) == vort(s1(1)-2*nd,time_col));
s2 = size(vort(row,time_col));

np = 64; colors = zeros(np);

xcol_num = 7;
ycol_num = 8;

min1 = min(vort(:,xcol_num));
max1 = max(vort(:,xcol_num));
min2 = min(vort(:,ycol_num));
max2 = max(vort(:,ycol_num));

d1 = (max1 - min1)/(np-1);
d2 = (max2 - min2)/(np-1);

for i = 1:s1(1)
    x = floor((vort(i,xcol_num)-min1)/d1)+1;
    y = np-floor((vort(i,ycol_num)-min2)/d2);
    % the plotting with imagesc is weird (x-y are switched)
    colors(y,x) = colors(y,x) + 1;
end
colors = colors/sd;

plot(vort(:,xcol_num)+invRo,vort(:,ycol_num), 'o')
xlim([min1, max1]+2);

figure('color','white')

%xlabel("vorticity")
%ylabel("Thermal Dissipation")
imagesc(log(colors))
xtickformat('%.5f');
ytickformat('%.5f');
num_ticks = 10;
xticks(linspace(0,64,num_ticks+1));
yticks(linspace(0,64,num_ticks+1));
xticklabels(round(linspace(min1,max1,num_ticks+1)+invRo,2));
yticklabels(round(flip(linspace(min2,max2,num_ticks+1)),2));
colormap sky
colorbar

figure('color','white')

x1 = 1:766;
y1 = 1:766;

z1 = reshape(vort(row,7),[766, 766]);
z2 = abs(invFr./(z1+invRo));
z3 = reshape(vort(row,8),[766, 766])/2;
z4 = reshape(vort(row,9),[766, 766]);
z5 = reshape(vort(row,12),[766, 766])/(4*600);


surf(x1, y1, z1+invRo)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
%colormap(redwhiteblue)
colormap(redwhiteblue(min(z1+invRo,[],'all'),max(z1+invRo,[],'all')))
colorbar
%clim([-35,35])
title("$\widehat{\omega_z} + 2\Omega$",'Interpreter','latex','FontSize',22)
view(0,90)
saveas(gcf,'Om3B30_vortz_bar.pdf')

figure('color','white')

surf(x1, y1, z2)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
clim([0 1])
colorbar
title("$\frac{Ro_L}{Fr}$",'Interpreter','latex','FontSize',22)
saveas(gcf,'Om3B30_RC.pdf')

figure('color','white')

surf(x1, y1, z3)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
colorbar
%clim([0,3.5])
title("$\frac{\widehat{|\nabla T|^2}}{Fr^2Pe}$",'Interpreter','latex','FontSize',22)
view(0,90)
saveas(gcf,'Om3B30_tdisp_bar.pdf')

figure('color','white')

surf(x1, y1, z4)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
colorbar
%clim([0,1.2])
title("$\widehat{w^2}$", 'Interpreter','latex','FontSize',22)
saveas(gcf,'Om3B30_uzrms_bar.pdf')

figure('color','white')

surf(x1, y1, z5)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
colorbar
%clim([0,1.2])
title("$\frac{\widehat{|\nabla \mathbf{u}|^2}}{Re}$",'Interpreter','latex','FontSize',22)
saveas(gcf,'Om3B30_mdisp_bar.pdf')

figure('color','white') 

nd = 766^2;

% this reads all data tim from multiple timestamps into 1 array (it ignores the lineskips in the .dat file)
vort = readmatrix("Om10B30.dat");
invRo = 10;

s1 = size(vort);
sd = s1(1)/nd

% in order to plot this by individual timestep, you have to restrict the plot by column 3 which is the timestep. 

time_col = 3;
[row,col] = find(vort(:,time_col) == vort(s1(1),time_col));
s2 = size(vort(row,time_col));

np = 64; colors = zeros(np);

xcol_num = 7;
ycol_num = 8;

min1 = min(vort(:,xcol_num));
max1 = max(vort(:,xcol_num));
min2 = min(vort(:,ycol_num));
max2 = max(vort(:,ycol_num));

d1 = (max1 - min1)/(np-1);
d2 = (max2 - min2)/(np-1);

for i = 1:s1(1)
    x = floor((vort(i,xcol_num)-min1)/d1)+1;
    y = np-floor((vort(i,ycol_num)-min2)/d2);
    % the plotting with imagesc is weird (x-y are switched)
    colors(y,x) = colors(y,x) + 1;
end
colors = colors/sd;

plot(vort(:,xcol_num)+invRo,vort(:,ycol_num), 'o')
xlim([min1, max1]+2);

figure('color','white')

%xlabel("vorticity")
%ylabel("Thermal Dissipation")
imagesc(log(colors))
xtickformat('%.5f');
ytickformat('%.5f');
num_ticks = 10;
xticks(linspace(0,64,num_ticks+1));
yticks(linspace(0,64,num_ticks+1));
xticklabels(round(linspace(min1,max1,num_ticks+1)+invRo,2));
yticklabels(round(flip(linspace(min2,max2,num_ticks+1)),2));
colormap sky
colorbar

figure('color','white')

x1 = 1:766;
y1 = 1:766;

z1 = reshape(vort(row,7),[766, 766]);
z2 = abs(invFr./(z1+invRo));
z3 = reshape(vort(row,8),[766, 766])/2;
z4 = reshape(vort(row,9),[766, 766]);
z5 = reshape(vort(row,12),[766, 766])/(4*600);


surf(x1, y1, z1+invRo)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
%colormap(redwhiteblue)
colormap(redwhiteblue(min(z1+invRo,[],'all'),max(z1+invRo,[],'all')))
colorbar
%clim([-35,35])
title("$\widehat{\omega_z} + 2\Omega$",'Interpreter','latex','FontSize',22)
view(0,90)
saveas(gcf,'Om10B30_vortz_bar.pdf')

figure('color','white')

surf(x1, y1, z2)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
clim([0 1])
colorbar
title("$\frac{Ro_L}{Fr}$",'Interpreter','latex','FontSize',22)
saveas(gcf,'Om10B30_RC.pdf')

figure('color','white')

surf(x1, y1, z3)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
colorbar
%clim([0,3.5])
title("$\frac{\widehat{|\nabla T|^2}}{Fr^2Pe}$",'Interpreter','latex','FontSize',22)
view(0,90)
saveas(gcf,'Om10B30_tdisp_bar.pdf')

figure('color','white')

surf(x1, y1, z4)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
colorbar
%clim([0,1.2])
title("$\widehat{w^2}$", 'Interpreter','latex','FontSize',22)
saveas(gcf,'Om10B30_uzrms_bar.pdf')

figure('color','white')

surf(x1, y1, z5)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
colorbar
%clim([0,1.2])
title("$\frac{\widehat{|\nabla \mathbf{u}|^2}}{Re}$",'Interpreter','latex','FontSize',22)
saveas(gcf,'Om10B30_mdisp_bar.pdf')

vort = readmatrix("Om1B100.dat");
invRo = 1;
invFr = sqrt(100);

s1 = size(vort)
sd = s1(1)/nd

% in order to plot this by individual timestep, you have to restrict the plot by column 3 which is the timestep. 
num_timesteps = 766^2;

time_col = 3;
[row,col] = find(vort(:,time_col) == vort(s1(1),time_col));
s2 = size(vort(row,time_col));

np = 64;
colors = zeros(np);

xcol_num = 7;
ycol_num = 8;

min1 = min(vort(:,xcol_num));
max1 = max(vort(:,xcol_num));
min2 = min(vort(:,ycol_num));
max2 = max(vort(:,ycol_num));
d1 = (max1 - min1)/(np-1);
d2 = (max2 - min2)/(np-1);

for i = 1:s1(1)
    x = floor((vort(i,xcol_num)-min1)/d1)+1;
    y = np-floor((vort(i,ycol_num)-min2)/d2);
    % the plotting with imagesc is weird (x-y are switched)
    colors(y,x) = colors(y,x) + 1;
end
colors = colors/sd;

figure('color','white')

plot(vort(:,xcol_num)+invRo,vort(:,ycol_num), 'o')
xlim([min1, max1]+2);

figure('color','white')

imagesc(log(colors))
xtickformat('%.5f');
ytickformat('%.5f');
num_ticks = 10;
xticks(linspace(0,64,num_ticks+1));
yticks(linspace(0,64,num_ticks+1));
xticklabels(round(linspace(min1,max1,num_ticks+1)+invRo,2));
yticklabels(round(flip(linspace(min2,max2,num_ticks+1)),2));
colormap sky
xlabel("$\overline{\omega_z} + 2\Omega$", 'Interpreter','latex','Fontsize',22)
ylabel("$\overline{|\nabla T|^2}$", 'Interpreter','latex', 'Rotation', 0, 'FontSize',22)

figure('color','white')

x1 = 1:766;
y1 = 1:766;

z1 = reshape(vort(row,7),[766, 766]);
z2 = abs(invFr./(z1+invRo));
z3 = reshape(vort(row,8),[766, 766])/2;
z4 = reshape(vort(row,9),[766, 766]);
z5 = reshape(vort(row,12),[766, 766])/(4*600);


surf(x1, y1, z1+invRo)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap(redwhiteblue(min(z1+invRo,[],'all'),max(z1+invRo,[],'all')))
colorbar
title("$\widehat{\omega_z} + 2\Omega$", 'Interpreter','latex','FontSize',22)
%clim([-35,35])
saveas(gcf,'Om1B100_vortz_bar.pdf')

figure('color','white')

surf(x1, y1, z2)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
clim([0 1])
colorbar
title("$\frac{Ro_L}{Fr}$",'Interpreter','latex','FontSize',22)
saveas(gcf,'Om1B100_RC.pdf')

figure('color','white')

surf(x1, y1, z3)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
colorbar
%clim([0,3.5])
title("$\frac{\widehat{|\nabla T|^2}}{Fr^2Pe}$",'Interpreter','latex','FontSize',22)
saveas(gcf,'Om1B100_tdisp_bar.pdf')

figure('color','white')

surf(x1, y1, z4)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
colorbar
%clim([0,1.2])
title("$\widehat{w^2}$",'Interpreter','latex','FontSize',22)
saveas(gcf,'Om1B100_uzrms_bar.pdf')

figure('color','white')

surf(x1, y1, z5)
xlim([0,766])
ylim([0,766])
set(gca,'XTick',[],'YTick',[])
view(0,90)
shading interp
colormap parula
colorbar
%clim([0,1.2])
title("$\frac{\widehat{|\nabla \mathbf{u}|^2}}{Re}$",'Interpreter','latex','FontSize',22)
saveas(gcf,'Om1B100_mdisp_bar.pdf')

figure('color','white') 

