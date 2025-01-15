close all; clear; clc;

format shortE

np = 766^2;

% this reads all data tim from multiple timestamps into 1 array (it ignores the lineskips in the .dat file)
vort = readmatrix("new_Om2B100.dat");
invRo = 2;

s1 = size(vort);

% in order to plot this by individual timestep, you have to restrict the plot by column 3 which is the timestep. 
time_col = 3
last_timestep = vort(s1(1),time_col);

[row, col] = find(vort(:,time_col) == last_timestep);

s2 = size(vort(row,time_col));

np = 64;
colors = zeros(np);

xcol_num = 7
ycol_num = 8

min1 = min(vort(row,xcol_num));
max1 = max(vort(row,xcol_num));
min2 = min(vort(row,ycol_num));
max2 = max(vort(row,ycol_num));

d1 = (max1 - min1)/(np-1);
d2 = (max2 - min2)/(np-1);

for i = 1:s2(1)
    x = floor((vort(row(i),xcol_num)-min1)/d1)+1;
    y = np-floor((vort(row(i),ycol_num)-min2)/d2);
    % the plotting with imagesc is weird (x-y are switched)
    colors(y,x) = colors(y,x) + 1;
end

plot(vort(row,xcol_num)+invRo,vort(row,ycol_num), 'o')
xlim([min1, max1]+2);

figure

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

figure

x1 = 1:766;
y1 = 1:766;

z1 = reshape(vort(row,7),[766, 766]);
z2 = reshape(vort(row,11),[766, 766]);
z3 = reshape(vort(row,8),[766, 766]);
z4 = reshape(vort(row,9),[766, 766]);


contour(x1, y1, z1,20)
colormap parula
colorbar
title("$\frac{Ro_L}{Fr}$", "Interpreter", "Latex")

%figure

%contour(x1, y1, z2,20)
%colormap parula
%colorbar
%title("Barotropic Energy Ratio")

figure

contour(x1, y1, z3,20)
colormap parula
colorbar
title("Thermal Dissipation")

figure

contour(x1, y1, z4,20)
colormap parula
colorbar
title("uz_rms")

figure 

vort = readmatrix("new_Om1B100.dat");
invRo = 1;

s1 = size(vort);

% in order to plot this by individual timestep, you have to restrict the plot by column 3 which is the timestep. 
time_col = 3
last_timestep = vort(s1(1),time_col);

[row, col] = find(vort(:,time_col) == last_timestep);

s2 = size(vort(row,time_col));

np = 64;
colors = zeros(np);

xcol_num = 7
ycol_num = 8

min1 = min(vort(row,xcol_num));
max1 = max(vort(row,xcol_num));
min2 = min(vort(row,ycol_num));
max2 = max(vort(row,ycol_num));

d1 = (max1 - min1)/(np-1);
d2 = (max2 - min2)/(np-1);

for i = 1:s2(1)
    x = floor((vort(row(i),xcol_num)-min1)/d1)+1;
    y = np-floor((vort(row(i),ycol_num)-min2)/d2);
    % the plotting with imagesc is weird (x-y are switched)
    colors(y,x) = colors(y,x) + 1;
end

plot(vort(row,xcol_num)+invRo,vort(row,ycol_num), 'o')
xlim([min1, max1]+2);

figure

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

z1 = reshape(vort(row,7),[766, 766]);
z2 = reshape(vort(row,11),[766, 766]);
z3 = reshape(vort(row,8),[766, 766]);
z4 = reshape(vort(row,9),[766, 766]);


contour(x1, y1, z1,20)
colormap parula
colorbar
title("Vorticity")

figure

contour(x1, y1, z2,20)
colormap parula
colorbar
title("Barotropic Energy Ratio")

figure

contour(x1, y1, z3,20)
colormap parula
colorbar
title("Thermal Dissipation")

figure

contour(x1, y1, z4,20)
colormap parula
colorbar
title("uz_rms")

