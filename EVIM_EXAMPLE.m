%EVIM_EXAMPLE.m 
%
% An example of using Evim to invert landslide thickness based on simulated
% landslide scenarios. In the simulation, we assume that the deformation of 
% the landslide is parallel to the sliding surface. For more details, 
% please refer to the  Kang et al., 2025.
%
% Copyright (C) 2025 Ya Kang, [kangya@njupt.edu.cn].
%
%-------------------------------------------------------------------------%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%1
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>.
%-------------------------------------------------------------------------%
clc,clear
ddp=0.03 ;  % filtering(damping) factor
load ./evim_test                           %  u, v, dzdt are coordinates in meters,we can use the gamma_GBI.m to get this *mat
% where u is east-west deformation field
%       v is north-south deformation field
%       dzdt is vertical deformation field
%       surf_terrain is DEM
%       mask is the mask file, 1 indicates the landslide area, 0 indicates the non-landslide area
%       
%% Setting the parameter
dx = 5;   % pixel spacing of the image
up_bod=150;  % up boundary depth for inversion of the slip surface
low_bod=0;   % low boundary depth for inversion of the slip surface
curl=0;  % whether eastimate the optimization filtering parameters - ddp, 1 > Yes, 0 > no need

%% Data preparation
mask(isnan(mask)) = 0;          %Replaces no data with zero at
u(isnan(u)) = 0;                %the edges of the grids.
v(isnan(v)) = 0;
dzdt(isnan(dzdt)) = 0;

% plot the horizontal deformation field
[m,n]=size(u);
xx=zeros(m,n);
 for i=1:m
xx(i,:)=-(n-1)/2:(n-1)/2;
 end
x=dx*xx;
yy=zeros(m,n);
 for i=1:m
yy(:,i)=-(m-1)/2:(m-1)/2;
 end
y=-dx*yy;
figure; 
quiver(x,y,u,v);               % plot the horizontal deformation
axis equal tight;

%Creat the design matrix for inversion
[A] = evim_A(u,v,dx);   %  the design matrix for inversion
fs1=fspecial('average',[5 5]);  % DEM filtering 
surf_terrain_f=imfilter(surf_terrain, fs1, 'replicate'); % DEM filtering 

 % aspan is the smooth factor, dzdt is the vertical deformtaion, dx is the pixel spacing,
if curl==1
aspan =  10.^linspace(-2,1,30); % The range of filtering parameters could be set according to the actual situation,e.g. 10.^linspace(-3,1,40); 10.^linspace(-4,1,50)
Lres_a = zeros(length(aspan),1);          %Pre-allocate vectors for the L-
Lbeta_a = zeros(length(aspan),1); 
parfor i = 1:length(aspan)
        [htemp(:,:,i),Lres_a(i),Lbeta_a(i)] = ...
       evim_smooth(A,dzdt,mask,dx,aspan(i),surf_terrain_f,up_bod,low_bod,u,v);
end
figure, plot(log10(Lres_a),log10(Lbeta_a));  %plot the L-curve 
% [idx_opt_l,lambda_opt_l]=find_corner_2025_single(Lres_a,Lbeta_a,aspan); 
[idx_opt_l,lambda_opt_l]=find_corner_2025_test(Lres_a,Lbeta_a,aspan); % L-curve to find the optimal damping factor
ddp=lambda_opt_l; 
end

%% Inversion
[h,Lres,Lbeta] = evim_smooth(A,dzdt,mask,dx,ddp,surf_terrain_f,up_bod,low_bod,u,v);
% mask_new=mask_bmp(az_start:az_end,r_start:r_end);
h(mask==0)=0;
slip_ele=surf_terrain-h;  % locate the slip-surface
slip_ele(mask==0)=0;      % slip-surface elevation

%% Results drawing
figure,imagesc(h)
title('Inverted depth (m)');colorbar
set(gca,'xtick',[0:25:100])
set(gca,'xticklabel',[0:25*dx:100*dx])
set(gca,'ytick',[0:20:100])
set(gca,'yticklabel',[100*dx:-20*dx:0])
set(gca,'FontSize',36, 'Fontname','Times New Roman')
ylabel('y-coordinate (m)','fontsize',36)
xlabel('x-coordinate (m)','fontsize',36)
set(gcf, 'position', [100, 100, 1100, 900]);
caxis([-70,70])
% colormap(flipud(parula))
saveas(gcf, 'Inverted depth.jpg')

figure,imagesc(slip_ele)
title('Slip surface elevation (m)');colorbar
set(gca,'xtick',[0:25:100])
set(gca,'xticklabel',[0:25*dx:100*dx])
set(gca,'ytick',[0:20:100])
set(gca,'yticklabel',[100*dx:-20*dx:0])
set(gca,'FontSize',36, 'Fontname','Times New Roman')
ylabel('y-coordinate (m)','fontsize',36)
xlabel('x-coordinate (m)','fontsize',36)
set(gcf, 'position', [100, 100, 1100, 900]);
% caxis([-70,70])
% colormap(flipud(parula))
saveas(gcf, 'Slip surface elevation.jpg')

%% Comparison with true value (only for the Simulated landslides)
% load the Simulated true value of the landslide thickness for comparison, for the real
% 3D data (e.g. SAR, InSAR 3D deformation) the following operations cannot be performed.  
load ./real_depth  % Simulated true value of the landslide thickness

figure,imagesc(depth_l)
title('Real depth (m)');colorbar
set(gca,'xtick',[0:25:100])
set(gca,'xticklabel',[0:25*dx:100*dx])
set(gca,'ytick',[0:20:100])
set(gca,'yticklabel',[100*dx:-20*dx:0])
set(gca,'FontSize',36, 'Fontname','Times New Roman')
ylabel('y-coordinate (m)','fontsize',36)
xlabel('x-coordinate (m)','fontsize',36)
set(gcf, 'position', [100, 100, 1100, 900]);
caxis([-70,70])
% colormap(flipud(parula))
saveas(gcf, 'Real depth.jpg')

difference=depth_l-h; % the Residual of the inverted landslide thickness

figure,imagesc(difference)
title('Residual(m)');colorbar
set(gca,'xtick',[0:25:100])
set(gca,'xticklabel',[0:25*dx:100*dx])
set(gca,'ytick',[0:20:100])
set(gca,'yticklabel',[100*dx:-20*dx:0])
set(gca,'FontSize',36, 'Fontname','Times New Roman')
ylabel('y-coordinate (m)','fontsize',36)
xlabel('x-coordinate (m)','fontsize',36)
set(gcf, 'position', [100, 100, 1100, 900]);
caxis([-70,70])
% colormap(flipud(parula))
saveas(gcf, 'Residual.jpg')

res_tmp=difference(31:70,11:90);
res_kk=res_tmp(:);
rmsek=sqrt(mean(res_kk.^2));
disp('The RMSE of the results is:')
disp(rmsek)

