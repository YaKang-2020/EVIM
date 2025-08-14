 function [A] = evim_A(u,v,dx)

%[A] = evim_A(u,v,dx)
%
%Generates the matrix A for inversion of the continuity equation to
%determine landslide thicknesses

% u = west-east deformation
% v = north-south deformation
% dx = pixel spacing
% A = design matrix for inversion 
% 
%
%Copyright (C) 2025 Ya Kang, [kangya@njupt.edu.cn].
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
%
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>.
%-------------------------------------------------------------------------%

%Sizes,

[nrows,ncols] = size(u);
len = nrows*ncols;

%Srtain, in this version, we removed the calculation of strain 20250101,by set zeros
strain_u=zeros(nrows,ncols); % in this version, we removed the calculation of strain 20250101,by set zeros
strain_v=zeros(nrows,ncols); % in this version, we removed the calculation of strain 20250101,by set zeros
strain_s=(1+strain_u).*(1+strain_v); % in this version, we removed the calculation of strain 20250101,by set zeros
% Replace NaN entries with zeros in u and v,
nanmask = ones(nrows,ncols);
nanmask(isnan(u) == 1 | isnan(v) == 1) = NaN;
u = u.*nanmask;
v = v.*nanmask;
u(isnan(u)) = 0;
v(isnan(v)) = 0;

%Reshape arrays to vectors,
u = reshape(u,len,1);
v = reshape(v,len,1);

 uv0=2*dx*((1-strain_s)./strain_s);
 uv0 = reshape(uv0,len,1); % third term
%Diagonals -nrows and +nrows from u,
u(1:nrows) = 0;                             %
uup = zeros(len,1);
uup(nrows+1:len) = -u(1:len-nrows);
u(len-nrows+1:len) = 0;                     %
udn = zeros(len,1); %
udn(1:len-nrows) = 1*u(nrows+1:len);
%Diagonals -1 and 1 from v,
v(nrows:nrows:len) = 0;                     %
v(1:nrows:len-nrows+1,:) = 0;               %
vup = zeros(len,1);
vup(2:len) = 1*v(1:len-1);
vdn = zeros(len,1);
vdn(1:len-1) =- v(2:len);

% create A,
 A = spdiags([udn vdn uv0 vup uup],[-nrows -1 0 1 nrows],len,len);


