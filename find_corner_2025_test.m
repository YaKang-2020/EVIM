function [idx_opt,lambda_opt]=find_corner_2025_test(Lres_a,Lbeta_a,aspan)
% This program aims to find the maximum curvature of the L-curve
% input:
% Lres_a: misfit norm, norm(A*h - b)
% Lbeta: smoothness norm, norm(B*h);
% aspan: filtering parameters

% output:
%  idx_opt: The location of the inflection point in aspan
%  lambda_opt: Optimal filtering parameters
% Copyright (C) 2025 Ya Kang, [kangya@njupt.edu.cn].
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
Lres_all=Lres_a';
Lbeta_all=Lbeta_a';
Lres_all_n=log10(Lres_all);
Lbeta_all_n=log10(Lbeta_all);

nn=length(Lres_all_n);

Lres_int=Lres_all_n(1):0.001:Lres_all_n(nn);

Lbeta_int=interp1(Lres_all_n,Lbeta_all_n,Lres_int,'cubic');
aspan_int=interp1(Lres_all_n,aspan,Lres_int,'cubic');


lambdas=aspan_int;

dlog_rho = gradient(Lres_int); % 
dlog_eta = gradient(Lbeta_int); % 
ddlog_rho = gradient(dlog_rho); % 
ddlog_eta = gradient(dlog_eta); % 

numerator = dlog_rho .* ddlog_eta - dlog_eta .* ddlog_rho;
denominator = (dlog_rho.^2 + dlog_eta.^2).^(3/2);
% curvature = abs(numerator) ./ denominator; % |abs|
curvature = numerator ./ denominator; % |abs|
[max_curve, idx_opt] = max(curvature);

% Lres_int_opt=Lres_int(idx_opt)

lambda_opt = lambdas(idx_opt);


% plot(curvature)


%  Lres_all_n=log10(Lres_a);
% Lbeta_all_n=log10(Lbeta_a);

x_max=Lres_int(idx_opt)
y_max=Lbeta_int(idx_opt)
figure, plot(Lres_int,Lbeta_int)
hold on;
plot(x_max,y_max,'rp')


