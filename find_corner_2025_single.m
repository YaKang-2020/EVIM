function [idx_opt,lambda_opt]=find_corner_2025_single(Lres_a,Lbeta_a,aspan)
% This program aims to obtain the maximum curvature of the L-curve
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
lambdas=aspan;

% dlog_rho = gradient(Lres_all_n, log(lambdas)); % d(logρ)/d(logλ)
dlog_eta = gradient(Lbeta_all_n, Lres_all_n); % d(logη)/d(logλ)
% ddlog_rho = gradient(dlog_rho, log(lambdas)); % d²(logρ)/d(logλ)²
ddlog_eta = gradient(dlog_eta, Lres_all_n); % d²(logη)/d(logλ)²

% numerator = dlog_rho .* ddlog_eta - dlog_eta .* ddlog_rho;
denominator = (1 + dlog_eta.^2).^(3/2);
% curvature = abs(numerator) ./ denominator; % |abs|
curvature = ddlog_eta ./ denominator; % |abs|
[max_curve, idx_opt] = max(curvature);
lambda_opt = lambdas(idx_opt);

