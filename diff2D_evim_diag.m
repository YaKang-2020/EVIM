function [d2] = diff2D_evim_diag(nrows, ncols)

%Generates sparse matrix of -6, 1's, and 1/2's for calculating second 
%finite difference includian diagonals (laplacian) in 2D using a 9 point 
%scheme.

% [d] = diff2D_evim_diag(5, 15)
% Copyright (C) 2025 Ya Kang, [kangya@njupt.edu.cn].
%-------------------------------------------------------------------------%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>.
%-------------------------------------------------------------------------%

len = nrows * ncols;  % 
d2 = sparse(len, len); % 

d2 = d2 + sparse(1:len, 1:len, -6*ones(len,1), len, len);


d2 = d2 + spdiags(ones(len,1), -nrows, len, len);   % 
d2 = d2 + spdiags(ones(len,1), -1, len, len);       % 
d2 = d2 + spdiags(ones(len,1), 1, len, len);        % 
d2 = d2 + spdiags(ones(len,1), nrows, len, len);    % 

d2 = d2 + spdiags(0.5*ones(len,1), -nrows-1, len, len);  % 
d2 = d2 + spdiags(0.5*ones(len,1), -nrows+1, len, len);  % 
d2 = d2 + spdiags(0.5*ones(len,1), nrows-1, len, len);   % 
d2 = d2 + spdiags(0.5*ones(len,1), nrows+1, len, len);   % 

end