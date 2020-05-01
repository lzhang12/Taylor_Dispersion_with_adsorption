function [p] = root_tr_eq(Kp, num_k, mode)
% TR_ROOT_EQ roots of the transcendental equation for eq model
%   [p] = tr_root_eq(Kp, num_k) calculates the first num_k roots of the
%   transcendental equation for the equilibrium model.
%
% INPUT:
%   Kp -- 1*num_set array of retention factors.
%   num_k -- num of roots calculated
%
% OUTPUT:
%   p -- num_k*num_set matrix of roots for the transcendental equation
%

switch nargin
    case 1
        num_k = 2;
        mode = 'zero';    
    case 2
        mode = 'zero';
end

num_set = length(Kp);
p = zeros(num_k, num_set);      % include the trial zero roots as the first row

options = optimset('TolX',1e-25);

eps = 1e-14;

for i = 1:num_set
    for k = 2:num_k         % start from the first non-zero root
        lb = (k-2)*pi + pi/2;
        rb = (k-1)*pi + pi/2;
        p(k, i) = fzero(@(x) tan(x)/x + Kp(i), [eps + lb, rb - eps], options);
    end
end

if strcmpi(mode, 'nonzero')
    p = p(2:num_k, :);
end
    