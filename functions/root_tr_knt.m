function [p] = root_tr_knt(Ka, Kd, num_k, mode)
% ROOT_TR_KNT roots of the transcendental equation for kinetic model
%   [p] = tr_root_eq(Ka, Kd, num_k) calculates the first num_k roots of the
%   transcendental equation for the equilibrium model.
%
% INPUT:
%   Ka -- 1*num_set array of adsorption rate constants
%   Kd -- 1*num_set array of desorption rate constants
%   num_k -- num of roots calculated
%
% OUTPUT:
%   p -- num_k*num_set matrix of roots for the transcendental equation
%

switch nargin
    case 2
        num_k = 2;
        mode = 'zero';    
    case 3
        mode = 'zero';
end

num_set = length(Ka);
p = zeros(num_k, num_set);      % include the trial zero roots as the first row

options = optimset('TolX',1e-12);

eps = 1e-14;

for i = 1:num_set
    for k = 2:num_k
        lb = (k-2)*pi;
        rb = (k-1)*pi;
        
        if Kd(i)==0
            p(k, i) = fzero(@(x) (tan(x)*(x^2 - Kd(i))/x - Ka(i)), [rb - pi/2 + eps, rb + pi/2 - eps], options);            
        elseif Kd(i) < ((2*(k-1)-1)*pi/2)^2
            p(k, i) = fzero(@(x) (tan(x)*(x^2 - Kd(i))/x - Ka(i)), [lb + eps, lb + pi/2 - eps], options);
        else
            p(k, i) = fzero(@(x) (tan(x)*(x^2 - Kd(i))/x - Ka(i)), [rb - pi/2 + eps, rb - eps], options);
        end
    end 
end

if strcmpi(mode, 'nonzero')
    p = p(2:num_k, :);
end