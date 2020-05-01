function [M_LM, B, V, D] = moment_LM(moment, Ka, Para)
% MOMENT_LM normalized global moments for adsorption-only model
%   m_LM = moment_LM(moment, Ka, Para) calculates the zeroth, first, second
%   moments for transient transport between two plates with Poiseuille flow
%   and adsorption-only model (Lungu&Moffatt model) using the results from
%   LM to get the decay rate for zeroth moment, slope for first moment and
%   slope for second moment.
%
% INPUT:
%   moment -- 0 or 1 or 2, assigning up to which moment you want calculate.
%   Ka -- 1*num_set array of adsorption rate constants.
%   Para -- structure including the physical and computational parameters
%
% OUTPUT:
%   M_LM -- (num_t * num_set * moment+1) matrix for the moments.
%   B -- 1*num_set of the decay rate for zeroth moment
%   V -- 1*num_set of the dimensionless asymptotic velocity for the first
%   moment.
%   D -- 1*num_set of the dispersion coefficient for the second moment
%

num_set = length(Ka); 
num_t = length(Para.tD);
B = zeros(1, num_set);
V = zeros(1, num_set);
D = zeros(1, num_set);
M_LM = zeros(num_t, num_set, moment + 1);

q = zeros(num_set,1);
for i = 1:num_set
    q(i) = fzero(@(x)x*tan(x) - Ka(i), [0,pi/2]);
end

% ---- zeroth moment ----
if moment >= 0
    B(1, :) = q.^2;
    M_LM(:,:,1) = Para.m_ini*exp(-q.^2*Para.tD)';
end

% ---- first moment ----
if moment >= 1

    c = cos(2*q); s = sin(2*q);
    p01 =(8.*q.^3-6.*q.*c+3.*s)./(6.*q.^2.*(2.*q+s)).*10;

    V = p01/10.*3/2;
    
    M_LM(:,:,2) = (V * Para.Pe * Para.tD)';
    
end

% ---- second moment ----
if moment >= 2

    syms x y
    c = cos(2*x); s = sin(2*x);
    p = (8*x.^3-6*x.*c+3*s)./(6*x.^2.*(2*x+s));
    g = (2*x.^2+3*c+3*x.*s+3)/6./x.^2./(2*x+s);
    Y = (1-y.^2-p).*(-y.^2/4./x.^2.*cos(x.*y)+(g.*y-y.^3/6./x).*sin(x.*y)).*cos(x*y);

    M = -1./(1+sin(2*x)/2./x) .* int(Y,y,-1,1);
    M = matlabFunction(M);
    
    D = 9/4*M(q) * Para.Pe^2 + 1;

    M_LM(:,:,3) = 2 * (D * Para.tD)';

    
end