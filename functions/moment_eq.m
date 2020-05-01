function M_eq = moment_eq(moment, method, Kp, Para, num_k)
% MOMENT_EQ normalized global moments for equilibrium model
%   m_eq = moment_eq(moment, method, Kp, Para) calculates the zeroth,
%   first, second normalized global moments for transient transport between
%   two plates with Poiseuille flow and equilibrium surface reaction model
%   (isotherm model) using either numerical inverse Laplace transform, or
%   truncated series solution derived by residue theorem.
%
% INPUT:
%   moment -- 0 or 1 or 2, assigning up to which moment you want calculate.
%   method -- 'num' to use the numerical inverse Laplace transform, or
%             'ana' to use the truncated analytical solution.
%   Kp -- 1*num_set array of retention factors.
%   Para -- structure including the physical and computational parameters
%   num_k --number of roots used for series solution
%
% OUTPUT:
%   M_eq -- num_t * num_set * moment matrix for the moments
%

if(nargin < 5)
    num_k = 2;
end

num_set = length(Kp); 
num_t = length(Para.tD);
m0_eq = zeros(num_t, num_set);
m1_eq = zeros(num_t, num_set);
m2_eq = zeros(num_t, num_set);
M_eq = zeros(num_t, num_set, moment + 1);

Kp_mat = repmat(Kp, num_k, 1);
tD_mat = repmat(Para.tD, num_k, 1);

% % If mupad has not been evaluated, evaluate to the end
% Eval_Mark = getVar(nb, 'Eval_Mark');
% if ~strcmp(char(Eval_Mark), 'yes');
%     evaluateMuPADNotebook(nb);
% end

% ------ Numerical Inverse Laplace Transform ------
num_sum = 20;

if strcmpi(method,'num') == 1
    FN_mupad = './mupad/moment_eq.mn';
    nb = mupad(FN_mupad);
    
    % zeroth moment
    if moment >= 0
        Lm0_ana = getVar(nb, 'Lm0_ana');
        Lm0_ana = matlabFunction(Lm0_ana);
        for i = 1:num_set
            fun = @(s) Lm0_ana(Kp(i), Para.m_ini, s);
            m0_eq(:, i) = talbot_inversion(fun, Para.tD, num_sum);
        end
        
        M_eq(:, :, 1) = m0_eq./Para.m_ini;
    end
    
    % first moment
    if moment >= 1
        Lm1_ana = getVar(nb, 'Lm1_ana');
        Lm1_ana = matlabFunction(Lm1_ana);
        for i = 1:num_set
            fun = @(s) Lm1_ana(Kp(i), Para.Pe, Para.m_ini, s);
            m1_eq(:, i) = talbot_inversion(fun, Para.tD, num_sum);
        end
        
        M_eq(:, :, 2) = m1_eq./m0_eq;
    end
    
    % second moment
    if moment >= 2
        Lm2_ana = getVar(nb, 'Lm2_ana');
        Lm2_ana = matlabFunction(Lm2_ana);
        for i = 1:num_set
            fun = @(s) Lm2_ana(Kp(i), Para.Pe, Para.m_ini, s, Para.var0);
            m2_eq(:, i) = talbot_inversion(fun, Para.tD, num_sum);
        end
        
        M_eq(:, :, 3) = m2_eq./m0_eq - (m1_eq./m0_eq).^2;
    end

% ------ Truncated Analytical Solution ------
elseif strcmpi(method, 'series') == 1
    
    p = root_tr_eq(Kp, num_k);
    
    if moment >= 0
        
        a_eq = (2*Kp_mat.^2*Para.m_ini)./(Kp_mat.^2.*p.^2 + Kp_mat + 1);
        a_eq(1, :) = Para.m_ini./(Kp + 1);      % coefficients for zero root
        
        for i = 1:num_set
            m0_eq(:, i) = a_eq(:, i)' * exp(-p(:,i).^2*Para.tD);
        end
        
        M_eq(:, :, 1) = m0_eq./Para.m_ini;
    end
    
    if moment >= 1
        bk2_eq = (Kp_mat.^2*Para.Pe*Para.m_ini.*(4*Kp_mat.^2.*p.^4 + 3*Kp_mat.^2.*p.^2 - 3*Kp_mat + 4*p.^2 - 3))./(2*p.^2.*(Kp_mat.^2.*p.^2 + Kp_mat + 1).^2);
        bk1_eq = -(Kp_mat*Para.Pe*Para.m_ini.*(-6*Kp_mat.^5.*p.^6 + 15*Kp_mat.^5.*p.^4 + 16*Kp_mat.^4.*p.^4 + 51*Kp_mat.^4.*p.^2 - 8*Kp_mat.^3.*p.^4 + 93*Kp_mat.^3.*p.^2 + 30*Kp_mat.^3 + 40*Kp_mat.^2.*p.^2 + 84*Kp_mat.^2 - 2*Kp_mat.*p.^2 + 78*Kp_mat + 24))./...
               (2*p.^4.*(Kp_mat.^2.*p.^2 + Kp_mat + 1).^3);
        bk2_eq(1, :) = Para.Pe*Para.m_ini./(Kp + 1).^2;
        bk1_eq(1, :) = 2*Kp*Para.Pe*Para.m_ini.*(6*Kp+1)./(15*(Kp+1).^3);
        
        for i = 1:num_set
            m1_eq(:, i) = bk1_eq(:, i)' * exp(-p(:,i).^2*Para.tD) + bk2_eq(:, i)' * ( tD_mat .* exp(-p(:,i).^2*Para.tD));
        end
        
        M_eq(:, :, 2) = m1_eq./m0_eq;
    end
    
    if moment >= 2
        
        FN_mupad = './mupad/Res_limit_eq.mn';
        nb = mupad(FN_mupad);

        f_ck3_eq = getVar(nb, 'ck3_eq');
        f_ck3_eq = matlabFunction(f_ck3_eq);
        ck3_eq = f_ck3_eq(Kp_mat, Para.Pe, Para.m_ini, p);

        f_ck2_eq = getVar(nb, 'ck2_eq');
        f_ck2_eq = matlabFunction(f_ck2_eq);
        ck2_eq = f_ck2_eq(Kp_mat, Para.Pe, Para.m_ini, p);

        f_ck1_eq = getVar(nb, 'ck1_eq');
        f_ck1_eq = matlabFunction(f_ck1_eq);
        ck1_eq = f_ck1_eq(Kp_mat, Para.Pe, Para.m_ini, p);
        
        ck3_eq(1, :) = Para.Pe^2*Para.m_ini./((Kp+1).^3);
        ck2_eq(1, :) = 2*Para.m_ini*(135*Kp.^2*Para.Pe^2+105*Kp.^2+32*Kp*Para.Pe^2+210*Kp+2*Para.Pe^2+105)./(105*(Kp+1).^4);
        ck1_eq(1, :) = - Para.m_ini*(-1824*Kp.^4*Para.Pe^2-2100*Kp.^4+620*Kp.^3*Para.Pe^2-4200*Kp.^3+418*Kp.^2*Para.Pe^2-2100*Kp.^2+80*Kp*Para.Pe^2+6*Para.Pe^2)./(1575*(Kp+1).^5);
        
        for i = 1:num_set
            m2_eq(:, i) = ck1_eq(:, i)' * exp(-p(:,i).^2*Para.tD) ...
                        + ck2_eq(:, i)' * ( tD_mat .* exp(-p(:,i).^2*Para.tD)) ...
                        + ck3_eq(:, i)' * ( tD_mat.^2 .* exp(-p(:,i).^2*Para.tD));
        end

        M_eq(:, :, 3) = m2_eq./m0_eq - (m1_eq./m0_eq).^2;
    end
else
    
    error('Please input a method \n')
    
end