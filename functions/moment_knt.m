function M_knt = moment_knt(moment, method, Ka, Kd, Para, num_k)
% MOMENT_KNT normalized global moments for kinetic model
%   M_knt = moment_knt(moment, method, Ka, Kd, Para) calculates the zeroth,
%   first, second moments for transient transport between two plates with
%   Poiseuille flow and kinetic adsorption and desorption on the surface
%   using numerical inverse Laplace transform (num), or truncated series
%   solution derived by residue theorem (series).
%
% INPUT:
%   moment -- 0 or 1 or 2, assigning up to which moment you want calculate.
%   method -- 'num' to use the numerical inverse Laplace transform, or
%             'series' to use the truncated series solution.
%   Ka -- 1*num_set array of adsorption rate constants.
%   Kd -- 1*num_set array of desorption rate constants.
%   Para -- structure including the physical and computational parameters
%   num_k -- number of roots used for series solution
%
% OUTPUT:
%   M_knt -- num_t * num_set * moment matrix for the moments
%
% OTHER PARAMETERS:
%   num_sum -- number of terms used to sum up in talbot_inversion
%

if(nargin < 6)
    num_k = 2;
end

num_set = length(Ka); 
num_t = length(Para.tD);
m0_knt = zeros(num_t, num_set);
m1_knt = zeros(num_t, num_set);
m2_knt = zeros(num_t, num_set);
M_knt = zeros(num_t, num_set, moment + 1);

Ka_mat = repmat(Ka, num_k, 1);
Kd_mat = repmat(Kd, num_k, 1);
tD_mat = repmat(Para.tD, num_k, 1);

% % If mupad has not been evaluated, evaluate to the end
% Eval_Mark = getVar(nb, 'Eval_Mark');
% if ~strcmp(char(Eval_Mark), 'yes');
%     evaluateMuPADNotebook(nb);
% end

% ------ Numerical Inverse Laplace Transform ------
num_sum = 20;

if strcmpi(method,'num') == 1
    
    FN_mupad = './mupad/moment_knt.mn';
    nb = mupad(FN_mupad);
    
    % zeroth moment
    if moment >= 0
        Lm0_ana = getVar(nb, 'Lm0_ana');
        Lm0_ana = matlabFunction(Lm0_ana);
        for i = 1:num_set
            fun = @(s) Lm0_ana(Ka(i), Kd(i), Para.m_ini, s);
            m0_knt(:, i) = talbot_inversion(fun, Para.tD, num_sum);
        end
        
        M_knt(:, :, 1) = m0_knt./Para.m_ini;
        
    end
    
    % first moment
    if moment >= 1
        Lm1_ana = getVar(nb, 'Lm1_ana');
        Lm1_ana = matlabFunction(Lm1_ana);
        for i = 1:num_set
            fun = @(s) Lm1_ana(Ka(i), Kd(i), Para.Pe, Para.m_ini, s);
            m1_knt(:, i) = talbot_inversion(fun, Para.tD, num_sum);
        end
        
        M_knt(:, :, 2) = m1_knt./m0_knt;
    end
    
    % second moment
    if moment >= 2
        Lm2_ana = getVar(nb, 'Lm2_ana');
        Lm2_ana = matlabFunction(Lm2_ana);
        for i = 1:num_set
            fun = @(s) Lm2_ana(Ka(i), Kd(i), Para.Pe, Para.m_ini, s, Para.var0);
            m2_knt(:, i) = talbot_inversion(fun, Para.tD, num_sum);
        end
        
        M_knt(:, :, 3) = m2_knt./m0_knt - (m1_knt./m0_knt).^2;
    end

% ------ Truncated Analytical Solution ------
elseif strcmpi(method, 'series') == 1

    p = root_tr_knt(Ka, Kd, num_k);
    
    if moment >= 0
        
        ak = (2*Ka_mat.^2*Para.m_ini)./(Ka_mat.^2.*p.^2 + Ka_mat.*Kd_mat + Ka_mat.*p.^2 + Kd_mat.^2 - 2*Kd_mat.*p.^2 + p.^4);
        ak(1, :) = Para.m_ini*Kd./(Ka + Kd);      % coefficients for zero root
        
        for i = 1:num_set
            m0_knt(:, i) = ak(:, i)' * exp(-p(:,i).^2*Para.tD);
        end
        
        M_knt(:, :, 1) = m0_knt./Para.m_ini;
    end
    
    if moment >= 1
        bk2 = (Ka_mat.^2*Para.Pe*Para.m_ini.*(4*Ka_mat.^2.*p.^4 + 3*Ka_mat.^2.*p.^2 - 3*Ka_mat.*Kd_mat + 3*Ka_mat.*p.^2 + 4*Kd_mat.^2.*p.^2 - 3*Kd_mat.^2 - 8*Kd_mat.*p.^4 + 6*Kd_mat.*p.^2 + 4*p.^6 - 3*p.^4))./...
                 (2*p.^2.*(Ka_mat.^2.*p.^2 + Ka_mat.*Kd_mat + Ka_mat.*p.^2 + Kd_mat.^2 - 2*Kd_mat.*p.^2 + p.^4).^2);
        bk1 = -(Ka_mat*Para.Pe*Para.m_ini.*(-6*Ka_mat.^5.*p.^6 + 15*Ka_mat.^5.*p.^4 + 16*Ka_mat.^4.*Kd_mat.*p.^4 + 51*Ka_mat.^4.*Kd_mat.*p.^2 - 24*Ka_mat.^4.*p.^6 + 30*Ka_mat.^4.*p.^4 - 8*Ka_mat.^3.*Kd_mat.^2.*p.^4 ...
                + 93*Ka_mat.^3.*Kd_mat.^2.*p.^2 + 30*Ka_mat.^3.*Kd_mat.^2 + 24*Ka_mat.^3.*Kd_mat.*p.^6 - 96*Ka_mat.^3.*Kd_mat.*p.^4 + 51*Ka_mat.^3.*Kd_mat.*p.^2 - 16*Ka_mat.^3.*p.^8 ...
                + 3*Ka_mat.^3.*p.^6 + 15*Ka_mat.^3.*p.^4 + 40*Ka_mat.^2.*Kd_mat.^3.*p.^2 + 84*Ka_mat.^2.*Kd_mat.^3 -128*Ka_mat.^2.*Kd_mat.^2.*p.^4 - 51*Ka_mat.^2.*Kd_mat.^2.*p.^2 + 136*Ka_mat.^2.*Kd_mat.*p.^6 ...
                - 54*Ka_mat.^2.*Kd_mat.*p.^4 - 48*Ka_mat.^2.*p.^8 + 21*Ka_mat.^2.*p.^6 - 2*Ka_mat.*Kd_mat.^4.*p.^2 + 78*Ka_mat.*Kd_mat.^4 + 16*Ka_mat.*Kd_mat.^3.*p.^4 - 222*Ka_mat.*Kd_mat.^3.*p.^2 ...
                - 36*Ka_mat.*Kd_mat.^2.*p.^6 + 198*Ka_mat.*Kd_mat.^2.*p.^4 + 32*Ka_mat.*Kd_mat.*p.^8 -42*Ka_mat.*Kd_mat.*p.^6 - 10*Ka_mat.*p.^10 - 12*Ka_mat.*p.^8 +24*Kd_mat.^5 - 120*Kd_mat.^4.*p.^2 ...
                + 240*Kd_mat.^3.*p.^4 - 240*Kd_mat.^2.*p.^6 + 120*Kd_mat.*p.^8 - 24*p.^10))./...
               (2*p.^4.*(Ka_mat.^2.*p.^2 + Ka_mat.*Kd_mat + Ka_mat.*p.^2 + Kd_mat.^2 - 2*Kd_mat.*p.^2 + p.^4).^3);
        bk2(1, :) = Para.Pe*Para.m_ini.*Kd.^2./(Ka + Kd).^2;
        bk1(1, :) = 2*Ka.*Kd.*Para.Pe*Para.m_ini.*(6*Ka+Kd+15)./(15*(Ka+Kd).^3);
        
        for i = 1:num_set
            m1_knt(:, i) = bk1(:, i)' * exp(-p(:,i).^2*Para.tD) + bk2(:, i)' * ( tD_mat .* exp(-p(:,i).^2*Para.tD));
        end
        
        M_knt(:, :, 2) = m1_knt./m0_knt;
    end
    
    if moment >= 2
        FN_mupad = './mupad/Res_limit_knt.mn';
        nb = mupad(FN_mupad);

        f_ck3 = getVar(nb, 'ck3');
        f_ck3 = matlabFunction(f_ck3);
        ck3 = f_ck3(Ka_mat, Kd_mat, Para.Pe, Para.m_ini, p);

        f_ck2 = getVar(nb, 'ck2');
        f_ck2 = matlabFunction(f_ck2);
        ck2 = f_ck2(Ka_mat, Kd_mat, Para.Pe, Para.m_ini, p);      

        f_ck1 = getVar(nb, 'ck1');
        f_ck1 = matlabFunction(f_ck1);
        ck1 = f_ck1(Ka_mat, Kd_mat, Para.Pe, Para.m_ini, p);
        
        ck3(1, :) = Kd.^3*Para.Pe^2*Para.m_ini./((Ka+Kd).^3);
        ck2(1, :) = 2*Kd.^2*Para.m_ini.*(135*Ka.^2*Para.Pe^2+105*Ka.^2+32*Ka.*Kd.*Para.Pe^2+210*Ka.*Kd+315*Ka*Para.Pe^2+2*Kd.^2*Para.Pe^2+105*Kd.^2)./(105*(Ka+Kd).^4);
        ck1(1, :) = Kd*Para.m_ini.*(1824*Ka.^4*Para.Pe^2+2100*Ka.^4-620*Ka.^3.*Kd*Para.Pe^2+4200*Ka.^3.*Kd+8100*Ka.^3*Para.Pe^2+6300*Ka.^3-418*Ka.^2.*Kd.^2*Para.Pe^2+2100*Ka.^2.*Kd.^2 ...
                  - 5220*Ka.^2.*Kd*Para.Pe^2 + 12600*Ka.^2.*Kd + 9450*Ka.^2*Para.Pe^2 - 80*Ka.*Kd.^3*Para.Pe^2 - 720*Ka.*Kd.^2*Para.Pe^2 + 6300*Ka.*Kd.^2 - 9450*Ka.*Kd*Para.Pe^2 - 6*Kd.^4*Para.Pe^2) ...
                  ./(1575*(Ka+Kd).^5);
        
        for i = 1:num_set
            m2_knt(:, i) = ck1(:, i)' * exp(-p(:,i).^2*Para.tD) ...
                         + ck2(:, i)' * ( tD_mat .* exp(-p(:,i).^2*Para.tD)) ...
                         + ck3(:, i)' * ( tD_mat.^2 .* exp(-p(:,i).^2*Para.tD));
        end

        M_knt(:, :, 3) = m2_knt./m0_knt - (m1_knt./m0_knt).^2;
    end    
    
else
    
    error('Given method is not correct\n')
    
end