%% Moment Analysis
%
% This file computes the normalized global moments (zeroth, first, second)
% of a reactive solute in Poiseuille flow between plates with surface
% adsorption.
% 
% For more details, please refer to Li Zhang, Marc Hesse and Moran Wang
% (2017) Transient solute transport with sorption in Poiseuille flow.
% Journal of Fluid Mechanics, 828: 733-752.
%
% Input physical parameters : 
%       Pe -- Peclet number;
%       m_ini -- initial mass per volume
%       tD -- dimensionless time
%       var0 -- initial variance
%
% Surface adsorption can be
% (1) knt -- general kinetic model with linear adsorption and desorption.
% The moment can be obtained either using numerical inverse Laplace
% transform or truncated analytical solution. The truncated analytical
% solution is equivalent to using the result from chromatography by Khan to
% describe the late-time behavior.
% 
% Input:
%       Ka -- adsorption rate constant;
%       Kd -- desorption rate constant;
%
% ! Note that for method = 'num', 'moment_knt.mn' should be open and
% evaluated. For method = 'series' and num_k >=2, 'Res_limit.mn' should be
% open and evaluated.
%
% (2) Eq -- equilibrium or isotherm model. The moment can be obtained by
% numerical Laplace transform, or truncated analytical solution. The
% truncated analytical solution is equivalent to using the result from
% chromatography by Golay to describe the late-time behavior.
%
% Input:
%       Kp -- retention factor, or distribution coefficient
%
% ! Note that for method = 'num', 'moment_eq.mn' should be open and
% evaluated. For method = 'series' and num_k >=2, 'Res_limit_eq.mn' should be
% open and evaluated.
%
% (3) LM -- adsorption-only model, or Lungu & Moffatt (LM) model. We
% essentiall use the results from LM to get the decay rate for zeroth
% moment, slope for first moment and slope for second moment.
%
% Input:
%       Ka -- adsorption rate constant (sometimes called Gamma)
%
% Moments can be calculated by
% (1) method = 'num' : numerical inverse Laplace transform,
% (2) method = 'series' : truncated series solution derived by residue
% theorem.
%
% author : zl
% date : 2016/4/28
%

clear; close all; clc

addpath('./functions')

%% ============ Input Parameters ===============

% surface reaction parameters
Ka = 50; Kd = 1;
Kp = Ka./Kd;

% other dimensionless physical parameters
A = 1./(1 + Kp);
Para.m_ini = 1;
Para.Pe = 10;
Para.Da = Para.Pe./Kd;
Para.var0 = 0;
Para.t_diff = (3/sqrt(2)*Para.Pe)^(-2/3);
fprintf('Typical Diffusion Time Scale : %5.2f \n', Para.t_diff);

% time scale and time step
num_t = 1000;
min_t = 1e-3; max_t = 1e1;
dt = (max_t - min_t)/num_t;
Para.tD = logspace(log10(min_t), log10(max_t), num_t);

% moment level
moment = 2;

%% =============== Calculate model =================

%--- Kinetic Model ---
M_knt_num = moment_knt(moment, 'num', Ka, Kd, Para);
M_knt_series = moment_knt(moment, 'series', Ka, Kd, Para, 2);

%--- Equilibrium  Model ---
M_eq_num = moment_eq(moment, 'num', Kp, Para);
M_eq_series = moment_eq(moment, 'series', Kp, Para, 2);

%--- Lungu-Moffatt Model ---
[M_LM, B_lm, V_lm, D_lm] = moment_LM(moment, Ka, Para);

%% ==================  Plot  =====================
figure
set(groot, 'Units', 'centimeter')
scr_pos = get(groot, 'ScreenSize'); 
scr_w = scr_pos(3); scr_h = scr_pos(4);
fig_w = scr_w/1.5; fig_h = scr_h/2; fig_l = (scr_w - fig_w)/2; fig_b = (scr_h - fig_h)/2;
set(gcf, 'Units', 'centimeter', 'Position', [fig_l fig_b fig_w fig_h], 'PaperPositionMode', 'Auto')

skip = 1:10:num_t;

subplot 121
plot(Para.tD, M_knt_num(:,:,2), 'k', 'linewidth', 2); hold on
plot(Para.tD(skip), M_knt_series(skip,:,2), 'bs');
plot(Para.tD, M_eq_num(:,:,2), 'r', 'linewidth', 2);
plot(Para.tD(skip), M_eq_series(skip,:,2), 'go');
plot(Para.tD, M_LM(:,:,2), 'm', 'linewidth', 2);
xlabel('t'); ylabel('M_1');
legend('knt-num','knt-series','eq-num', 'eq-series', 'LM','location','best')
set(gca, 'fontsize', 14, 'fontname', 'times')

subplot 122
plot(Para.tD, M_knt_num(:,:,3), 'k', 'linewidth', 2); hold on
plot(Para.tD(skip), M_knt_series(skip,:,3), 'bs');
plot(Para.tD, M_eq_num(:,:,3), 'r', 'linewidth', 2);
plot(Para.tD(skip), M_eq_series(skip,:,3), 'go');
plot(Para.tD, M_LM(:,:,3), 'm', 'linewidth', 2);
xlabel('t'); ylabel('M_2')
legend('knt-num','knt-series','eq-num', 'eq-series', 'LM','location','best')
set(gca, 'fontsize', 14, 'fontname', 'times')