clear
clc
close all
format long g
% This is the folder where the VFI toolkit files are saved
folder1 = 'C:\Users\aledi\OneDrive\Desktop\vfi_toolkit\VFIToolkit-matlab';
addpath(genpath(folder1))

%% Set options
me_shock = 'tauchen'; % 'rouwenhorst' vs 'tauchen'

%% Set parameters

par.beta  = 0.96;
par.crra  = 2.0;
par.delta = 0.08;
par.alpha = 0.36;
% With progressive taxes
par.lam_hsv = 0.9;
par.tau_hsv = 0.1;
% With proportional taxes, tax_rate=1-lambda
%par.lam_hsv = 0.9;
%par.tau_hsv = 0;
% No taxes %TODO distrib gives error unless I turn off Tan improv
%par.lam_hsv = 1;
%par.tau_hsv = 0;

%% Guess for prices
% In the absence of idiosyncratic risk, the steady state equilibrium is given by
r_ss=1/par.beta-1;
K_ss=((r_ss+par.delta)/par.alpha)^(1/(par.alpha-1)); %The steady state capital in the absence of aggregate uncertainty.

%% Discretize productivity shock
n_z   = 7;
ave_z = 0.0;
rho_z = 0.6;
sig_z = sqrt(0.04*(1-rho_z^2)); 
switch me_shock
    case 'rouwenhorst'
        [z_grid_log,pi_z] = rouwenhorst(n_z,ave_z,rho_z,sig_z);
    case 'tauchen'
        Tauchen_q=3; 
        [z_grid_log,pi_z]=discretizeAR1_Tauchen(0,rho_z,sig_z,n_z,Tauchen_q);
    otherwise
        error('Selected me_shock not available')
end
z_grid = exp(z_grid_log);
aux = pi_z^1000;
z_distrib = aux(1,:);
Expectation_l = dot(z_grid,z_distrib);

%% Setup grids for a and d variables
n_d = 0;
d_grid = [];
n_a   = 1000;
a_min = 0;
a_max = 15*K_ss;
a_curve = 3;
a_grid = a_min+(a_max-a_min)*(linspace(0,1,n_a).^a_curve)';

%% Setup toolkit inputs
DiscountFactorParamNames={'beta'};

ReturnFn = @(aprime,a,z,K_to_L,alpha,delta,crra,lam_hsv,tau_hsv) Model_ReturnFn(aprime,a,z,K_to_L,alpha,delta,crra,lam_hsv,tau_hsv);

%% Solve general equil 

% Create functions to be evaluated
FnsToEvaluate.K = @(aprime,a,z) a; % A, assets or capital
FnsToEvaluate.L = @(aprime,a,z) z; % L, labor in efficiency units
GeneralEqmEqns.CapitalMarket = @(K_to_L,K,L) K_to_L-K/L; 

GEPriceParamNames={'K_to_L'};
% Set initial value for K/L
par.K_to_L = K_ss/Expectation_l;

% %% My stationary distribution function
% fprintf('Start distribution... \n')
% simoptions.verbose=1;
% tic
% StatDist=StationaryDist_Case1_ale(Policy,n_d,n_a,n_z,pi_z, simoptions);
% toc

%% Do general equilibrium
vfoptions.verbose   = 0;
vfoptions.lowmemory = 1;
simoptions.verbose  = 0;
heteroagentoptions.verbose=1; % verbose means that you want it to give you feedback on what is going on
heteroagentoptions.toleranceGEprices=10^(-6); % Accuracy of general eqm prices
heteroagentoptions.toleranceGEcondns=10^(-6); % Accuracy of general eqm eqns
%heteroagentoptions.fminalgo=0;
fprintf('Calculating price corresponding to the stationary general eqm \n')
[p_eqm,~,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(n_d,n_a,n_z,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,par,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions, vfoptions);
disp(GeneralEqmCondn)

%% Recompute at equil prices
par.K_to_L = p_eqm.K_to_L; % Put the equilibrium interest rate into Params so we can use it to calculate things based on equilibrium parameters

% par.r = par.alpha*par.K_to_L^(par.alpha-1)-par.delta;
% par.w = (1-par.alpha)*par.K_to_L^par.alpha;

[par.r,par.w] = fun_prices(par.K_to_L,par.alpha,par.delta);

% VFI
vfoptions.verbose = 1;
fprintf('Start VFI... \n')
[V1,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, par, DiscountFactorParamNames, [], vfoptions);
% Stationary distribution
fprintf('Start distribution... \n')
simoptions.verbose=1;
simoptions.tanimprovement=1; 
StatDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);

%% Calculate distributional statistics
fprintf('Distributional statistics... \n')

% Functions to be evaluated
FnsToEvaluate.K = @(aprime,a,z) a; 
FnsToEvaluate.Kprime = @(aprime,a,z) aprime;
FnsToEvaluate.L = @(aprime,a,z) z;
%FnsToEvaluate.H = @(d,aprime,a,z) d;
FnsToEvaluate.cons = @(aprime,a,z,K_to_L,alpha,delta,lam_hsv,tau_hsv) Model_ConsFn(aprime,a,z,K_to_L,alpha,delta,lam_hsv,tau_hsv);
%FnsToEvaluate.hours    = @(d,aprime,a,z,w) d;
FnsToEvaluate.earnings = @(aprime,a,z,w) w*z;
FnsToEvaluate.income   = @(aprime,a,z,r,w) w*z+r*a;
FnsToEvaluate.wealth   = @(aprime,a,z,r,w) a;
FnsToEvaluate.taxes   = @(aprime,a,z,K_to_L,alpha,delta,lam_hsv,tau_hsv) Model_TaxesFn(aprime,a,z,K_to_L,alpha,delta,lam_hsv,tau_hsv);

AllStats = EvalFnOnAgentDist_AllStats_Case1(StatDist, Policy, FnsToEvaluate, par, [], n_d, n_a, n_z, d_grid, a_grid, z_grid,simoptions);

% Aggregate moments
agg.KK = AllStats.K.Mean;
agg.LL = AllStats.L.Mean;
agg.CC = AllStats.cons.Mean;
agg.YY = (agg.KK^par.alpha)*(agg.LL^(1-par.alpha));
agg.II = par.delta*agg.KK;
agg.TT = AllStats.taxes.Mean;

resid_walras = abs(agg.CC+agg.II+agg.TT-agg.YY);

% Calculate the wealth gini
gini.wealth   = AllStats.wealth.Gini;
%gini_hours    = AllStats.hours.Gini;
gini.income   = AllStats.income.Gini;
gini.earnings = AllStats.earnings.Gini;
gini.cons     = AllStats.cons.Gini;

% Coefficient of variation for hours, earnings, income and wealth
%mom.cv.hours    = AllStats.H.StdDeviation/AllStats.H.Mean;
cv.earnings = AllStats.earnings.StdDeviation/AllStats.earnings.Mean;
cv.income   = AllStats.income.StdDeviation/AllStats.income.Mean;
cv.wealth   = AllStats.wealth.StdDeviation/AllStats.wealth.Mean;
cv.cons     = AllStats.cons.StdDeviation/AllStats.cons.Mean;

%% Display results

disp('==================================================================')
disp('PARAMETERS')
fprintf('beta  : %f \n',par.beta)
fprintf('crra  : %f \n',par.crra)
fprintf('alpha : %f \n',par.alpha)
fprintf('delta : %f \n',par.delta)
fprintf('rho_z : %f \n',rho_z)
fprintf('sig_z : %f \n',sig_z)
fprintf('lambda : %f \n',par.lam_hsv)
fprintf('tau    : %f \n',par.tau_hsv)
disp('------------------------------------')
disp('AGGREGATE QUANTITIES AND PRICES')
% fprintf('Corr(h,z)  : %f \n',mom.corr_h_z)
% fprintf('CV(h)      : %f \n',mom.cv.hours)
% fprintf('Hours      : %f \n',agg.HH)
fprintf('K/Y        : %f \n',agg.KK/agg.YY)
fprintf('Int. rate  : %f \n',par.r)
fprintf('Wage       : %f \n',par.w)
fprintf('C/Y        : %f \n',agg.CC/agg.YY)
fprintf('I/Y        : %f \n',agg.II/agg.YY)
fprintf('res GE     : %f \n',GeneralEqmCondn)
fprintf('res_walras : %f \n',resid_walras)

% fprintf('w*L/Y      : %f \n',p_eq.w*agg.LL/agg.YY)
% fprintf('I/Y        : %f \n',agg.II/agg.YY)
disp('------------------------------------')
disp('CV')
%fprintf('CV(Hours)   : %f \n',mom.cv.hours)
fprintf('CV(Earnings): %f \n',cv.earnings)
fprintf('CV(Income)  : %f \n',cv.income)
fprintf('CV(Consum)  : %f \n',cv.cons)
fprintf('CV(Wealth)  : %f \n',cv.wealth)
disp('------------------------------------')
disp('GINI')
%fprintf('Gini(Hours)   : %f \n',mom.gini.hours)
fprintf('Gini(Earnings): %f \n',gini.earnings)
fprintf('Gini(Income)  : %f \n',gini.income)
fprintf('Gini(Consum)  : %f \n',gini.cons)
fprintf('Gini(Wealth)  : %f \n',gini.wealth)
disp('------------------------------------')
% disp('CORR')
% fprintf('corr(Hours,z) : %f \n',mom.corr_h_z)
% fprintf('corr(Wealth,z): %f \n',mom.corr_k_z)
% disp('------------------------------------')
% disp('SHARES EARNINGS')
% fprintf('q1 earnings: %f \n',mom.shares.earnings(1))
% fprintf('q2 earnings: %f \n',mom.shares.earnings(2))
% fprintf('q3 earnings: %f \n',mom.shares.earnings(3))
% fprintf('q4 earnings: %f \n',mom.shares.earnings(4))
% fprintf('q5 earnings: %f \n',mom.shares.earnings(5))
% disp('------------------------------------')
% disp('SHARES WEALTH')
% fprintf('q1 wealth: %f \n',mom.shares.wealth(1))
% fprintf('q2 wealth: %f \n',mom.shares.wealth(2))
% fprintf('q3 wealth: %f \n',mom.shares.wealth(3))
% fprintf('q4 wealth: %f \n',mom.shares.wealth(4))
% fprintf('q5 wealth: %f \n',mom.shares.wealth(5))

%% Plots

Policy1 = squeeze(gather(Policy));

pol_ap = a_grid(Policy1);

figure
plot(a_grid,a_grid,'--','LineWidth',2)
hold on
plot(a_grid,pol_ap(:,1),'LineWidth',2)
hold on
plot(a_grid,pol_ap(:,n_z),'LineWidth',2)
legend('45 line', 'low z', 'high z')

kkk = 300;
figure
plot(a_grid(1:kkk),a_grid(1:kkk),'--','LineWidth',2)
hold on
plot(a_grid(1:kkk),pol_ap(1:kkk,1),'LineWidth',2)
hold on
plot(a_grid(1:kkk),pol_ap(1:kkk,n_z),'LineWidth',2)
legend('45 line', 'low z', 'high z')

mu_a = sum(StatDist,2);

kkk = 500;
figure
plot(a_grid(1:kkk),mu_a(1:kkk),'LineWidth',2)

figure
plot(a_grid,cumsum(mu_a),'LineWidth',2)


