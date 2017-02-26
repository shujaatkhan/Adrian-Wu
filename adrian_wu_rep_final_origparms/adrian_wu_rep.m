clear all
close all
warning off all;
clc
workdir = pwd;
% mkdir('data')
rundate = clock;
%--------------------------------------------------------------------------
nx = 5;            % # of latent factors
tau = [3:10];      % List of maturities
K = length(tau);   % # of maturities
plots_on = 0;      % 1 to make plots; 0 otherwise
fignum = 0;        % Initializing the figure number
% daily2weekly;      % Convert daily to weekly data and estimate BEKK GARCH
%%
data = csvread([workdir,'/data/data_all_2017_2_8.csv']);
data_vix = csvread([workdir,'/data/data_vix_2017_2_8.csv']);
data(:,1) = x2mdate(data(:,1));
data_vix(:,1) = x2mdate(data_vix(:,1));
data_inp = data(:,2:end);
mdates = data(:,1);
yields = data(:,2:18);
% mdates_yields = csvread([workdir,'/data/mdates_yields_2017_2_5.csv']);
% mdates = mdates_yields(:,1); yields = mdates_yields(:,2:end);
%  data = xlsread('data_all_0310.xls','data');
%  data_y = xlsread('data_all_0310.xls','mdates_yields');
%  mdates = data_y(:,1); yields = data_y(:,2:end);
%--------------------------------------------------------------------------
% Plotting the 10-yr nominal yield, real yield, and breakeven inflation
%--------------------------------------------------------------------------
if plots_on == 1
    fignum = fignum+1; figure(fignum)
    hold on
    plot(mdates,yields(:,16),'LineWidth',1)
    plot(mdates,yields(:,17),'--','LineWidth',1)
    plot(mdates,yields(:,16)-yields(:,17),'LineWidth',1)
    hold off
    datetick('x','YYYY');
    ylabel('percent annual')
    legend('Nominal 10-yr', 'Real 10-yr', 'Breakeven 10-yr', ...
        'Location', 'southwest')
    print(['Fig',num2str(fignum)],'-dpdf')
end
%--------------------------------------------------------------------------
% Plotting the 10-yr nominal variance, real variance, and covariance
%--------------------------------------------------------------------------
if plots_on == 1
    fignum = fignum+1; figure(fignum);
    hold on
    plot(mdates,52*data(:,end-3),'LineWidth',1)
    plot(mdates,52*data(:,end-2),'--','LineWidth',1)
    plot(mdates,52*data(:,end),'LineWidth',1)
    hold off
    datetick('x','YYYY');
    ylabel('percent annual')
    legend('Nominal variance 10-yr', 'Covaraince 10-yr', 'Real variance 10-yr', ...
        'Location', 'northwest')
    print(['Fig',num2str(fignum)],'-dpdf')
end
%--------------------------------------------------------------------------
% Creating optimzation options structure
%--------------------------------------------------------------------------
tol = 1e-4;
maxit = 5e5;
fminopt = optimset; 
fminopt.Display = 'iter'; fminopt.MaxIter = maxit; fminopt.MaxFunEvals = maxit;
fminopt.TolX = tol; fminopt.TolFun = tol; fminopt.LargeScale = 'off';
%--------------------------------------------------------------------------
% Parameter estimation using maximum likelihood
%--------------------------------------------------------------------------
% 102615 - Using data from 2004 to 2009
% parms = csvread([workdir,'/data/parms_est_2015_10_26.csv']);
% parms = [0.006541231 -0.005557388 -0.004745913 -0.034383616 -0.004803461 -0.003535709 -0.000133338 -0.00403841 4.8543E-05 0.008569288 0.055464176 0.021321538 -0.001946605 -0.036521116 -0.001920892 -0.05673472 -0.004187146 0.098166737 -0.018555409 -0.754950909 -0.195037734 1.003236581 -0.065371609 0.055849867 -0.82769416 -0.197693299 0.171806735 -9.197346488 0.555022582 0.028456815 -0.002389565 0.018442464 0.005477629 0.00410086 -0.003792662 -0.005704861 0.007927754]';
% 102815 - Using data from 2003 to 2011
% parms = csvread([workdir,'/data/parms_est_2015_10_28.csv']);
% parms = [0.00503427514920085,-0.00721005429439332,-0.00378630694905654,-0.0341154863764201,-0.00500248281780726,-0.00286989738790494,-0.000118712251416265,-0.00281967993789441,-0.000204671569541372,0.00827650073422688,0.0732432322640664,0.0226448638558452,-0.00232053613932428,-0.0426898196313971,-0.00193742265437202,-0.0603726522066159,-0.00185349665326880,0.0810118785733868,-0.0241867748920772,-0.830121087039134,-0.227257167798093,0.950015489702371,-0.108475553014499,0.0564214310368789,-0.809012633299193,-0.282589494549213,0.159498437360001,-11.3760637165020,0.616547493493698,0.0370132885111431,-0.00160211207067715,0.0200637672892669,0.00528380057058203,0.00326025412627588,-0.00257124192494099,-0.00681688076196265,0.00776021576697711]';
% Original parameter estimates from Adrian and Wu on data from 03-09
% parms = [0.0072552,-0.0068665,-0.0038404,-0.033143,-0.00448,-0.0034581,-0.00012567,-0.0031061,-0.00042656,0.008186,0.071966,0.021114,-0.0024013,-0.041596,-0.0017874,-0.060856,-0.0025354,0.082803,-0.017733,-0.83531,-0.23082,0.97032,-0.096882,0.059836,-0.77935,-0.26843,0.16712,-12.833,0.56487,0.034964,-0.0018776,0.016675,0.0040866,0.0035249,-0.0025024,-0.0064531,0.0065505]';
% Using the original parameters from Adrian & Wu, we re-estimate the model
% and get the following parameter estimates
parms = csvread([workdir,'/data/parms_est_0103_12092017_2_8.csv']);
%--------------------------------------------------------------------------
% Uncomment for estimation of parameters
% Current run time on Macbook Pro 5668seconds
%--------------------------------------------------------------------------
% tic
% [parms_est, fval, exitflag, output] = fminsearch('f_loglik_adwu', parms, fminopt, data_inp, nx, tau);
% toc
% parms = parms_est;
% csvwrite(['./data/parms_est_0103_1209',num2str(rundate(1)),'_',num2str(rundate(2)),'_',num2str(rundate(3)),'.csv'],parms_est);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
tic
[av_log_lik, L_t, Y_err_mat, X_tgt_mat, Y_fit_mat] = f_loglik_adwu(parms, data_inp, nx, tau);
toc
display(-av_log_lik)
%--------------------------------------------------------------------------
% Reassigning values to the parameters
%--------------------------------------------------------------------------
mu = zeros(nx,1); mu(1)= parms(1); mu(2) = parms(2);

Phi = diag([parms(3),parms(4),parms(5),parms(6),parms(7)]);
Phi(2,1) = parms(8);
Phi(3,1) = parms(9);
Phi(4,1:3) = [parms(10) parms(11) parms(12)];
Phi(5,1:2) = [parms(13) parms(14)]; Phi(5,4) = parms(15);
Phi = Phi + eye(5);

S_0tilde = diag([parms(16) 0 parms(17) parms(18) parms(19)]);
S_0mat = S_0tilde*S_0tilde';
S_0 = S_0mat(:);

S_1tilde = diag([0 parms(20) parms(21) parms(22)]);
S_1tilde(1,2:4) = [parms(23) parms(24) parms(25)];
S_1tilde(2,4) = parms(26);
S_1tilde(3,4) = parms(27);
S_1temp = S_1tilde*S_1tilde';
S_1mat = zeros(nx); S_1mat(1:nx-1,1:nx-1) = S_1temp;
S_1 = zeros(nx^2,nx); S_1(:,end) = S_1mat(:);

del_0 = 0;

del_1 = zeros(nx,1); del_1(2,1) = parms(28); del_1(4,1) = parms(29);

lambda_0 = zeros(nx,1);

lambda_1 = diag([0 parms(30) 0  parms(31) 0]);
lambda_1(1,2) = parms(32);
lambda_1(2,3:5) = [parms(33) parms(34) parms(35)];
lambda_1(3,5) = parms(36);
lambda_1(4,5) = [parms(37)];
%--------------------------------------------------------------------------
% Plotting the factors
%--------------------------------------------------------------------------
if plots_on==1
    fignum = fignum+1; figure(fignum);
    hold on
    plot(mdates,X_tgt_mat(:,1)+X_tgt_mat(:,2),'LineWidth',1)
    plot(mdates,X_tgt_mat(:,3)+X_tgt_mat(:,4),'LineWidth',1)
    hold off
    datetick('x','YYYY');
    ylabel('percent annual')
    legend('\pi', 'm', 'Location', 'northwest')
    print(['Fig',num2str(fignum)],'-dpdf')
    
    fignum = fignum+1; figure(fignum);
    plot(mdates,X_tgt_mat(:,5)*52,'LineWidth',1);
    datetick('x','YYYY');
    ylabel('percent annual')
    legend('\sigma^2', 'Location', 'northwest')
    print(['Fig',num2str(fignum)],'-dpdf')
end
%--------------------------------------------------------------------------
% Calculating Expected Inflation
%--------------------------------------------------------------------------
N = 520;            % Max. maturity = 10yrs; # of weeks in a yr. = 52
phi = [1 1 0 0 0];
I5 = eye(5);
C_0tilde = cell(N,1); C_0tilde{1} = 0 + (I5 + 0)*mu;
C_1tilde = cell(N,1); C_1tilde{1} = (I5 + 0)*Phi;
for t=2:N
    C_0tilde{t} = C_0tilde{t-1} + (I5 + C_1tilde{t-1})*mu;
    C_1tilde{t} = (I5 + C_1tilde{t-1})*Phi;
end

c_0tilde_tau = cell(K,1);
c_1tilde_tau = cell(K,1);
inf_exp = cell(1,K);
for i=1:K
    t = i+2;
    w = t*52;
    c_0tilde_tau{i} = C_0tilde{w}/w;
    c_1tilde_tau{i} = C_1tilde{w}/w;
    inf_exp{i} = (phi*(repmat(c_0tilde_tau{i},1,size(X_tgt_mat,1)) + c_1tilde_tau{i}*X_tgt_mat'))';
end
inf_exp = cell2mat(inf_exp);
%--------------------------------------------------------------------------
% Plot to compare 10-yr breakeven inflation with model's estimate for
% 10-yr expected inflation
%--------------------------------------------------------------------------
if plots_on==1
    fignum=fignum+1; figure(fignum);
    be10 = data(:,17)-data(:,18);
    be10_fit = Y_fit_mat(:,16)-Y_fit_mat(:,17);
    be10_err = be10-be10_fit;
    riskprem_10 = be10 - inf_exp(:,8);
    riskprem_10_fit = be10_fit - inf_exp(:,8);
    hold on
    plot(mdates, be10,'LineWidth',1)
    plot(mdates, inf_exp(:,8),'LineWidth',1)
    %plot(mdates, riskprem_10,'LineWidth',1)
    plot(mdates, riskprem_10_fit,'LineWidth',1)
    plot(mdates, be10_err,'LineWidth',1)
    hold off
    datetick('x','YYYY');
    ylabel('percent annual')
    %legend('Breakeven 10-yr', 'Expected inflation 10-yr', 'Risk premium 10-yr', ...
    %    'Risk premium 10-yr - FIT', 'Breakeven 10-yr fitting error', 'Location', 'southwest')
    legend('Breakeven 10-yr', 'Expected inflation 10-yr', ...
        'Risk premium 10-yr - FIT', 'Breakeven 10-yr fitting error', 'Location', 'southwest')
    print(['Fig',num2str(fignum)],'-dpdf')
end
%--------------------------------------------------------------------------
% Plot forward 5-10 yr
%--------------------------------------------------------------------------
if plots_on==1
    fignum=fignum+1; figure(fignum);
    be5 = data(:,7)-data(:,8);
    be5_fit = Y_fit_mat(:,6)-Y_fit_mat(:,7);
    riskprem_5 = be5 - inf_exp(:,3);
    riskprem_5_fit = be5_fit - inf_exp(:,3);
    be_fwd_5_10 = 2*be10 - be5;
    inf_exp_fwd_5_10 = 2*inf_exp(:,end) - inf_exp(:,3);
    riskprem_fwd_5_10 = 2*riskprem_10 - riskprem_5;
    riskprem_fwd_5_10_fit = 2*riskprem_10_fit - riskprem_5_fit;
    hold on
    plot(mdates, be_fwd_5_10,'LineWidth',1)
    plot(mdates, inf_exp_fwd_5_10,'LineWidth',1)
    %plot(mdates, riskprem_fwd_5_10,'LineWidth',1)
    plot(mdates, riskprem_fwd_5_10_fit,'LineWidth',1)
    hold off
    datetick('x','YYYY');
    ylabel('percent annual')
    %legend('Forward Breakeven 5-10 yr', 'Forward Expected inflation 5-10 yr', 'Forward Risk premium 5-10 yr', ...
    %    'Forward Risk premium 5-10 yr - FIT', 'Location', 'southwest')
    legend('Forward Breakeven 5-10 yr', 'Forward Expected inflation 5-10 yr', ...
        'Forward Risk premium 5-10 yr - FIT', 'Location', 'southwest')
    print(['Fig',num2str(fignum)],'-dpdf')
end
%--------------------------------------------------------------------------
% Plot forward 5-10 yr risk and VIX
%--------------------------------------------------------------------------
if plots_on==1
    fignum=fignum+1; figure(fignum);
    yyaxis left
    hold on
    %plot(mdates, riskprem_fwd_5_10,'LineWidth',1)
    plot(mdates, riskprem_fwd_5_10_fit,'LineWidth',1)
    ylabel('percent annual')
    yyaxis right
    plot(mdates, data_vix(:,2),'LineWidth',1)
    hold off
    datetick('x','YYYY');
    ylabel('percent annual')
    %legend('Forward Risk premium 5-10 yr', 'Forward Risk premium 5-10 yr - FIT', 'VIX', ...
    %    'Location', 'southwest')
    legend('Forward Risk premium 5-10 yr - FIT', 'VIX index (right axis)', 'Location', 'northwest')
    print(['Fig',num2str(fignum)],'-dpdf')
end



