clear all
close all
warning off all;
clc
workdir = pwd;
% mkdir('data')
rundate = clock;
s = {}; s1 = {};
s{1} = [1:2]; s1{1} = 'mu';
s{2} = [3:15]; s1{2} = 'Phi';
s{3} = [16:19]; s1{3} = 'S_0';
s{4} = [20:27]; s1{4} = 'S_1';
s{5} = [28:29]; s1{5} = 'del_1';
s{6} = [30:37]; s1{6} = 'lambda_1';

fignum = 0;
for p=1:6
    %--------------------------------------------------------------------------
    nx = 5;            % # of latent factors
    tau = [3:10];      % List of maturities
    K = length(tau);   % # of maturities
    plots_on = 1;      % 1 to make plots; 0 otherwise
    %daily2weekly;      % Convert daily to weekly data and estimate BEKK GARCH
    %%
    data = csvread([workdir,'/data/data_all_2017_2_7.csv']);
    data_vix = csvread([workdir,'/data/data_vix_2017_2_7.csv']);
    data(:,1) = x2mdate(data(:,1));
    data_vix(:,1) = x2mdate(data_vix(:,1));
    data_inp = data(:,2:end);
    mdates = data(:,1);
    yields = data(:,2:18);
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
    % Original parameter estimates from Adrian and Wu on data from 03-09
    parms_orig = [0.0072552,-0.0068665,-0.0038404,-0.033143,-0.00448,-0.0034581,-0.00012567,-0.0031061,-0.00042656,0.008186,0.071966,0.021114,-0.0024013,-0.041596,-0.0017874,-0.060856,-0.0025354,0.082803,-0.017733,-0.83531,-0.23082,0.97032,-0.096882,0.059836,-0.77935,-0.26843,0.16712,-12.833,0.56487,0.034964,-0.0018776,0.016675,0.0040866,0.0035249,-0.0025024,-0.0064531,0.0065505]';
    % Estimates with the large data set 2003-2016
    parms = csvread([workdir,'/data/parms_est_2017_2_7.csv']);   
    parms(s{p})=parms_orig(s{p});
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
        legend('Breakeven 10-yr', 'Expected inflation 10-yr', ...
            'Risk premium 10-yr - FIT', 'Breakeven 10-yr fitting error', 'Location', 'southwest')
        title([s1{p},' fixed to orig parms'])
        print(['Fig',num2str(fignum)],'-dpdf')
    end
    
end
