function [av_log_lik, L_t, Y_err_mat, X_tgt_mat, Y_fit_mat] = f_loglik_adwu(parms, data, nx, tau)


%--------------------------------------------------------------------------
% Initializing the parameters of the model
%--------------------------------------------------------------------------
K = length(tau);
T = size(data,1);
ny = size(data,2);
%{
% Using the parms from the 5-factor model and some modifications
mu = zeros(nx,1); mu(1)= parms(1);

Phi = diag(parms(2:5));
Phi(2:4,1) = parms(6:8);
Phi(3,2) = parms(9);
Phi(4,3) = parms(10);
Phi = Phi + eye(4);

S_0tilde = diag(parms(11:14));
S_0mat = S_0tilde*S_0tilde';
S_0 = S_0mat(:);

S_1tilde = diag(parms(15:17));
S_1tilde(1,2:end) = parms(18:19);
S_1tilde(2,3) = parms(20);
S_1temp = S_1tilde*S_1tilde';
S_1mat = zeros(nx); S_1mat(1:nx-1,1:nx-1) = S_1temp;
S_1 = zeros(nx^2,nx); S_1(:,end) = S_1mat(:);

del_0 = 0;

del_1 = zeros(nx,1); del_1(1,1) = parms(21); del_1(3,1) = parms(22);

lambda_0 = zeros(nx,1);

lambda_1 = diag([parms(23) 0  parms(24) 0]);
lambda_1(1,2:4) = parms(25:27);
lambda_1(2,4) = parms(28);
lambda_1(3,4) = parms(29);
%}

% Using the parms from the 8-factor model and some modifications
% This model has pretty much no restrictions on the parameters
%{
mu = zeros(nx,1); mu(1)= parms(1);

Phi = diag(parms(2:5));
Phi(2:4,1) = parms(6:8);
Phi(3,2) = parms(9);
Phi(4,3) = parms(10);
Phi(1,2:3)=parms(11:12);
Phi = Phi + eye(4);
% TEST2: TAKE OUT Phi(1,2:3) to make Phi lower-triangular
Phi(1,2:3) = [0 0];

S_0tilde = diag(parms(13:16));
S_0mat = S_0tilde*S_0tilde';
S_0 = S_0mat(:);

%%% S_1tilde not in v1
S_1tilde = diag(parms(15+2:17+2));
S_1tilde(1,2:end) = parms(18+2:19+2);
S_1tilde(2,3) = parms(20+2);
S_1temp = S_1tilde*S_1tilde';
S_1mat = zeros(nx); S_1mat(1:nx-1,1:nx-1) = S_1temp;
S_1 = zeros(nx^2,nx); S_1(:,end) = S_1mat(:);

del_0 = 0;

del_1 = zeros(nx,1); del_1(1:3) = parms(23:25);

lambda_0 = zeros(nx,1);
% TEST3: lambda_0 non-zero in v1, so make lambda_0 non-zero
%lambda_0(1:3) = parms(35:37);

lambda_1 = zeros(nx);
lambda_1(1,1:3) = parms(26:28);
lambda_1(2,1:3) = parms(29:31);
lambda_1(3,1:3) = parms(32:34);
% TEST4: lambda_1 upper triangular in v2, make give it the values from v2 setup1
lambda_1(2:3,1) = [0 0]';
lambda_1(3,2) = 0;
% TEST5: lambda_1's right-most column in v1 is 0. Change this to non-zero
% as in v2 setup1
lambda_1(:,4) = [-0.0025 -0.0065 0.0066 0]';
%}
getting_good_parms
%--------------------------------------------------------------------------
% From the older version of the paper
R0=0.001;   
R1=zeros(1,16)+ 0.00045905;
R2=zeros(1,32)+ 0.05e-006;
nr=0.1e-005;
R2([4:4:32]) = nr;
%%
%--------------------------------------------------------------------------
% Solving the term-structure model for coefficients of the measurement
% equation. The coefficients follow a recursive process.
%--------------------------------------------------------------------------
N = 520;            % Max. maturity = 10yrs; # of weeks in a yr. = 52
gamma = 5200;       % Scaling factor for the quadratic term when the yield 
                    % data are weekly and expressed as % per annum
phi = [1 0 0 0]';

% This initialization structure is necessary because it will help in 
% stacking across maturities:
C_0r_temp = zeros(N,1); C_1r_temp = cell(N,1);
C_0r_temp(1) = del_0; C_1r_temp{1} = del_1;
C_0y_temp = zeros(N,1); C_1y_temp = cell(N,1);
C_0y_temp(1) = (mu - lambda_0)'*phi - (1/gamma)*(1/2)*kron(phi,phi)'*S_0 + del_0;
C_1y_temp{1} = (Phi - lambda_1)'*phi - (1/gamma)*(1/2)*(kron(phi,phi)'*S_1)' + del_1;

for i=2:N
    C_1r_temp{i} = (Phi - lambda_1)'*C_1r_temp{i-1} - (1/gamma)*(1/2)*(kron(C_1r_temp{i-1},C_1r_temp{i-1})'*S_1)' + del_1;
    C_0r_temp(i) = C_0r_temp(i-1) + (mu - lambda_0)'*C_1r_temp{i-1} - (1/gamma)*(1/2)*(kron(C_1r_temp{i-1},C_1r_temp{i-1})'*S_0) + del_0;
    
    C_1y_temp{i} = (Phi - lambda_1)'*(C_1y_temp{i-1} + phi) - (1/gamma)*(1/2)*(kron(C_1y_temp{i-1} + phi,C_1y_temp{i-1} + phi)'*S_1)' + del_1;
    C_0y_temp(i) = C_0y_temp(i-1) + (mu - lambda_0)'*(C_1y_temp{i-1} + phi) - (1/gamma)*(1/2)*(kron(C_1y_temp{i-1} + phi,C_1y_temp{i-1} + phi)'*S_0) + del_0;
end
% Matching the short end of the yield curve using the 3-month T-Bill
% The thirteenth week gives the 3-months maturity yield
C_0y_3m = C_0y_temp(13);   
C_1y_3m = C_1y_temp{13};
c_0y_3m = C_0y_3m/(13);
c_1y_3m = C_1y_3m/(13);

C_0r = zeros(K,1); C_1r = cell(K,1);
C_0y = zeros(K,1); C_1y = cell(K,1);
C_0t = cell(K,1); C_1t = cell(K,1);
c_0t = cell(K,1); c_1t = cell(K,1);
d_0t = cell(K,1); d_1t = cell(K,1);
for i=1:K
    t = i+2;
    w = t*52;
    C_0r(i) = C_0r_temp(w);         C_1r{i} = C_1r_temp{w};
    C_0y(i) = C_0y_temp(w);         C_1y{i} = C_1y_temp{w};
    C_0t{i} = [C_0y(i); C_0r(i)];   C_1t{i} = [C_1y{i}';C_1r{i}'];
    c_0t{i} = C_0t{i}/w;            c_1t{i} = C_1t{i}/w;
    
    d_0t{i} = kron(c_1t{i},c_1t{i})*S_0;
    d_1t{i} = kron(c_1t{i},c_1t{i})*S_1;
end

C_0 = cell2mat([c_0y_3m; c_0t; d_0t]);
C_1 = cell2mat([c_1y_3m'; c_1t; d_1t]);


%%
%--------------------------------------------------------------------------
% Using the Unscented Kalman Filter to get a smoothed estimate for the 
% state variables. Due to the state dependence of the variance of the state
% variable, the assumption of a linear process do not hold. Hence, a 
% standard Kalman Filter cannot be used.
%--------------------------------------------------------------------------
X_init = [0 0 0 0.01]';
P_init = 1/52*eye(nx);
R = diag([R0 R1 R2]);
Y = data;

X_tgt_mat = nan(T,nx);
Y_err_mat = nan(T,ny);
Y_fit_mat = nan(T,ny);
L_t = nan(T,1);

for t=1:T
    Y_t = Y(t,:)';
    [X_tgt, P_tgt, Y_err, omega_y, K_t] = f_ukf(Y_t, X_init, P_init, R, C_0, C_1, mu, Phi, S_0, S_1);

    if isfinite(det(P_tgt)) && isfinite(det(omega_y)) && isfinite(1/det(omega_y))
        L_t(t) = -0.5*(log(abs(det(omega_y))) + Y_err'*inv(omega_y)*Y_err);
        Y_err_mat(t,:) = Y_err';
        X_tgt_mat(t,:) = X_tgt';
        Y_fit_mat(t,:) = (C_0 + C_1*X_tgt)';
    else
        L_t = -1e10;
        disp('bad parameter');
        break; return;
    end
    X_init = X_tgt;
    P_init = P_tgt;
end

av_log_lik = -mean(L_t)/100;

return