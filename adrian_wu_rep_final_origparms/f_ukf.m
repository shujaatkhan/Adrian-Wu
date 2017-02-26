function [X_tgt, P_tgt, Y_err, omega_y, K_t] = f_ukf(Y_t, X_init, P_init, R, C_0, C_1, mu, Phi, S_0, S_1)
%--------------------------------------------------------------------------
% The Unscented Kalman Filter (UKF)
% The algorithm is based on the original paper on UKF:
% Julier, Simon J. and Jeffrey K. Uhlmann (1997), 'A New Extension of the
% Kalman Filter to Nonlinear Systems.'
%--------------------------------------------------------------------------
% Inputs:
% Y_t      = Vector of observation data
% X_init   = X_t-1|t-1
% P_init   = P_t-1|t-1
% R        = Pricing error's constant, diagonal covariance matrix R
% C_0, C_1 = Estimated coefficients for the term structure model (Eq. 21a)
% mu, Phi  = Estimated coefficients for the term structure model (Eq. 21b)
% S_0, S_1 = Estimated coefficients for the term structure model (Eq. 21c)
%
% Outputs:
% X_tgt    = X_t|t
% P_tgt    = P_t|t
% Y_err    = Y_t - Y_t|t-1
% omega_y  = Omega_t^y = C_1 Var(X_t|Y_t-1) C_1 + R
% K_t      = Kalman gain
%--------------------------------------------------------------------------
% Last edited: 10/28/2015
%--------------------------------------------------------------------------
nx = length(X_init);
ny = length(Y_t);
nsig = 2*nx + 1;
kappa = 0;          % Typically set to 0

% X_init = X_t-1|t-1
P_S = S_0 + S_1*X_init; P_S = reshape(P_S,5,5);
P_Q = P_init + P_S;

[P_Q_sqrt, p]=chol(P_Q);
% if not positive definite
if (p~=0)||(~isfinite(det(P_Q_sqrt)))
    P_Q_sqrt=.01*eye(nx);
end
%--------------------------------------------------------------------------
% Sigma points for X (Chi_tm1|tm1)
X_SigPts = repmat(X_init,1,nsig) + ...
           [zeros(nx,1) sqrt(nx+kappa)*P_Q_sqrt' -sqrt(nx+kappa)*P_Q_sqrt'];

% Weights for the sigma points
W_SigPts = 1/(nx+kappa) * [kappa 0.5*ones(1,nsig-1)];

% Chi_t|tm1
X_SigPts_tgtm1 = repmat(mu,1,nsig) + Phi*X_SigPts;
% ScY_t|tm1
Y_SigPts_tgtm1 = repmat(C_0,1,nsig) + C_1*X_SigPts_tgtm1;

% Predicted state X_t|t-1
X_tgtm1 = X_SigPts_tgtm1*W_SigPts';
% Predicted covariance of state P_t|t-1
W_XSigPts_mat = repmat(W_SigPts,nx,1);
diffX = X_SigPts_tgtm1 - repmat(X_tgtm1,1,nsig);
P_tgtm1 = (W_XSigPts_mat.*diffX)*diffX';

% Predicted observation Y_t|t-1
Y_tgtm1 = Y_SigPts_tgtm1*W_SigPts';
% Predicted covariance of observation Py_t|t-1
W_YSigPts_mat = repmat(W_SigPts,ny,1);
diffY = Y_SigPts_tgtm1 - repmat(Y_tgtm1,1,nsig);
Py_tgtm1_tmp = (W_YSigPts_mat.*diffY)*diffY';

% Predicted state-measurement cross-covariance Pxy_t|t-1
Pxy_tgtm1 = (W_XSigPts_mat.*diffX) * diffY';

Py_tgtm1 = Py_tgtm1_tmp + R;
[Py_tgtm1_chol,p] = chol(Py_tgtm1);
% if not positive definite
if (p~=0)||(~isfinite(det(Py_tgtm1)))
    Py_tgtm1 = 0.01*eye(ny);
end
%--------------------------------------------------------------------------
% Kalman Gain
K_t = Pxy_tgtm1*pinv(Py_tgtm1);

% Y_t - Y_t|t-1
Y_err = Y_t - Y_tgtm1;

% X_t|t
X_tgt = X_tgtm1 + K_t*Y_err;

% P_t|t
P_tgt = P_tgtm1 - K_t*Py_tgtm1*K_t';

% omega_y
omega_y = Py_tgtm1;
%--------------------------------------------------------------------------

return