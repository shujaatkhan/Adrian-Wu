
mu = [ 0.0077   0  0   0]';
     
Phi=[0.9990        0         0         0
    0.0031    0.9998         0         0
    0.0036    0.0006    0.9884         0
   -0.0200         0   -0.0018    0.9999];

S_0tilde=[0.0616         0         0         0
               0   -0.0037         0         0
               0         0    0.0049         0
               0         0         0   -0.0177];
S_0mat = S_0tilde*S_0tilde'; S_0 = S_0mat(:);
           
S_1tilde=[   -0.7353    0.0500   -0.5000
                   0   -0.2308    0.1670
                   0         0    0.9703];
S_1temp = S_1tilde*S_1tilde'; S_1mat = zeros(nx); S_1mat(1:nx-1,1:nx-1) = S_1temp; S_1 = zeros(nx^2,nx); S_1(:,end) = S_1mat(:);
               
del_0 = 0;

del_1 = [ 0.0807  1.2922  0.9952  0]';
           
lambda_0 = [  0.0084 -0.0039  -0.0015   0]';
               
lambda_1 = [   -0.0058   -0.0029    0.0048   -0.0025
                     0   -0.0020    0.0144   -0.0065
                     0         0    0.0028    0.0066
                     0         0         0         0];