tau_vect=[1:520]';
C_0y_temp_tau = C_0y_temp./tau_vect;
C_0r_temp_tau = C_0r_temp./tau_vect;

C_1y_temp_mat = C_1y_temp';
C_1y_temp_mat = cell2mat(C_1y_temp_mat);
tau_vec = [1:520];
tau_mat = repmat(tau_vec,5,1);
C_1y_temp_mat_tau = C_1y_temp_mat./tau_mat;
C_1r_temp_mat = C_1r_temp';
C_1r_temp_mat = cell2mat(C_1r_temp_mat);
C_1r_temp_mat_tau = C_1r_temp_mat./tau_mat;

for tau=1:N
    C_0t{tau} = [C_0y_temp(tau); C_0r_temp(tau)];
    C_1t{tau} = [C_1y_temp{tau}';C_1r_temp{tau}'];
    
    D_0t{tau} = (kron(C_1t{tau},C_1t{tau})*S_0/tau^2);
    D_1t{tau} = (kron(C_1t{tau},C_1t{tau})*S_1/tau^2);
    D_1t_tran{tau} = D_1t{tau}(:,end);
end
D_0t_mat = cell2mat(D_0t);
D_1t_mat = cell2mat(D_1t_tran);

D_0y_temp_tau = D_0t_mat(1,:);
D_0r_temp_tau = D_0t_mat(4,:);
D_0yr_temp_tau = D_0t_mat(2,:);

D_1y_temp_tau = D_1t_mat(1,:);
D_1r_temp_tau = D_1t_mat(4,:);
D_1yr_temp_tau = D_1t_mat(3,:);
%{
% plots
figure(1)
subplot(3,2,1)
plot(C_0y_temp_tau)
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('$C_0$','Interpreter','LaTex','FontSize',15)
for i=1:5
    subplot(3,2,i+1)
    plot(C_1y_temp_mat_tau(i,:))
    xlim([0 520])
    xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
    if i==1
        title('$C_1$ loading for $\pi_t^1$','Interpreter','LaTex','FontSize',15)
    elseif i==2
        title('$C_1$ loading for $\pi_t^2$','Interpreter','LaTex','FontSize',15)
    elseif i==3
        title('$C_1$ loading for $m_t^1$','Interpreter','LaTex','FontSize',15)
    elseif i==4
        title('$C_1$ loading for $m_t^2$','Interpreter','LaTex','FontSize',15)
    elseif i==5
        title('$C_1$ loading for $\sigma_t^2$','Interpreter','LaTex','FontSize',15)
    end
end
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\textbf{For} \boldmath$y_t$','Interpreter','LaTex','FontSize',15,'HorizontalAlignment','center','VerticalAlignment', 'top')    
print(['Fig_y'],'-dpdf')
%--------------------------------------------------------------------------
figure(2)
subplot(3,2,1)
plot(C_0r_temp_tau)
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('$C_0$','Interpreter','LaTex','FontSize',15)
for i=1:5
    subplot(3,2,i+1)
    plot(C_1r_temp_mat_tau(i,:))
    xlim([0 520])
    xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
    if i==1
        title('$C_1$ loading for $\pi_t^1$','Interpreter','LaTex','FontSize',15)
    elseif i==2
        title('$C_1$ loading for $\pi_t^2$','Interpreter','LaTex','FontSize',15)
    elseif i==3
        title('$C_1$ loading for $m_t^1$','Interpreter','LaTex','FontSize',15)
    elseif i==4
        title('$C_1$ loading for $m_t^2$','Interpreter','LaTex','FontSize',15)
    elseif i==5
        title('$C_1$ loading for $\sigma_t^2$','Interpreter','LaTex','FontSize',15)
    end
end
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\textbf{For} \boldmath$r_t$','Interpreter','LaTex','FontSize',15,'HorizontalAlignment','center','VerticalAlignment', 'top')    
print(['Fig_r'],'-dpdf')
%--------------------------------------------------------------------------
figure(3)
subplot(3,2,1)
plot(D_0y_temp_tau)
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('For $Var_t(y_{t+1})$: $C_0$ ','Interpreter','LaTex','FontSize',13)
subplot(3,2,2)
plot(D_1y_temp_tau)
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('For $Var_t(y_{t+1})$: $C_1$ loading for $\sigma_t^2$','Interpreter','LaTex','FontSize',13)

subplot(3,2,3)
plot(D_0r_temp_tau)
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('For $Var_t(r_{t+1})$: $C_0$ ','Interpreter','LaTex','FontSize',13)
subplot(3,2,4)
plot(D_1r_temp_tau)
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('For $Var_t(r_{t+1})$: $C_1$ loading for $\sigma_t^2$','Interpreter','LaTex','FontSize',13)

subplot(3,2,5)
plot(D_0yr_temp_tau)
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('For $Var_t(r_{t+1})$: $C_0$ ','Interpreter','LaTex','FontSize',13)
subplot(3,2,6)
plot(D_1yr_temp_tau)
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('For $Cov_t(y_{t+1},r_{t+1})$: $C_1$ loading for $\sigma_t^2$','Interpreter','LaTex','FontSize',13)
%}
%--------------------------------------------------------------------------
%{
data_run.C_0y_temp_tau_orig = C_0y_temp_tau;
data_run.C_0r_temp_tau_orig = C_0r_temp_tau;
data_run.C_1y_temp_mat_tau_orig = C_1y_temp_mat_tau;
data_run.C_1r_temp_mat_tau_orig = C_1r_temp_mat_tau;
data_run.D_0y_temp_tau_orig = D_0y_temp_tau;
data_run.D_1y_temp_tau_orig = D_1y_temp_tau;
data_run.D_0r_temp_tau_orig = D_0r_temp_tau;
data_run.D_1r_temp_tau_orig = D_1r_temp_tau;
data_run.D_0yr_temp_tau_orig = D_0yr_temp_tau;
data_run.D_1yr_temp_tau_orig = D_1yr_temp_tau;
save('origdata.mat', '-struct', 'data_run');
%}

data_run.C_0y_temp_tau_origreest = C_0y_temp_tau;
data_run.C_0r_temp_tau_origreest = C_0r_temp_tau;
data_run.C_1y_temp_mat_tau_origreest = C_1y_temp_mat_tau;
data_run.C_1r_temp_mat_tau_origreest = C_1r_temp_mat_tau;
data_run.D_0y_temp_tau_origreest = D_0y_temp_tau;
data_run.D_1y_temp_tau_origreest = D_1y_temp_tau;
data_run.D_0r_temp_tau_origreest = D_0r_temp_tau;
data_run.D_1r_temp_tau_origreest = D_1r_temp_tau;
data_run.D_0yr_temp_tau_origreest = D_0yr_temp_tau;
data_run.D_1yr_temp_tau_origreest = D_1yr_temp_tau;
save('origdatareest.mat', '-struct', 'data_run');
