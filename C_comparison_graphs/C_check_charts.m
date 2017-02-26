clear all

% load data
load origdata
load origdatareest
load largedata
load s1data


% plots
figure(1)
subplot(3,2,1)
hold on
plot(C_0y_temp_tau_orig,'LineWidth',2)
plot(C_0y_temp_tau_origreest,':','LineWidth',2)
plot(C_0y_temp_tau_large,'--','LineWidth',2)
plot(C_0y_temp_tau_s1,'-.','LineWidth',2)
hold off
legend('Orig','Orig re-est','Full','Full Set1')
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('$C_0$','Interpreter','LaTex','FontSize',15)
for i=1:5
    subplot(3,2,i+1)
    hold on
    plot(C_1y_temp_mat_tau_orig(i,:),'LineWidth',2)
    plot(C_1y_temp_mat_tau_origreest(i,:),':','LineWidth',2)
    plot(C_1y_temp_mat_tau_large(i,:),'--','LineWidth',2)
    plot(C_1y_temp_mat_tau_s1(i,:),'-.','LineWidth',2)
    hold off
    legend('Orig','Orig re-est','Full','Full Set1')
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
print('-fillpage','Fig_y','-dpdf')
%--------------------------------------------------------------------------
figure(2)
subplot(3,2,1)
hold on
plot(C_0r_temp_tau_orig,'LineWidth',2)
plot(C_0r_temp_tau_origreest,':','LineWidth',2)
plot(C_0r_temp_tau_large,'--','LineWidth',2)
plot(C_0r_temp_tau_s1,'-.','LineWidth',2)
hold off
legend('Orig','Orig re-est','Full','Full Set1')
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('$C_0$','Interpreter','LaTex','FontSize',15)
for i=1:5
    subplot(3,2,i+1)
    hold on
    plot(C_1r_temp_mat_tau_orig(i,:),'LineWidth',2)
    plot(C_1r_temp_mat_tau_origreest(i,:),':','LineWidth',2)
    plot(C_1r_temp_mat_tau_large(i,:),'--','LineWidth',2)
    plot(C_1r_temp_mat_tau_s1(i,:),'-.','LineWidth',2)
    hold off
    legend('Orig','Orig re-est','Full','Full Set1')
    xlim([0 520])
    xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
    legend('Orig','Orig re-est','Full','Full Set1')
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
print('-fillpage','Fig_r','-dpdf')
%--------------------------------------------------------------------------
figure(3)
subplot(3,2,1)
hold on
plot(D_0y_temp_tau_orig,'LineWidth',2)
plot(D_0y_temp_tau_origreest,':','LineWidth',2)
plot(D_0y_temp_tau_large,'--','LineWidth',2)
plot(D_0y_temp_tau_s1,'-.','LineWidth',2)
hold off
legend('Orig','Orig re-est','Full','Full Set1')
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('For $Var_t(y_{t+1})$: $C_0$ ','Interpreter','LaTex','FontSize',13)
subplot(3,2,2)
hold on
plot(D_1y_temp_tau_orig,'LineWidth',2)
plot(D_1y_temp_tau_origreest,':','LineWidth',2)
plot(D_1y_temp_tau_large,'--','LineWidth',2)
plot(D_1y_temp_tau_s1,'-.','LineWidth',2)
hold off
legend('Orig','Orig re-est','Full','Full Set1')
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('For $Var_t(y_{t+1})$: $C_1$ loading for $\sigma_t^2$','Interpreter','LaTex','FontSize',13)

subplot(3,2,3)
hold on
plot(D_0r_temp_tau_orig,'LineWidth',2)
plot(D_0r_temp_tau_origreest,':','LineWidth',2)
plot(D_0r_temp_tau_large,'--','LineWidth',2)
plot(D_0r_temp_tau_s1,'-.','LineWidth',2)
hold off
legend('Orig','Orig re-est','Full','Full Set1')
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('For $Var_t(r_{t+1})$: $C_0$ ','Interpreter','LaTex','FontSize',13)
subplot(3,2,4)
hold on
plot(D_1r_temp_tau_orig,'LineWidth',2)
plot(D_1r_temp_tau_origreest,':','LineWidth',2)
plot(D_1r_temp_tau_large,'--','LineWidth',2)
plot(D_1r_temp_tau_s1,'-.','LineWidth',2)
hold off
legend('Orig','Orig re-est','Full','Full Set1')
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('For $Var_t(r_{t+1})$: $C_1$ loading for $\sigma_t^2$','Interpreter','LaTex','FontSize',13)

subplot(3,2,5)
hold on
plot(D_0yr_temp_tau_orig,'LineWidth',2)
plot(D_0yr_temp_tau_origreest,':','LineWidth',2)
plot(D_0yr_temp_tau_large,'--','LineWidth',2)
plot(D_0yr_temp_tau_s1,'-.','LineWidth',2)
hold off
legend('Orig','Orig re-est','Full','Full Set1')
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('For $Var_t(r_{t+1})$: $C_0$ ','Interpreter','LaTex','FontSize',13)
subplot(3,2,6)
hold on
plot(D_1yr_temp_tau_orig,'LineWidth',2)
plot(D_1yr_temp_tau_origreest,':','LineWidth',2)
plot(D_1yr_temp_tau_large,'--','LineWidth',2)
plot(D_1yr_temp_tau_s1,'-.','LineWidth',2)
hold off
legend('Orig','Orig re-est','Full','Full Set1')
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
title('For $Cov_t(y_{t+1},r_{t+1})$: $C_1$ loading for $\sigma_t^2$','Interpreter','LaTex','FontSize',13)

print('-fillpage','Fig_VCovs','-dpdf')
