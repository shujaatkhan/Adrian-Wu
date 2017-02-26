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

% plots
figure(1)
subplot(1,2,1)
plot(C_0y_temp_tau)
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
ylabel('\boldmath$\frac{1}{\tau}C_{0y}^\tau$','Interpreter','LaTex','FontSize',15)

subplot(1,2,2)
plot(C_0r_temp_tau)
xlim([0 520])
xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
ylabel('\boldmath$\frac{1}{\tau}C_{0r}^\tau$','Interpreter','LaTex','FontSize',15)
print(['FigC_0'],'-dpdf')

figure(2)
for i=1:5
    subplot(3,2,i)
    plot(C_1y_temp_mat_tau(i,:))
    xlim([0 520])
    xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
    if i==1
        ylabel('$$\pi_t^1$','Interpreter','LaTex','FontSize',15)
    elseif i==2
        ylabel('$\pi_t^2$','Interpreter','LaTex','FontSize',15)
    elseif i==3
        ylabel('$m_t^1$','Interpreter','LaTex','FontSize',15)
    elseif i==4
        ylabel('$m_t^2$','Interpreter','LaTex','FontSize',15)
    elseif i==5
        ylabel('$\sigma_t^2$','Interpreter','LaTex','FontSize',15)
    end
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\boldmath$\frac{1}{\tau}C_{1y}^\tau$ loadings for','Interpreter','LaTex','FontSize',15,'HorizontalAlignment','center','VerticalAlignment', 'top')    
end
print(['FigC_1y'],'-dpdf')

figure(3)
for i=1:5
    subplot(3,2,i)
    plot(C_1r_temp_mat_tau(i,:))
    xlim([0 520])
    xlabel('$\tau$','Interpreter','LaTex','FontSize',15)
    if i==1
        ylabel('$$\pi_t^1$','Interpreter','LaTex','FontSize',15)
    elseif i==2
        ylabel('$\pi_t^2$','Interpreter','LaTex','FontSize',15)
    elseif i==3
        ylabel('$m_t^1$','Interpreter','LaTex','FontSize',15)
    elseif i==4
        ylabel('$m_t^2$','Interpreter','LaTex','FontSize',15)
    elseif i==5
        ylabel('$\sigma_t^2$','Interpreter','LaTex','FontSize',15)
    end
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\boldmath$\frac{1}{\tau}C_{1r}^\tau$ loadings for','Interpreter','LaTex','FontSize',15,'HorizontalAlignment','center','VerticalAlignment', 'top')    
end
print(['FigC_1r'],'-dpdf')
