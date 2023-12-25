txt = {strcat('theta\_1=',string(theta_1)),strcat('theta\_2=',string(theta_2)),strcat('mu\_1=',string(mu_1)),strcat('mu\_2=',string(mu_2)),strcat('sigma\_1=',string(sigma_1)),strcat('sigma\_2=',string(sigma_2)),strcat('k\_1=',string(k_1)),strcat('k\_2=',string(k_2)),strcat('omega=',string(omega)),strcat('L=',string(L))};

for x=1:N
    figure(1)
    lgd_1(1:2) = plot(traj_phi_pure_1(:,x),tt_phi_pure,traj_phi_pure_2(:,x),tt_phi_pure,'Color',[0 0.4470 0.7410]);
    drawnow
    hold on;
    lgd_1(3:4) = plot(traj_phi_pure_shifted_1(:,x),tt_phi_pure_shifted,traj_phi_pure_shifted_2(:,x),tt_phi_pure_shifted,'Color',[0.8500 0.3250 0.0980]);
    drawnow
    hold on;

    figure(2)
    lgd_2(1:2) = plot(traj_phi_entangled_1(:,x),tt_phi_entangled,traj_phi_entangled_2(:,x),tt_phi_entangled,'Color',[0 0.4470 0.7410]);
    drawnow
    hold on;
    lgd_2(3:4) = plot(traj_phi_entangled_shifted_1(:,x),tt_phi_entangled_shifted,traj_phi_entangled_shifted_2(:,x),tt_phi_entangled_shifted,'Color',[0.8500 0.3250 0.0980]);
    drawnow
    hold on;

    figure(3)
    lgd_3(1:2) = plot(traj_psi_pure_1(:,x),tt_psi_pure,traj_psi_pure_2(:,x),tt_psi_pure,'Color',[0 0.4470 0.7410]);
    drawnow
    hold on;
    lgd_3(3:4) = plot(traj_psi_pure_shifted_1(:,x),tt_psi_pure_shifted,traj_psi_pure_shifted_2(:,x),tt_psi_pure_shifted,'Color',[0.8500 0.3250 0.0980]);
    drawnow
    hold on;

    figure(4)
    lgd_4(1:2) = plot(traj_psi_entangled_1(:,x),tt_psi_entangled,traj_psi_entangled_2(:,x),tt_psi_entangled,'Color',[0 0.4470 0.7410]);
    drawnow
    hold on;
    lgd_4(3:4) = plot(traj_psi_entangled_shifted_1(:,x),tt_psi_entangled_shifted,traj_psi_entangled_shifted_2(:,x),tt_psi_entangled_shifted,'Color',[0.8500 0.3250 0.0980]);
    drawnow
    hold on;
end

figure(1)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for pure product without boundary condition');
legend(lgd_1([2,4]),{'normal','shifted'},'Location','southwest');
xlim([-L,L])
ylim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);
hold off;

figure(2)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for entangled pair without boundary condition');
legend(lgd_2([2,4]),{'normal','shifted'},'Location','southwest');
xlim([-L,L])
ylim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);
hold off;

figure(3)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for pure product with boundary condition');
legend(lgd_3([2,4]),{'normal','shifted'},'Location','southwest');
xlim([-L,L])
ylim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);
hold off;

figure(4)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for entangled pair with boundary condition');
legend(lgd_4([2,4]),{'normal','shifted'},'Location','southwest');
xlim([-L,L])
ylim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);
hold off;

figure(5)
scatter(yy_phi_pure_A,yy_phi_pure_B);
hold on;
scatter(yy_phi_pure_A_shifted,yy_phi_pure_B_shifted);
hold on;
xlabel('Arrival time at Alice','FontSize',20);
ylabel('Arrival time at Bob','FontSize',20);
title('Arrival times for pure product without boundary condition');
legend({'normal','shifted'},'Location','southwest');
xlim([0 time])
ylim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(6)
scatter(yy_phi_entangled_A,yy_phi_entangled_B);
hold on;
scatter(yy_phi_entangled_A_shifted,yy_phi_entangled_B_shifted);
hold on;
xlabel('Arrival time at Alice','FontSize',20);
ylabel('Arrival time at Bob','FontSize',20);
title('Arrival times for entangled pair without boundary condition');
legend({'normal','shifted'},'Location','southwest');
xlim([0 time])
ylim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(7)
scatter(yy_psi_pure_A,yy_psi_pure_B);
hold on;
scatter(yy_psi_pure_A_shifted,yy_psi_pure_B_shifted);
hold on;
xlabel('Arrival time at Alice','FontSize',20);
ylabel('Arrival time at Bob','FontSize',20);
title('Arrival times for pure product with boundary condition');
legend({'normal','shifted'},'Location','southwest');
xlim([0 time])
ylim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(8)
scatter(yy_psi_entangled_A,yy_psi_entangled_B);
hold on;
scatter(yy_psi_entangled_A_shifted,yy_psi_entangled_B_shifted);
hold on;
xlabel('Arrival time at Alice','FontSize',20);
ylabel('Arrival time at Bob','FontSize',20);
title('Arrival times for entangled pair with boundary condition');
legend({'normal','shifted'},'Location','southwest');
xlim([0 time])
ylim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(9)
cdfplot(yy_phi_pure_A);
hold on
cdfplot(yy_phi_pure_B);
hold on
cdfplot(yy_phi_pure_A_shifted);
hold on
cdfplot(yy_phi_pure_B_shifted);
hold on
xlabel('Time','Fontsize',20)
title('Cumulative distribution of arrival time for pure product without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast');
xlim([0 time])
ylim([0 1])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(10)
cdfplot(yy_phi_entangled_A);
hold on
cdfplot(yy_phi_entangled_B);
hold on
cdfplot(yy_phi_entangled_A_shifted);
hold on
cdfplot(yy_phi_entangled_B_shifted);
hold on
xlabel('Time','Fontsize',20)
title('Cumulative distribution of arrival time for entangled pair without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast');
xlim([0 time])
ylim([0 1])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(11)
cdfplot(yy_psi_pure_A);
hold on
cdfplot(yy_psi_pure_B);
hold on
cdfplot(yy_psi_pure_A_shifted);
hold on
cdfplot(yy_psi_pure_B_shifted);
hold on
xlabel('Time','Fontsize',20)
title('Cumulative distribution of arrival time for pure product with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast');
xlim([0 time])
ylim([0 1])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(12)
cdfplot(yy_psi_entangled_A);
hold on
cdfplot(yy_psi_entangled_B);
hold on
cdfplot(yy_psi_entangled_A_shifted);
hold on
cdfplot(yy_psi_entangled_B_shifted);
hold on
xlabel('Time','Fontsize',20)
title('Cumulative distribution of arrival time for entangled pair with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast');
xlim([0 time])
ylim([0 1])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(13)
histogram(yy_phi_pure_A);
hold on
histogram(yy_phi_pure_B);
hold on
histogram(yy_phi_pure_A_shifted);
hold on
histogram(yy_phi_pure_B_shifted);
hold on
xlabel('Time','Fontsize',20)
title('Histogram of arrival time for pure product without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,ylimits(2)*5/8,txt);
hold off

figure(14)
histogram(yy_phi_entangled_A);
hold on
histogram(yy_phi_entangled_B);
hold on
histogram(yy_phi_entangled_A_shifted);
hold on
histogram(yy_phi_entangled_B_shifted);
hold on
xlabel('Time','Fontsize',20)
title('Histogram of arrival time for entangled pair without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,ylimits(2)*5/8,txt);
hold off

figure(15)
histogram(yy_psi_pure_A);
hold on
histogram(yy_psi_pure_B);
hold on
histogram(yy_psi_pure_A_shifted);
hold on
histogram(yy_psi_pure_B_shifted);
hold on
xlabel('Time','Fontsize',20)
title('Histogram of arrival time for pure product with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,ylimits(2)*5/8,txt);
hold off

figure(16)
histogram(yy_psi_entangled_A);
hold on
histogram(yy_psi_entangled_B);
hold on
histogram(yy_psi_entangled_A_shifted);
hold on
histogram(yy_psi_entangled_B_shifted);
hold on
xlabel('Time','Fontsize',20)
title('Histogram of arrival time for entangled pair with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,ylimits(2)*5/8,txt);
hold off

figure(17)
plot(times,mu_phi_pure_A);
hold on
plot(times,mu_phi_pure_B);
hold on
plot(times,mu_phi_pure_A_shifted);
hold on
plot(times,mu_phi_pure_B_shifted);
hold on
xlabel('Time','Fontsize',20)
title('Probability density function of arrival times for pure product without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(18)
plot(times,mu_phi_entangled_A);
hold on
plot(times,mu_phi_entangled_B);
hold on
plot(times,mu_phi_entangled_A_shifted);
hold on
plot(times,mu_phi_entangled_B_shifted);
hold on
xlabel('Time','Fontsize',20)
title('Probability density function of arrival times for entangled pair without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(19)
plot(times,mu_phi_pure_A);
hold on
plot(times,mu_phi_pure_B);
hold on
plot(times,mu_phi_pure_A_shifted);
hold on
plot(times,mu_phi_pure_B_shifted);
hold on
xlabel('Time','Fontsize',20)
title('Probability density function of arrival times for pure product with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(20)
plot(times,mu_phi_entangled_A);
hold on
plot(times,mu_phi_entangled_B);
hold on
plot(times,mu_phi_entangled_A_shifted);
hold on
plot(times,mu_phi_entangled_B_shifted);
hold on
xlabel('Time','Fontsize',20)
title('Probability density function of arrival times for entangled pair with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off