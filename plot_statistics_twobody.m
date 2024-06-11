txt = {strcat('theta\_1=',string(theta_1)),strcat('theta\_2=',string(theta_2)),strcat('mu\_1=',string(mu_1)),strcat('mu\_2=',string(mu_2)),strcat('alpha\_1=',string(sigma_1)),strcat('alpha\_2=',string(sigma_2)),strcat('k\_1=',string(k_1)),strcat('k\_2=',string(k_2)),strcat('omega=',string(omega)),strcat('L=',string(L))};

figure(1)
hold on;
xlim([-L,L])
ylim([0 time])

figure(2)
hold on;
xlim([-L,L])
ylim([0 time])

figure(3)
hold on;
xlim([-L,L])
ylim([0 time])

figure(4)
hold on;
xlim([-L,L])
ylim([0 time])

for x=1:N
    figure(1)
    lgd_1(1:2) = plot(traj_phi_pure_1(:,x),tt_phi_pure,traj_phi_pure_2(:,x),tt_phi_pure,'Color',[0 0.4470 0.7410]);
    drawnow
    lgd_1(3:4) = plot(traj_phi_pure_shifted_1(:,x),tt_phi_pure_shifted,traj_phi_pure_shifted_2(:,x),tt_phi_pure_shifted,'Color',[0.8500 0.3250 0.0980]);
    drawnow

    figure(2)
    lgd_2(1:2) = plot(traj_phi_entangled_1(:,x),tt_phi_entangled,traj_phi_entangled_2(:,x),tt_phi_entangled,'Color',[0 0.4470 0.7410]);
    drawnow
    lgd_2(3:4) = plot(traj_phi_entangled_shifted_1(:,x),tt_phi_entangled_shifted,traj_phi_entangled_shifted_2(:,x),tt_phi_entangled_shifted,'Color',[0.8500 0.3250 0.0980]);
    drawnow

    figure(3)
    lgd_3(1:2) = plot(traj_psi_pure_1(:,x),tt_psi_pure,traj_psi_pure_2(:,x),tt_psi_pure,'Color',[0 0.4470 0.7410]);
    drawnow
    lgd_3(3:4) = plot(traj_psi_pure_shifted_1(:,x),tt_psi_pure_shifted,traj_psi_pure_shifted_2(:,x),tt_psi_pure_shifted,'Color',[0.8500 0.3250 0.0980]);
    drawnow

    figure(4)
    lgd_4(1:2) = plot(traj_psi_entangled_1(:,x),tt_psi_entangled,traj_psi_entangled_2(:,x),tt_psi_entangled,'Color',[0 0.4470 0.7410]);
    drawnow
    lgd_4(3:4) = plot(traj_psi_entangled_shifted_1(:,x),tt_psi_entangled_shifted,traj_psi_entangled_shifted_2(:,x),tt_psi_entangled_shifted,'Color',[0.8500 0.3250 0.0980]);
    drawnow
end

figure(1)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for pure product without boundary condition');
legend(lgd_1([2,4]),{'normal','shifted'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);
hold off

figure(2)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for entangled pair without boundary condition');
legend(lgd_2([2,4]),{'normal','shifted'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);
hold off

figure(3)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for pure product with boundary condition');
legend(lgd_3([2,4]),{'normal','shifted'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);
hold off

figure(4)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for entangled pair with boundary condition');
legend(lgd_4([2,4]),{'normal','shifted'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);
hold off

figure(5)
hold on;
scatter(yy_phi_pure_A,yy_phi_pure_B);
scatter(yy_phi_pure_A_shifted,yy_phi_pure_B_shifted);
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
hold on;
scatter(yy_phi_entangled_A,yy_phi_entangled_B);
scatter(yy_phi_entangled_A_shifted,yy_phi_entangled_B_shifted);
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
hold on;
scatter(yy_psi_pure_A,yy_psi_pure_B);
scatter(yy_psi_pure_A_shifted,yy_psi_pure_B_shifted);
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
hold on;
scatter(yy_psi_entangled_A,yy_psi_entangled_B);
scatter(yy_psi_entangled_A_shifted,yy_psi_entangled_B_shifted);
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
hold on;
cdfplot(yy_phi_pure_A);
cdfplot(yy_phi_pure_B);
cdfplot(yy_phi_pure_A_shifted);
cdfplot(yy_phi_pure_B_shifted);
plot(times,cumtrapz(times,mu_phi_pure_A));
plot(times,cumtrapz(times,mu_phi_pure_B));
plot(times,cumtrapz(times,mu_phi_pure_A_shifted));
plot(times,cumtrapz(times,mu_phi_pure_B_shifted));
xlabel('Time','Fontsize',20)
title('Cumulative distributions of arrival time for pure product without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Alice','Bob','Alice shifted','Bob shifted','Location','southeast');
xlim([0 time])
ylim([0 1])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(10)
hold on;
cdfplot(yy_phi_entangled_A);
cdfplot(yy_phi_entangled_B);
cdfplot(yy_phi_entangled_A_shifted);
cdfplot(yy_phi_entangled_B_shifted);
plot(times,cumtrapz(times,mu_phi_entangled_A));
plot(times,cumtrapz(times,mu_phi_entangled_B));
plot(times,cumtrapz(times,mu_phi_entangled_A_shifted));
plot(times,cumtrapz(times,mu_phi_entangled_B_shifted));
xlabel('Time','Fontsize',20)
title('Cumulative distributions of arrival time for entangled pair without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Alice','Bob','Alice shifted','Bob shifted','Location','southeast');
xlim([0 time])
ylim([0 1])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(11)
hold on;
cdfplot(yy_psi_pure_A);
cdfplot(yy_psi_pure_B);
cdfplot(yy_psi_pure_A_shifted);
cdfplot(yy_psi_pure_B_shifted);
plot(times,cumtrapz(times,mu_psi_pure_A));
plot(times,cumtrapz(times,mu_psi_pure_B));
plot(times,cumtrapz(times,mu_psi_pure_A_shifted));
plot(times,cumtrapz(times,mu_psi_pure_B_shifted));
xlabel('Time','Fontsize',20)
title('Cumulative distributions of arrival time for pure product with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Alice','Bob','Alice shifted','Bob shifted','Location','southeast');
xlim([0 time])
ylim([0 1])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(12)
hold on;
cdfplot(yy_psi_entangled_A);
cdfplot(yy_psi_entangled_B);
cdfplot(yy_psi_entangled_A_shifted);
cdfplot(yy_psi_entangled_B_shifted);
plot(times,cumtrapz(times,mu_psi_entangled_A));
plot(times,cumtrapz(times,mu_psi_entangled_B));
plot(times,cumtrapz(times,mu_psi_entangled_A_shifted));
plot(times,cumtrapz(times,mu_psi_entangled_B_shifted));
xlabel('Time','Fontsize',20)
title('Cumulative distributions of arrival time for entangled pair with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Alice','Bob','Alice shifted','Bob shifted','Location','southeast');
xlim([0 time])
ylim([0 1])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(13)
hold on;
histogram(yy_phi_pure_A);
histogram(yy_phi_pure_B);
histogram(yy_phi_pure_A_shifted);
histogram(yy_phi_pure_B_shifted);
plot(times,mu_phi_pure_A);
plot(times,mu_phi_pure_B);
plot(times,mu_phi_pure_A_shifted);
plot(times,mu_phi_pure_B_shifted);
xlabel('Time','Fontsize',20)
title('Histogram and pdf of arrival time for pure product without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,ylimits(2)*5/8,txt);
hold off

figure(14)
hold on;
histogram(yy_phi_entangled_A);
histogram(yy_phi_entangled_B);
histogram(yy_phi_entangled_A_shifted);
histogram(yy_phi_entangled_B_shifted);
plot(times,mu_phi_entangled_A);
plot(times,mu_phi_entangled_B);
plot(times,mu_phi_entangled_A_shifted);
plot(times,mu_phi_entangled_B_shifted);
xlabel('Time','Fontsize',20)
title('Histogram and pdf of arrival time for entangled pair without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,ylimits(2)*5/8,txt);
hold off

figure(15)
hold on;
histogram(yy_psi_pure_A);
histogram(yy_psi_pure_B);
histogram(yy_psi_pure_A_shifted);
histogram(yy_psi_pure_B_shifted);
plot(times,mu_psi_pure_A);
plot(times,mu_psi_pure_B);
plot(times,mu_psi_pure_A_shifted);
plot(times,mu_psi_pure_B_shifted);
xlabel('Time','Fontsize',20)
title('Histogram and pdf of arrival time for pure product with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,ylimits(2)*5/8,txt);
hold off

figure(16)
hold on;
histogram(yy_psi_entangled_A);
histogram(yy_psi_entangled_B);
histogram(yy_psi_entangled_A_shifted);
histogram(yy_psi_entangled_B_shifted);
plot(times,mu_psi_entangled_A);
plot(times,mu_psi_entangled_B);
plot(times,mu_psi_entangled_A_shifted);
plot(times,mu_psi_entangled_B_shifted);
xlabel('Time','Fontsize',20)
title('Histogram and pdf of arrival time for entangled pair with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,ylimits(2)*5/8,txt);
hold off