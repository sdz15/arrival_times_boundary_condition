txt = {strcat('theta=',string(theta)),strcat('mu=',string(mu)),strcat('sigma=',string(sigma)),strcat('k=',string(k)),strcat('omega=',string(omega)),strcat('L=',string(L))};

for x=1:N
    figure(1)
    plot(traj_free(x,:),tt_free,'Color',[0 0.4470 0.7410]);
    drawnow
    hold on;
    plot(traj_boundary(x,:),tt_boundary,'Color',[0.8500 0.3250 0.0980]);
    drawnow
    hold on;
end
hold off

figure(1)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories');
xlim([-L,L])
ylim([0 time])
legend({'free evolution','with boundary condition'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);

figure(2)
scatter(q0,yy_free);
hold on;
scatter(q0,yy_boundary);
hold on;
xlabel('Initial position of particle','FontSize',20);
ylabel('Time of detection','FontSize',20);
title('Detection time with respect to initial position');
xlim([-L,L])
ylim([0 time])
legend({'free evolution','with boundary condition'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(3)
cdfplot(yy_free);
hold on;
cdfplot(yy_boundary);
hold on;
title('Cumulative distribution of detection times');
xlim([0 time])
ylim([0 1])
legend({'free evolution','with boundary condition'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(4)
histogram(yy_free);
hold on;
histogram(yy_boundary);
hold on;
title('Histogram of detection times');
xlim([0 time])
legend({'free evolution','with boundary condition'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(5)
plot(times,mu_phi);
hold on;
plot(times,mu_psi);
hold on;
xlim([0 time])
title('Probability density function of detection time')
legend({'free evolution','with boundary condition'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(6)
plot(times,cumtrapz(times,mu_phi));
hold on;
plot(times,cumtrapz(times,mu_psi));
hold on;
xlim([0 time])
title('Cumulative distribution function of detection time')
legend({'free evolution','with boundary condition'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off