txt = {strcat('theta=',string(theta)),strcat('mu=',string(mu)),strcat('sigma=',string(sigma)),strcat('k=',string(k)),strcat('omega=',string(omega)),strcat('L=',string(L))};

figure(1)
hold on;
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories');
xlim([-L,L])
ylim([0 time])

for x=1:N
    figure(1)
    plot(traj_free(x,:),tt_free,'Color',[0 0.4470 0.7410]);
    drawnow
    plot(traj_boundary(x,:),tt_boundary,'Color',[0.8500 0.3250 0.0980]);
    drawnow
end

legend({'free evolution','with boundary condition'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);
hold off

figure(2)
hold on;
scatter(q0,yy_free);
scatter(q0,yy_boundary);
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
hold on;
cdfplot(yy_free);
cdfplot(yy_boundary);
plot(times,cumtrapz(times,mu_phi));
plot(times,cumtrapz(times,mu_psi));
xlabel('Time','Fontsize',20)
title('Cumulative distribution of detection times');
xlim([0 time])
ylim([0 1])
legend({'free evolution','with boundary condition','free evolution','with boundary condition'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(4)
hold on;
histogram(yy_free,N/10);
histogram(yy_boundary,N/10);
plot(times,mu_phi*10);
plot(times,mu_psi*10);
xlabel('Time','Fontsize',20)
title('Histogram and pdf of detection times');
xlim([0 time])
legend({'free evolution','with boundary condition','free evolution','with boundary condition'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off