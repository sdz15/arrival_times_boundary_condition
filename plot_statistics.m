theta = pi/2;
mu = 0;
sigma = .1;
k = 0;
omega = 2;
time = 2;
N = 1;
mesh = .01;
step = 1/mesh;
L = 1;

times = (1e-1:mesh:time)+1e-6;
yy1 = Inf(1,N);
yy2 = Inf(1,N);

ode_num = length(times);
traj1 = Inf(N,ode_num);
traj2 = Inf(N,ode_num);

% Generating initial data
q0 = initialvals(sigma,mu,N,1);

for x=1:N
    x
    % [tt,qq] = ode45(@(t,q) velocity(abs(phiMinus(t,q,theta,sigma,k,mu,omega,step)).^2,abs(phiPlus(t,q,theta,sigma,k,mu,omega,step)).^2),times,q0(x)); % No boundary condition
    [tt1,qq1] = ode45(@(t,q) velocity(t,q,theta,sigma,k,mu,omega,L,step),times,q0(x)); % Boundary condition

    % f1 = min([find(qq>=L,1) find(qq<=-L,1)]); % Time of detection when no boundary condition
    % f2 = min([find(qq1>=L,1) find(qq1<=-L,1)]); %Time of detection with boundary condition

    % if ~isempty(f1)
    %     yy1(x) = tt(f1);
    % end
    % if ~isempty(f2)
    %     yy2(x) = tt1(f2);
    % end
end

% space = (-L:mesh:L);
% times = (0:mesh:time);
% 
% mu_phi = zeros(1,time/mesh+1);
% mu_psi = zeros(1,time/mesh+1);
% 
% for t = 1:time/mesh+1
%     mu_phi(t) = j0(abs(phiMinus(times(t),-L,theta,sigma,k,mu,omega,step)).^2,abs(phiPlus(times(t),-L,theta,sigma,k,mu,omega,step)).^2)+j0(abs(phiMinus(times(t),L,theta,sigma,k,mu,omega,step)).^2,abs(phiPlus(times(t),L,theta,sigma,k,mu,omega,step)).^2);
%     mu_psi(t) = j0(abs(psiMinus(times(t),-L,theta,sigma,k,mu,omega,L,step)).^2,abs(psiPlus(times(t),-L,theta,sigma,k,mu,omega,L,step)).^2)+j0(abs(psiMinus(times(t),L,theta,sigma,k,mu,omega,L,step)).^2,abs(psiPlus(times(t),L,theta,sigma,k,mu,omega,L,step)).^2);
% end
% 
% txt = {strcat('theta=',string(theta)),strcat('mu=',string(mu)),strcat('sigma=',string(sigma)),strcat('k=',string(k)),strcat('omega=',string(omega))};
% 
% for x=1:N
%     % figure(1)
%     % plot(qq,tt);
%     % xlim([-L,L])
%     % drawnow
%     % hold on;
% 
%     figure(2)
%     plot(qq1,tt1);
%     xlim([-L,L])
%     drawnow
%     hold on;
% end
% hold off
% 
% % figure(1)
% % xlabel('Position','FontSize',20);
% % ylabel('Time','FontSize',20);
% % title('Trajectories without boundary condition');
% % xlim([-L,L])
% % ylim([0 time])
% % xlimits=xlim;
% % ylimits=ylim;
% % text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);
% 
% figure(2)
% xlabel('Position','FontSize',20);
% ylabel('Time','FontSize',20);
% title('Trajectories with boundary condition');
% xlim([-L,L])
% ylim([0 time])
% xlimits=xlim;
% ylimits=ylim;
% text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);
% 
% % figure(3)
% % scatter(q0,yy1);
% % xlabel('Initial position of particle','FontSize',20);
% % ylabel('Time of detection','FontSize',20);
% % title('Detection time without boundary condition');
% % xlim([-L,L])
% % ylim([0 time])
% % xlimits=xlim;
% % ylimits=ylim;
% % text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
% 
% figure(4)
% scatter(q0,yy2);
% xlabel('Initial position of particle','FontSize',20);
% ylabel('Time of detection','FontSize',20);
% title('Detection time with boundary condition');
% xlim([-L,L])
% ylim([0 time])
% xlimits=xlim;
% ylimits=ylim;
% text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
% 
% % figure(5)
% % cdfplot(yy1);
% % title('Cumulative distribution of detection time without boundary condition');
% % xlim([0 time])
% % xlimits=xlim;
% % ylimits=ylim;
% % text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
% 
% figure(6)
% cdfplot(yy2);
% title('Cumulative distribution of detection time with boundary condition');
% xlim([0 time])
% xlimits=xlim;
% ylimits=ylim;
% text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
% 
% % figure(7)
% % histogram(yy1,N);
% % title('Histogram of detection time without boundary condition');
% % xlim([0 time])
% % xlimits=xlim;
% % ylimits=ylim;
% % text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
% 
% figure(8)
% histogram(yy2,N);
% title('Histogram of detection time with boundary condition');
% xlim([0 time])
% xlimits=xlim;
% ylimits=ylim;
% text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
% 
% % figure(9)
% % xlim([0 time])
% % plot(times,mu_phi);
% % title('Probability density function of detection time without boundary condition')
% % xlimits=xlim;
% % ylimits=ylim;
% % text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
% % 
% % figure(10)
% % xlim([0 time])
% % plot(times,mu_psi);
% % title('Probability density function of detection time with boundary condition')
% % xlimits=xlim;
% % ylimits=ylim;
% % text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
% % 
% % figure(11)
% % xlim([0 time])
% % plot(times,cumtrapz(times,mu_phi));
% % title('Cumulative distribution function of detection time without boundary condition')
% % xlimits=xlim;
% % ylimits=ylim;
% % text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
% % 
% % figure(12)
% % xlim([0 time])
% % plot(times,cumtrapz(times,mu_psi));
% % title('Cumulative distribution function of detection time with boundary condition')
% % xlimits=xlim;
% % ylimits=ylim;
% % text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);