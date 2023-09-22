theta = 0;
mu = 0;
sigma = .1;
k = 0;
omega = 2;
time = 2;
N = 5;
bound = 1;
mesh = .01;
L = 1;

times = (1e-1:mesh:time)+1e-6;

yy1 = Inf(1,N);
yy2 = Inf(1,N);

cmap = colormap;

% Generating initial data
q0 = initialvals(sigma,mu,N,1);

for x=1:N
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [tt,qq] = ode45(@(t,q) velocity(abs(phiMinus(t,q,theta,sigma,k,mu,omega)).^2,abs(phiPlus(t,q,theta,sigma,k,mu,omega)).^2),times,q0(x),opts); % No boundary condition
    [tt1,qq1] = ode45(@(t,q) velocity(abs(psiMinus(t,q,theta,sigma,k,mu,omega,L)).^2,abs(psiPlus(t,q,theta,sigma,k,mu,omega,L)).^2),times,q0(x),opts); % Boundary condition
    
    c = cmap(randi(size(cmap,1)),:);

    figure(1)
    plot(qq,tt,'Color',c);
    xlim([-L,L])
    drawnow
    hold on;

    figure(2)
    plot(qq1,tt1,'Color',c);
    xlim([-L,L])
    drawnow
    hold on;

    f1 = min([find(qq>=L,1) find(qq<=-L,1)]); % Time of detection when no boundary condition
    f2 = min([find(qq1>=L,1) find(qq1<=-L,1)]); %Time of detection with boundary condition

    if ~isempty(f1)
        yy1(x) = tt(f1);
    end
    if ~isempty(f2)
        yy2(x) = tt1(f2);
    end

end

txt = {strcat('theta=',string(theta)),strcat('mu=',string(mu)),strcat('sigma=',string(sigma)),strcat('k=',string(k)),strcat('omega=',string(omega))};

figure(1)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories without boundary condition');
xlim([-L,L])
ylim([0 time])
text(L-.5,.75,txt);

figure(2)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories with boundary condition');
xlim([-L,L])
ylim([0 time])
text(L-.5,.75,txt);

figure(3)
scatter(q0,yy1);
xlabel('Initial position of particle','FontSize',20);
ylabel('Time of detection','FontSize',20);
title('Detection time without boundary condition');
xlim([-L,L])
ylim([0 time])
text(L-.5,.75,txt);

figure(4)
scatter(q0,yy2);
xlabel('Initial position of particle','FontSize',20);
ylabel('Time of detection','FontSize',20);
title('Detection time with boundary condition');
xlim([-L,L])
ylim([0 time])
text(L-.5,.75,txt);

figure(5)
cdfplot(yy1);
title('Cumulative distribution of detection time without boundary condition');
xlim([0 time])
text(time-.5,.5,txt);

figure(6)
cdfplot(yy2);
title('Cumulative distribution of detection time with boundary condition');
xlim([0 time])
text(time-.5,.5,txt);

figure(7)
histogram(yy1,N);
title('Histogram of detection time without boundary condition');
xlim([0 time])
text(time-.5,.5,txt);

figure(8)
histogram(yy2,N);
title('Histogram of detection time with boundary condition');
xlim([0 time])
text(time-.5,.5,txt);

space = (-L:mesh:L);
times = (0:mesh:time);

fun_phi = zeros(time/mesh+1,2/mesh*L+1);
fun_psi = zeros(time/mesh+1,2/mesh*L+1);
mu_phi = zeros(1,time/mesh+1);
mu_psi = zeros(1,time/mesh+1);

for t = 1:time/mesh+1
    mu_phi(t) = j0(abs(phiMinus(times(t),-L,theta,sigma,k,mu,omega)).^2,abs(phiPlus(times(t),-L,theta,sigma,k,mu,omega)).^2)+j0(abs(phiMinus(times(t),L,theta,sigma,k,mu,omega)).^2,abs(phiPlus(times(t),L,theta,sigma,k,mu,omega)).^2);
    mu_psi(t) = j0(abs(psiMinus(times(t),-L,theta,sigma,k,mu,omega,L)).^2,abs(psiPlus(times(t),-L,theta,sigma,k,mu,omega,L)).^2)+j0(abs(psiMinus(times(t),L,theta,sigma,k,mu,omega,L)).^2,abs(psiPlus(times(t),L,theta,sigma,k,mu,omega,L)).^2);
end

figure(9)
xlim([0 time])
plot(times,mu_phi);
title('Probability density function of detection time without boundary condition')
text(time-.5,.5,txt);

figure(10)
xlim([0 time])
plot(times,mu_psi);
title('Probability density function of detection time with boundary condition')
text(time-.5,.5,txt);

figure(11)
xlim([0 time])
plot(times,cumtrapz(times,mu_phi));
title('Cumulative distribution function of detection time without boundary condition')
text(time-.5,.5,txt);

figure(12)
xlim([0 time])
plot(times,cumtrapz(times,mu_psi));
title('Cumulative distribution function of detection time with boundary condition')
text(time-.5,.5,txt);

save('single_particle_variables.mat')