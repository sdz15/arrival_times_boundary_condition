theta = pi/2;
mu = 0;
sigma = .1;
k = 0;
omega = 2;
time = 2;
N = 2;
mesh = .1;
step = 1/mesh;
L = 1;

times = (1e-1:mesh:time)+1e-6;
yy_free = Inf(1,N);
yy_boundary = Inf(1,N);

ode_num = length(times);
traj_free = Inf(N,ode_num);
traj_boundary = Inf(N,ode_num);

% Generating initial data
q0 = initialvals(sigma,mu,N,1);

for x=1:N
    x
    [tt_free,qq_free] = ode45(@(t,q) velocity_single_electron(abs(phiMinus(t,q,theta,sigma,k,mu,omega,step)).^2,abs(phiPlus(t,q,theta,sigma,k,mu,omega,step)).^2),times,q0(x)); % Free evolution
    [tt_boundary,qq_boundary] = ode45(@(t,q) velocity_single_electron(abs(psiMinus(t,q,theta,sigma,k,mu,omega,L,step)).^2,abs(psiPlus(t,q,theta,sigma,k,mu,omega,L,step)).^2),times,q0(x)); % Boundary condition

    traj_free(x,:) = qq_free;
    traj_boundary(x,:) = qq_boundary;

    f_free = min([find(traj_free(x,:)>=L,1) find(traj_free(x,:)<=-L,1)]); % Time of detection with free evolution
    f_boundary = min([find(traj_boundary(x,:)>=L,1) find(traj_boundary(x,:)<=-L,1)]); % Time of detection with boundary condition

    if ~isempty(f_free)
        yy_free(x) = tt_free(f_free);
    end
    if ~isempty(f_boundary)
        yy_boundary(x) = tt_boundary(f_boundary);
    end
end

space = (-L:mesh:L);
times = (0:mesh:time);

mu_phi = zeros(1,time/mesh+1);
mu_psi = zeros(1,time/mesh+1);

for t = 1:time/mesh+1
    mu_phi(t) = j0(abs(phiMinus(times(t),-L,theta,sigma,k,mu,omega,step)).^2,abs(phiPlus(times(t),-L,theta,sigma,k,mu,omega,step)).^2)+j0(abs(phiMinus(times(t),L,theta,sigma,k,mu,omega,step)).^2,abs(phiPlus(times(t),L,theta,sigma,k,mu,omega,step)).^2);
    mu_psi(t) = j0(abs(psiMinus(times(t),-L,theta,sigma,k,mu,omega,L,step)).^2,abs(psiPlus(times(t),-L,theta,sigma,k,mu,omega,L,step)).^2)+j0(abs(psiMinus(times(t),L,theta,sigma,k,mu,omega,L,step)).^2,abs(psiPlus(times(t),L,theta,sigma,k,mu,omega,L,step)).^2);
end

% TO PLOT THE STATISTICS, GO TO PLOT_STATISTICS_SINGLE.M