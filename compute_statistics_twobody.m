theta_1 = 0;
theta_2 = pi/2;
mu_1 = 0;
mu_2 = 0;
sigma_1 = .005;
sigma_2 = .005;
k_1 = 0;
k_2 = 0;
omega = 2;
L = 1;
time = 2;
N = 200;
mesh = .01;
step = min(10/mesh,10);

delta = .1;
L1 = L+delta/2;
mu_3 = mu_1-delta/2;
mu_4 = mu_2-delta/2;

times = (0:mesh:time)+1e-6;

initvals = initialvals([sigma_1,sigma_2],[mu_1,mu_2],N,2);
initvals_shifted = initialvals([sigma_1,sigma_2],[mu_3,mu_4],N,2);

yy_phi_pure_A = Inf(1,N);
yy_phi_pure_B = Inf(1,N); 
yy_phi_entangled_A = Inf(1,N);
yy_phi_entangled_B = Inf(1,N); 
yy_psi_pure_A = Inf(1,N);
yy_psi_pure_B = Inf(1,N); 
yy_psi_entangled_A = Inf(1,N);
yy_psi_entangled_B = Inf(1,N); 

yy_phi_pure_A_shifted = Inf(1,N);
yy_phi_pure_B_shifted = Inf(1,N); 
yy_phi_entangled_A_shifted = Inf(1,N);
yy_phi_entangled_B_shifted = Inf(1,N); 
yy_psi_pure_A_shifted = Inf(1,N);
yy_psi_pure_B_shifted = Inf(1,N); 
yy_psi_entangled_A_shifted = Inf(1,N);
yy_psi_entangled_B_shifted = Inf(1,N);

ode_num = length(times);
traj_phi_pure_1 = Inf(ode_num,N);
traj_phi_pure_2 = Inf(ode_num,N);
traj_phi_entangled_1 = Inf(ode_num,N);
traj_phi_entangled_2 = Inf(ode_num,N);
traj_psi_pure_1 = Inf(ode_num,N);
traj_psi_pure_2 = Inf(ode_num,N);
traj_psi_entangled_1 = Inf(ode_num,N);
traj_psi_entangled_2 = Inf(ode_num,N); 
traj_phi_pure_shifted_1 = Inf(ode_num,N);
traj_phi_pure_shifted_2 = Inf(ode_num,N);
traj_phi_entangled_shifted_1 = Inf(ode_num,N);
traj_phi_entangled_shifted_2 = Inf(ode_num,N);
traj_psi_pure_shifted_1 = Inf(ode_num,N);
traj_psi_pure_shifted_2 = Inf(ode_num,N);
traj_psi_entangled_shifted_1 = Inf(ode_num,N);
traj_psi_entangled_shifted_2 = Inf(ode_num,N); 

for x=1:N
    x
    q0 = initvals(:,x);
    q0_shifted = initvals_shifted(:,x);

    [tt_phi_pure,qq_phi_pure] = ode45(@(t,q) velocity_twobody(phiMinus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega,step),phiMinus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega,step),phiMinus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega,step),phiMinus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega,step),phiPlus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega,step),phiPlus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega,step),phiPlus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega,step),phiPlus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega,step),false),times,q0);
    [tt_phi_entangled,qq_phi_entangled] = ode45(@(t,q) velocity_twobody(phiMinus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega,step),phiMinus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega,step),phiMinus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega,step),phiMinus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega,step),phiPlus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega,step),phiPlus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega,step),phiPlus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega,step),phiPlus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega,step),true),times,q0);

    [tt_psi_pure,qq_psi_pure] = ode45(@(t,q) velocity_twobody(psiMinus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega,L,step),psiMinus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega,L,step),psiMinus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega,L,step),psiMinus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega,L,step),psiPlus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega,L,step),psiPlus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega,L,step),psiPlus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega,L,step),psiPlus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega,L,step),false),times,q0);
    [tt_psi_entangled,qq_psi_entangled] = ode45(@(t,q) velocity_twobody(psiMinus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega,L,step),psiMinus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega,L,step),psiMinus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega,L,step),psiMinus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega,L,step),psiPlus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega,L,step),psiPlus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega,L,step),psiPlus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega,L,step),psiPlus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega,L,step),true),times,q0);

    % SAME STATISTICS ABOVE FOR WHEN L IS SHIFTED BY DELTA/2
    [tt_phi_pure_shifted,qq_phi_pure_shifted] = ode45(@(t,q) velocity_twobody(phiMinus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega,step),phiMinus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega,step),phiMinus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega,step),phiMinus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega,step),phiPlus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega,step),phiPlus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega,step),phiPlus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega,step),phiPlus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega,step),false),times,q0_shifted);
    [tt_phi_entangled_shifted,qq_phi_entangled_shifted] = ode45(@(t,q) velocity_twobody(phiMinus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega,step),phiMinus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega,step),phiMinus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega,step),phiMinus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega,step),phiPlus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega,step),phiPlus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega,step),phiPlus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega,step),phiPlus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega,step),true),times,q0_shifted);

    [tt_psi_pure_shifted,qq_psi_pure_shifted] = ode45(@(t,q) velocity_twobody(psiMinus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega,L1,step),psiMinus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega,L1,step),psiMinus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega,L1,step),psiMinus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega,L1,step),psiPlus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega,L1,step),psiPlus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega,L1,step),psiPlus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega,L1,step),psiPlus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega,L1,step),false),times,q0_shifted);
    [tt_psi_entangled_shifted,qq_psi_entangled_shifted] = ode45(@(t,q) velocity_twobody(psiMinus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega,L1,step),psiMinus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega,L1,step),psiMinus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega,L1,step),psiMinus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega,L1,step),psiPlus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega,L,step),psiPlus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega,L1,step),psiPlus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega,L1,step),psiPlus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega,L1,step),true),times,q0_shifted);

    traj_phi_pure_1(:,x) = qq_phi_pure(:,1);
    traj_phi_pure_2(:,x) = qq_phi_pure(:,2);
    traj_phi_entangled_1(:,x) = qq_phi_entangled(:,1);
    traj_phi_entangled_2(:,x) = qq_phi_entangled(:,2);
    traj_psi_pure_1(:,x) = qq_psi_pure(:,1);
    traj_psi_pure_2(:,x) = qq_psi_pure(:,2);
    traj_psi_entangled_1(:,x) = qq_psi_entangled(:,1);
    traj_psi_entangled_2(:,x) = qq_psi_entangled(:,2); 
    traj_phi_pure_shifted_1(:,x) = qq_phi_pure(:,1);
    traj_phi_pure_shifted_2(:,x) = qq_phi_pure(:,2);
    traj_phi_entangled_shifted_1(:,x) = qq_phi_entangled_shifted(:,1);
    traj_phi_entangled_shifted_2(:,x) = qq_phi_entangled_shifted(:,2);
    traj_psi_pure_shifted_1(:,x) = qq_psi_pure(:,1);
    traj_psi_pure_shifted_2(:,x) = qq_psi_pure(:,2);
    traj_psi_entangled_shifted_1(:,x) = qq_psi_entangled_shifted(:,1);
    traj_psi_entangled_shifted_2(:,x) = qq_psi_entangled_shifted(:,2); 

    f_phi_pure_A = min([find(traj_phi_pure_1(:,x)<=-L,1) find(traj_phi_pure_2(:,x)<=-L,1)]);
    f_phi_pure_B = min([find(traj_phi_pure_1(:,x)>=L,1) find(traj_phi_pure_2(:,x)>=L,1)]); 
    f_phi_entangled_A = min([find(traj_phi_entangled_1(:,x)<=-L,1) find(traj_phi_entangled_2(:,x)<=-L,1)]);
    f_phi_entangled_B = min([find(traj_phi_entangled_1(:,x)>=L,1) find(traj_phi_entangled_2(:,x)>=L,1)]); 
    f_psi_pure_A = min([find(traj_psi_pure_1(:,x)<=-L,1) find(traj_psi_pure_2(:,x)<=-L,1)]);
    f_psi_pure_B = min([find(traj_psi_pure_1(:,x)>=L,1) find(traj_psi_pure_2(:,x)>=L,1)]); 
    f_psi_entangled_A = min([find(traj_psi_entangled_1(:,x)<=-L,1) find(traj_psi_entangled_2(:,x)<=-L,1)]);
    f_psi_entangled_B = min([find(traj_psi_entangled_1(:,x)>=L,1) find(traj_psi_entangled_2(:,x)>=L,1)]); 

    f_phi_pure_A_shifted = min([find(traj_phi_pure_shifted_1(:,x)<=-L,1) find(traj_phi_pure_shifted_2(:,x)<=-L,1)]);
    f_phi_pure_B_shifted = min([find(traj_phi_pure_shifted_1(:,x)>=L,1) find(traj_phi_pure_shifted_2(:,x)>=L,1)]); 
    f_phi_entangled_A_shifted = min([find(traj_phi_entangled_shifted_1(:,x)<=-L,1) find(traj_phi_entangled_shifted_2(:,x)<=-L,1)]);
    f_phi_entangled_B_shifted = min([find(traj_phi_entangled_shifted_1(:,x)>=L,1) find(traj_phi_entangled_shifted_2(:,x)>=L,1)]); 
    f_psi_pure_A_shifted = min([find(traj_psi_pure_shifted_1(:,x)<=-L,1) find(traj_psi_pure_shifted_2(:,x)<=-L,1)]);
    f_psi_pure_B_shifted = min([find(traj_psi_pure_shifted_1(:,x)>=L,1) find(traj_psi_pure_shifted_2(:,x)>=L,1)]); 
    f_psi_entangled_A_shifted = min([find(traj_psi_entangled_shifted_1(:,x)<=-L,1) find(traj_psi_entangled_shifted_2(:,x)<=-L,1)]);
    f_psi_entangled_B_shifted = min([find(traj_psi_entangled_shifted_1(:,x)>=L,1) find(traj_psi_entangled_shifted_2(:,x)>=L,1)]); 

    if ~(isempty(f_phi_pure_A) || isempty(f_phi_pure_B))
        yy_phi_pure_A(x) = tt_phi_pure(f_phi_pure_A);
        yy_phi_pure_B(x) = tt_phi_pure(f_phi_pure_B);
    end
    if ~(isempty(f_phi_entangled_A) || isempty(f_phi_entangled_B))
        yy_phi_entangled_A(x) = tt_phi_entangled(f_phi_entangled_A);
        yy_phi_entangled_B(x) = tt_phi_entangled(f_phi_entangled_B);
    end
    if ~(isempty(f_psi_pure_A) || isempty(f_psi_pure_B))
        yy_psi_pure_A(x) = tt_phi_pure(f_psi_pure_A);
        yy_psi_pure_B(x) = tt_phi_pure(f_psi_pure_B);
    end
    if ~(isempty(f_psi_entangled_A) || isempty(f_psi_entangled_B))
        yy_psi_entangled_A(x) = tt_psi_entangled(f_psi_entangled_A);
        yy_psi_entangled_B(x) = tt_psi_entangled(f_psi_entangled_B);
    end

    if ~(isempty(f_phi_pure_A_shifted) || isempty(f_phi_pure_B_shifted))
        yy_phi_pure_A_shifted(x) = tt_phi_pure(f_phi_pure_A_shifted);
        yy_phi_pure_B_shifted(x) = tt_phi_pure(f_phi_pure_B_shifted);
    end
    if ~(isempty(f_phi_entangled_A_shifted) || isempty(f_phi_entangled_B_shifted))
        yy_phi_entangled_A_shifted(x) = tt_phi_entangled(f_phi_entangled_A_shifted);
        yy_phi_entangled_B_shifted(x) = tt_phi_entangled(f_phi_entangled_B_shifted);
    end
    if ~(isempty(f_psi_pure_A_shifted) || isempty(f_psi_pure_B_shifted))
        yy_psi_pure_A_shifted(x) = tt_phi_pure(f_psi_pure_A_shifted);
        yy_psi_pure_B_shifted(x) = tt_phi_pure(f_psi_pure_B_shifted);
    end
    if ~(isempty(f_psi_entangled_A_shifted) || isempty(f_psi_entangled_B_shifted))
        yy_psi_entangled_A_shifted(x) = tt_psi_entangled(f_psi_entangled_A_shifted);
        yy_psi_entangled_B_shifted(x) = tt_psi_entangled(f_psi_entangled_B_shifted);
    end
end

space = (-L:mesh:L);

fun_phi_pure_A = zeros(time/mesh+1,2/mesh*L+1);
fun_phi_pure_B = zeros(time/mesh+1,2/mesh*L+1); 
fun_phi_entangled_A = zeros(time/mesh+1,2/mesh*L+1);
fun_phi_entangled_B = zeros(time/mesh+1,2/mesh*L+1); 
fun_psi_pure_A = zeros(time/mesh+1,2/mesh*L+1);
fun_psi_pure_B = zeros(time/mesh+1,2/mesh*L+1); 
fun_psi_entangled_A = zeros(time/mesh+1,2/mesh*L+1);
fun_psi_entangled_B = zeros(time/mesh+1,2/mesh*L+1); 

fun_phi_pure_A_shifted = zeros(time/mesh+1,2/mesh*L+1);
fun_phi_pure_B_shifted = zeros(time/mesh+1,2/mesh*L+1); 
fun_phi_entangled_A_shifted = zeros(time/mesh+1,2/mesh*L+1);
fun_phi_entangled_B_shifted = zeros(time/mesh+1,2/mesh*L+1); 
fun_psi_pure_A_shifted = zeros(time/mesh+1,2/mesh*L+1);
fun_psi_pure_B_shifted = zeros(time/mesh+1,2/mesh*L+1); 
fun_psi_entangled_A_shifted = zeros(time/mesh+1,2/mesh*L+1);
fun_psi_entangled_B_shifted = zeros(time/mesh+1,2/mesh*L+1);

mu_phi_pure_A = zeros(1,time/mesh+1);
mu_phi_pure_B = zeros(1,time/mesh+1); 
mu_phi_entangled_A = zeros(1,time/mesh+1);
mu_phi_entangled_B = zeros(1,time/mesh+1); 
mu_psi_pure_A = zeros(1,time/mesh+1);
mu_psi_pure_B = zeros(1,time/mesh+1); 
mu_psi_entangled_A = zeros(1,time/mesh+1);
mu_psi_entangled_B = zeros(1,time/mesh+1); 

mu_phi_pure_A_shifted = zeros(1,time/mesh+1);
mu_phi_pure_B_shifted = zeros(1,time/mesh+1); 
mu_phi_entangled_A_shifted = zeros(1,time/mesh+1);
mu_phi_entangled_B_shifted = zeros(1,time/mesh+1); 
mu_psi_pure_A_shifted = zeros(1,time/mesh+1);
mu_psi_pure_B_shifted = zeros(1,time/mesh+1); 
mu_psi_entangled_A_shifted = zeros(1,time/mesh+1);
mu_psi_entangled_B_shifted = zeros(1,time/mesh+1);

for t = 1:1/mesh*time+1
    t
    phiMminL1 = phiMinus(times(t),-L,theta_1,sigma_1,k_1,mu_1,omega,step);
    phiMminL2 = phiMinus(times(t),-L,theta_2,sigma_2,k_2,mu_2,omega,step);
    phiPminL1 = phiPlus(times(t),-L,theta_1,sigma_1,k_1,mu_1,omega,step);
    phiPminL2 = phiPlus(times(t),-L,theta_2,sigma_2,k_2,mu_2,omega,step);
    phiMposL1 = phiMinus(times(t),L,theta_1,sigma_1,k_1,mu_1,omega,step);
    phiMposL2 = phiMinus(times(t),L,theta_2,sigma_2,k_2,mu_2,omega,step);
    phiPposL1 = phiPlus(times(t),L,theta_1,sigma_1,k_1,mu_1,omega,step);
    phiPposL2 = phiPlus(times(t),L,theta_2,sigma_2,k_2,mu_2,omega,step);
    psiMminL1 = psiMinus(times(t),-L,theta_1,sigma_1,k_1,mu_1,omega,L,step);
    psiMminL2 = psiMinus(times(t),-L,theta_2,sigma_2,k_2,mu_2,omega,L,step);
    psiPminL1 = psiPlus(times(t),-L,theta_1,sigma_1,k_1,mu_1,omega,L,step);
    psiPminL2 = psiPlus(times(t),-L,theta_2,sigma_2,k_2,mu_2,omega,L,step);
    psiMposL1 = psiMinus(times(t),L,theta_1,sigma_1,k_1,mu_1,omega,L,step);
    psiMposL2 = psiMinus(times(t),L,theta_2,sigma_2,k_2,mu_2,omega,L,step);
    psiPposL1 = psiPlus(times(t),L,theta_1,sigma_1,k_1,mu_1,omega,L,step);
    psiPposL2 = psiPlus(times(t),L,theta_2,sigma_2,k_2,mu_2,omega,L,step);

    phiMminL11 = phiMinus(times(t),-L1,theta_1,sigma_1,k_1,mu_3,omega,step);
    phiMminL12 = phiMinus(times(t),-L1,theta_2,sigma_2,k_2,mu_4,omega,step);
    phiPminL11 = phiPlus(times(t),-L1,theta_1,sigma_1,k_1,mu_3,omega,step);
    phiPminL12 = phiPlus(times(t),-L1,theta_2,sigma_2,k_2,mu_4,omega,step);
    phiMposL11 = phiMinus(times(t),L1,theta_1,sigma_1,k_1,mu_3,omega,step);
    phiMposL12 = phiMinus(times(t),L1,theta_2,sigma_2,k_2,mu_4,omega,step);
    phiPposL11 = phiPlus(times(t),L1,theta_1,sigma_1,k_1,mu_3,omega,step);
    phiPposL12 = phiPlus(times(t),L1,theta_2,sigma_2,k_2,mu_4,omega,step);

    psiMminL11 = psiMinus(times(t),-L1,theta_1,sigma_1,k_1,mu_3,omega,L1,step);
    psiMminL12 = psiMinus(times(t),-L1,theta_2,sigma_2,k_2,mu_4,omega,L1,step);
    psiPminL11 = psiPlus(times(t),-L1,theta_1,sigma_1,k_1,mu_3,omega,L1,step);
    psiPminL12 = psiPlus(times(t),-L1,theta_2,sigma_2,k_2,mu_4,omega,L1,step);
    psiMposL11 = psiMinus(times(t),L1,theta_1,sigma_1,k_1,mu_3,omega,L1,step);
    psiMposL12 = psiMinus(times(t),L1,theta_2,sigma_2,k_2,mu_4,omega,L1,step);
    psiPposL11 = psiPlus(times(t),L1,theta_1,sigma_1,k_1,mu_3,omega,L1,step);
    psiPposL12 = psiPlus(times(t),L1,theta_2,sigma_2,k_2,mu_4,omega,L1,step);

    for s = 1:2/mesh*L+1
        phiMsp1 = phiMinus(times(t),space(s),theta_1,sigma_1,k_1,mu_1,omega,step);
        phiMsp2 = phiMinus(times(t),space(s),theta_2,sigma_2,k_2,mu_2,omega,step);

        phiPsp1 = phiPlus(times(t),space(s),theta_1,sigma_1,k_1,mu_1,omega,step);
        phiPsp2 = phiPlus(times(t),space(s),theta_2,sigma_2,k_2,mu_2,omega,step);

        psiMsp1 = psiMinus(times(t),space(s),theta_1,sigma_1,k_1,mu_1,omega,L,step);
        psiMsp2 = psiMinus(times(t),space(s),theta_2,sigma_2,k_2,mu_2,omega,L,step);

        psiPsp1 = psiPlus(times(t),space(s),theta_1,sigma_1,k_1,mu_1,omega,L,step);
        psiPsp2 = psiPlus(times(t),space(s),theta_2,sigma_2,k_2,mu_2,omega,L,step);

        psiMspL11 = psiMinus(times(t),space(s),theta_1,sigma_1,k_1,mu_3,omega,L1,step);
        psiMspL12 = psiMinus(times(t),space(s),theta_2,sigma_2,k_2,mu_4,omega,L1,step);

        psiPspL11 = psiPlus(times(t),space(s),theta_1,sigma_1,k_1,mu_3,omega,L1,step);
        psiPspL12 = psiPlus(times(t),space(s),theta_2,sigma_2,k_2,mu_4,omega,L1,step);

        arr_phi_j10_A = velocity_twobody_helper(phiMminL1,phiMminL2,phiMsp1,phiMsp2,phiPminL1,phiPminL2,phiPsp1,phiPsp2);
        arr_phi_j01_A = velocity_twobody_helper(phiMsp1,phiMsp2,phiMminL1,phiMminL2,phiPsp1,phiPsp2,phiPminL1,phiPminL2);
        arr_phi_j10_B = velocity_twobody_helper(phiMposL1,phiMposL2,phiMsp1,phiMsp2,phiPposL1,phiPposL2,phiPsp1,phiPsp2);
        arr_phi_j01_B = velocity_twobody_helper(phiMsp1,phiMsp2,phiMposL1,phiMposL2,phiPsp1,phiPsp2,phiPposL1,phiPposL2);

        arr_psi_j10_A = velocity_twobody_helper(psiMminL1,psiMminL2,psiMsp1,psiMsp2,psiPminL1,psiPminL2,psiPsp1,psiPsp2);
        arr_psi_j01_A = velocity_twobody_helper(psiMsp1,psiMsp2,psiMminL1,psiMminL2,psiPsp1,psiPsp2,psiPminL1,psiPminL2);
        arr_psi_j10_B = velocity_twobody_helper(psiMposL1,psiMposL2,psiMsp1,psiMsp2,psiPposL1,psiPposL2,psiPsp1,psiPsp2);
        arr_psi_j01_B = velocity_twobody_helper(psiMsp1,psiMsp2,psiMposL1,psiMposL2,psiPsp1,psiPsp2,psiPposL1,psiPposL2);

        arr_phi_shifted_j10_A = velocity_twobody_helper(phiMminL11,phiMminL12,phiMsp1,phiMsp2,phiPminL11,phiPminL12,phiPsp1,phiPsp2);
        arr_phi_shifted_j01_A = velocity_twobody_helper(phiMsp1,phiMsp2,phiMminL11,phiMminL12,phiPsp1,phiPsp2,phiPminL11,phiPminL12);
        arr_phi_shifted_j10_B = velocity_twobody_helper(phiMposL11,phiMposL12,phiMsp1,phiMsp2,phiPposL11,phiPposL12,phiPsp1,phiPsp2);
        arr_phi_shifted_j01_B = velocity_twobody_helper(phiMsp1,phiMsp2,phiMposL11,phiMposL12,phiPsp1,phiPsp2,phiPposL11,phiPposL12);

        arr_psi_shifted_j10_A = velocity_twobody_helper(psiMminL11,psiMminL12,psiMspL11,psiMspL12,psiPminL11,psiPminL12,psiPspL11,psiPspL12);
        arr_psi_shifted_j01_A = velocity_twobody_helper(psiMspL11,psiMspL12,psiMminL11,psiMminL12,psiPspL11,psiPspL12,psiPminL11,psiPminL12);
        arr_psi_shifted_j10_B = velocity_twobody_helper(psiMposL11,psiMposL12,psiMspL11,psiMspL12,psiPposL11,psiPposL12,psiPspL11,psiPspL12);
        arr_psi_shifted_j01_B = velocity_twobody_helper(psiMspL11,psiMspL12,psiMposL11,psiMposL12,psiPspL11,psiPspL12,psiPposL11,psiPposL12);

        fun_phi_pure_A(t,s) = (j10(arr_phi_j10_A(5),arr_phi_j10_A(6),arr_phi_j10_A(7),arr_phi_j10_A(8))+j01(arr_phi_j01_A(5),arr_phi_j01_A(6),arr_phi_j01_A(7),arr_phi_j01_A(8)));
        fun_phi_pure_B(t,s) = (j10(arr_phi_j10_B(5),arr_phi_j10_B(6),arr_phi_j10_B(7),arr_phi_j10_B(8))+j01(arr_phi_j01_B(5),arr_phi_j01_B(6),arr_phi_j01_B(7),arr_phi_j01_B(8)));
        fun_phi_entangled_A(t,s) = (j10(arr_phi_j10_A(1),arr_phi_j10_A(2),arr_phi_j10_A(3),arr_phi_j10_A(4))+j01(arr_phi_j01_A(1),arr_phi_j01_A(2),arr_phi_j01_A(3),arr_phi_j01_A(4)));
        fun_phi_entangled_B(t,s) = (j10(arr_phi_j10_B(1),arr_phi_j10_B(2),arr_phi_j10_B(3),arr_phi_j10_B(4))+j01(arr_phi_j01_B(1),arr_phi_j01_B(2),arr_phi_j01_B(3),arr_phi_j01_B(4)));
        fun_psi_pure_A(t,s) = (j10(arr_psi_j10_A(5),arr_psi_j10_A(6),arr_psi_j10_A(7),arr_psi_j10_A(8))+j01(arr_psi_j01_A(5),arr_psi_j01_A(6),arr_psi_j01_A(7),arr_psi_j01_A(8)));
        fun_psi_pure_B(t,s) = (j10(arr_psi_j10_B(5),arr_psi_j10_B(6),arr_psi_j10_B(7),arr_psi_j10_B(8))+j01(arr_psi_j01_B(5),arr_psi_j01_B(6),arr_psi_j01_B(7),arr_psi_j01_B(8)));
        fun_psi_entangled_A(t,s) = (j10(arr_psi_j10_A(1),arr_psi_j10_A(2),arr_psi_j10_A(3),arr_psi_j10_A(4))+j01(arr_psi_j01_A(1),arr_psi_j01_A(2),arr_psi_j01_A(3),arr_psi_j01_A(4)));
        fun_psi_entangled_B(t,s) = (j10(arr_psi_j10_B(1),arr_psi_j10_B(2),arr_psi_j10_B(3),arr_psi_j10_B(4))+j01(arr_psi_j01_B(1),arr_psi_j01_B(2),arr_psi_j01_B(3),arr_psi_j01_B(4)));

        fun_phi_pure_A_shifted(t,s) = (j10(arr_phi_shifted_j10_A(5),arr_phi_shifted_j10_A(6),arr_phi_shifted_j10_A(7),arr_phi_shifted_j10_A(8))+j01(arr_phi_shifted_j01_A(5),arr_phi_shifted_j01_A(6),arr_phi_shifted_j01_A(7),arr_phi_shifted_j01_A(8)));
        fun_phi_pure_B_shifted(t,s) = (j10(arr_phi_shifted_j10_B(5),arr_phi_shifted_j10_B(6),arr_phi_shifted_j10_B(7),arr_phi_shifted_j10_B(8))+j01(arr_phi_shifted_j01_B(5),arr_phi_shifted_j01_B(6),arr_phi_shifted_j01_B(7),arr_phi_shifted_j01_B(8)));
        fun_phi_entangled_A_shifted(t,s) = (j10(arr_phi_shifted_j10_A(1),arr_phi_shifted_j10_A(2),arr_phi_shifted_j10_A(3),arr_phi_shifted_j10_A(4))+j01(arr_phi_shifted_j01_A(1),arr_phi_shifted_j01_A(2),arr_phi_shifted_j01_A(3),arr_phi_shifted_j01_A(4)));
        fun_phi_entangled_B_shifted(t,s) = (j10(arr_phi_shifted_j10_B(1),arr_phi_shifted_j10_B(2),arr_phi_shifted_j10_B(3),arr_phi_shifted_j10_B(4))+j01(arr_phi_shifted_j01_B(1),arr_phi_shifted_j01_B(2),arr_phi_shifted_j01_B(3),arr_phi_shifted_j01_B(4)));
        fun_psi_pure_A_shifted(t,s) = (j10(arr_psi_shifted_j10_A(5),arr_psi_shifted_j10_A(6),arr_psi_shifted_j10_A(7),arr_psi_shifted_j10_A(8))+j01(arr_psi_shifted_j01_A(5),arr_psi_shifted_j01_A(6),arr_psi_shifted_j01_A(7),arr_psi_shifted_j01_A(8)));
        fun_psi_pure_B_shifted(t,s) = (j10(arr_psi_shifted_j10_B(5),arr_psi_shifted_j10_B(6),arr_psi_shifted_j10_B(7),arr_psi_shifted_j10_B(8))+j01(arr_psi_shifted_j01_B(5),arr_psi_shifted_j01_B(6),arr_psi_shifted_j01_B(7),arr_psi_shifted_j01_B(8)));
        fun_psi_entangled_A_shifted(t,s) = (j10(arr_psi_shifted_j10_A(1),arr_psi_shifted_j10_A(2),arr_psi_shifted_j10_A(3),arr_psi_shifted_j10_A(4))+j01(arr_psi_shifted_j01_A(1),arr_psi_shifted_j01_A(2),arr_psi_shifted_j01_A(3),arr_psi_shifted_j01_A(4)));
        fun_psi_entangled_B_shifted(t,s) = (j10(arr_psi_shifted_j10_B(1),arr_psi_shifted_j10_B(2),arr_psi_shifted_j10_B(3),arr_psi_shifted_j10_B(4))+j01(arr_psi_shifted_j01_B(1),arr_psi_shifted_j01_B(2),arr_psi_shifted_j01_B(3),arr_psi_shifted_j01_B(4)));
    end

    mu_phi_pure_A(t) = trapz(-fun_phi_pure_A(t,:))*mesh;
    mu_phi_pure_B(t) = trapz(fun_phi_pure_B(t,:))*mesh;
    mu_phi_entangled_A(t) = trapz(-fun_phi_entangled_A(t,:))*mesh;
    mu_phi_entangled_B(t) = trapz(fun_phi_entangled_B(t,:))*mesh;
    mu_psi_pure_A(t) = trapz(-fun_psi_pure_A(t,:))*mesh;
    mu_psi_pure_B(t) = trapz(fun_psi_pure_B(t,:))*mesh;
    mu_psi_entangled_A(t) = trapz(-fun_psi_entangled_A(t,:))*mesh;
    mu_psi_entangled_B(t) = trapz(fun_psi_entangled_B(t,:))*mesh;

    mu_phi_pure_A_shifted(t) = trapz(-fun_phi_pure_A_shifted(t,:))*mesh;
    mu_phi_pure_B_shifted(t) = trapz(fun_phi_pure_B_shifted(t,:))*mesh;
    mu_phi_entangled_A_shifted(t) = trapz(-fun_phi_entangled_A_shifted(t,:))*mesh;
    mu_phi_entangled_B_shifted(t) = trapz(fun_phi_entangled_B_shifted(t,:))*mesh;
    mu_psi_pure_A_shifted(t) = trapz(-fun_psi_pure_A_shifted(t,:))*mesh;
    mu_psi_pure_B_shifted(t) = trapz(fun_psi_pure_B_shifted(t,:))*mesh;
    mu_psi_entangled_A_shifted(t) = trapz(-fun_psi_entangled_A_shifted(t,:))*mesh;
    mu_psi_entangled_B_shifted(t) = trapz(fun_psi_entangled_B_shifted(t,:))*mesh;
end

% TO PLOT THE STATISTICS, GO TO PLOT_STATISTICS_TWOBODY.M