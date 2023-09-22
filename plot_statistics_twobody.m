theta_1 = 0;
theta_2 = pi/2;
mu_1 = 0;
mu_2 = 1;
sigma_1 = .1;
sigma_2 = .1;
k_1 = 0;
k_2 = 0;
omega = 2;
L = 1;
time = 2;
N = 2;
bound = 1;
mesh = .01;

delta = .1;
L1 = L+delta/2;
mu_3 = mu_1-delta/2;
mu_4 = mu_2-delta/2;

times = (1e-1:mesh:time)+1e-6;

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

cmap = colormap;

for x=1:N
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    q0 = initvals(:,x);
    q0_shifted = initvals_shifted(:,x);

    [tt_phi_pure,qq_phi_pure] = ode45(@(t,q) velocity_twobody(phiMinus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega),phiMinus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega),phiMinus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega),phiMinus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega),phiPlus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega),phiPlus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega),phiPlus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega),phiPlus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega),false),times,q0,opts);
    [tt_phi_entangled,qq_phi_entangled] = ode45(@(t,q) velocity_twobody(phiMinus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega),phiMinus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega),phiMinus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega),phiMinus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega),phiPlus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega),phiPlus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega),phiPlus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega),phiPlus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega),true),times,q0,opts);

    [tt_psi_pure,qq_psi_pure] = ode45(@(t,q) velocity_twobody(psiMinus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega,L),psiMinus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega,L),psiMinus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega,L),psiMinus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega,L),psiPlus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega,L),psiPlus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega,L),psiPlus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega,L),psiPlus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega,L),false),times,q0,opts);
    [tt_psi_entangled,qq_psi_entangled] = ode45(@(t,q) velocity_twobody(psiMinus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega,L),psiMinus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega,L),psiMinus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega,L),psiMinus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega,L),psiPlus(t,q(1),theta_1,sigma_1,k_1,mu_1,omega,L),psiPlus(t,q(1),theta_2,sigma_2,k_2,mu_2,omega,L),psiPlus(t,q(2),theta_1,sigma_1,k_1,mu_1,omega,L),psiPlus(t,q(2),theta_2,sigma_2,k_2,mu_2,omega,L),true),times,q0,opts);

    % SAME STATISTICS ABOVE FOR WHEN L IS SHIFTED BY DELTA/2
    [tt_phi_pure_shifted,qq_phi_pure_shifted] = ode45(@(t,q) velocity_twobody(phiMinus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega),phiMinus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega),phiMinus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega),phiMinus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega),phiPlus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega),phiPlus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega),phiPlus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega),phiPlus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega),false),times,q0_shifted,opts);
    [tt_phi_entangled_shifted,qq_phi_entangled_shifted] = ode45(@(t,q) velocity_twobody(phiMinus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega),phiMinus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega),phiMinus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega),phiMinus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega),phiPlus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega),phiPlus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega),phiPlus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega),phiPlus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega),true),times,q0_shifted,opts);

    [tt_psi_pure_shifted,qq_psi_pure_shifted] = ode45(@(t,q) velocity_twobody(psiMinus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega,L),psiMinus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega,L),psiMinus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega,L),psiMinus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega,L),psiPlus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega,L),psiPlus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega,L),psiPlus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega,L),psiPlus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega,L),false),times,q0_shifted,opts);
    [tt_psi_entangled_shifted,qq_psi_entangled_shifted] = ode45(@(t,q) velocity_twobody(psiMinus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega,L),psiMinus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega,L),psiMinus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega,L),psiMinus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega,L),psiPlus(t,q(1),theta_1,sigma_1,k_1,mu_3,omega,L),psiPlus(t,q(1),theta_2,sigma_2,k_2,mu_4,omega,L),psiPlus(t,q(2),theta_1,sigma_1,k_1,mu_3,omega,L),psiPlus(t,q(2),theta_2,sigma_2,k_2,mu_4,omega,L),true),times,q0_shifted,opts);

    c = cmap(randi(size(cmap,1)),:);

    figure(1)
    plot(qq_phi_pure(:,1),tt_phi_pure,qq_phi_pure(:,2),tt_phi_pure,'Color',c);
    xlim([-L,L])
    ylim([0 time])
    drawnow
    hold on;

    figure(2)
    plot(qq_phi_entangled(:,1),tt_phi_entangled,qq_phi_entangled(:,2),tt_phi_entangled,'Color',c);
    xlim([-L,L])
    ylim([0 time])
    drawnow
    hold on;

    figure(3)
    plot(qq_psi_pure(:,1),tt_psi_pure,qq_psi_pure(:,2),tt_psi_pure,'Color',c);
    xlim([-L,L])
    ylim([0 time])
    drawnow
    hold on;

    figure(4)
    plot(qq_psi_entangled(:,1),tt_psi_entangled,qq_psi_entangled(:,2),tt_psi_entangled,'Color',c);
    xlim([-L,L])
    ylim([0 time])
    drawnow
    hold on;

    figure(5)
    plot(qq_phi_pure_shifted(:,1),tt_phi_pure_shifted,qq_phi_pure_shifted(:,2),tt_phi_pure_shifted,'Color',c);
    xlim([-L,L])
    ylim([0 time])
    drawnow
    hold on;

    figure(6)
    plot(qq_phi_entangled_shifted(:,1),tt_phi_entangled_shifted,qq_phi_entangled_shifted(:,2),tt_phi_entangled_shifted,'Color',c);
    xlim([-L,L])
    ylim([0 time])
    drawnow
    hold on;

    figure(7)
    plot(qq_psi_pure_shifted(:,1),tt_psi_pure_shifted,qq_psi_pure_shifted(:,2),tt_psi_pure_shifted,'Color',c);
    xlim([-L,L])
    ylim([0 time])
    drawnow
    hold on;

    figure(8)
    plot(qq_psi_entangled_shifted(:,1),tt_psi_entangled_shifted,qq_psi_entangled_shifted(:,2),tt_psi_entangled_shifted,'Color',c);
    xlim([-L,L])
    ylim([0 time])
    drawnow
    hold on;

    f_phi_pure_A = min([find(qq_phi_pure(:,1)<=-L,1) find(qq_phi_pure(:,2)<=-L,1)]);
    f_phi_pure_B = min([find(qq_phi_pure(:,1)>=L,1) find(qq_phi_pure(:,2)>=L,1)]); 
    f_phi_entangled_A = min([find(qq_phi_entangled(:,1)<=-L,1) find(qq_phi_entangled(:,2)<=-L,1)]);
    f_phi_entangled_B = min([find(qq_phi_entangled(:,1)>=L,1) find(qq_phi_entangled(:,2)>=L,1)]); 
    f_psi_pure_A = min([find(qq_psi_pure(:,1)<=-L,1) find(qq_psi_pure(:,2)<=-L,1)]);
    f_psi_pure_B = min([find(qq_psi_pure(:,1)>=L,1) find(qq_psi_pure(:,2)>=L,1)]); 
    f_psi_entangled_A = min([find(qq_psi_entangled(:,1)<=-L,1) find(qq_psi_entangled(:,2)<=-L,1)]);
    f_psi_entangled_B = min([find(qq_psi_entangled(:,1)>=L,1) find(qq_psi_entangled(:,2)>=L,1)]); 

    f_phi_pure_A_shifted = min([find(qq_phi_pure_shifted(:,1)<=-L,1) find(qq_phi_pure_shifted(:,2)<=-L,1)]);
    f_phi_pure_B_shifted = min([find(qq_phi_pure_shifted(:,1)>=L,1) find(qq_phi_pure_shifted(:,2)>=L,1)]); 
    f_phi_entangled_A_shifted = min([find(qq_phi_entangled_shifted(:,1)<=-L,1) find(qq_phi_entangled_shifted(:,2)<=-L,1)]);
    f_phi_entangled_B_shifted = min([find(qq_phi_entangled_shifted(:,1)>=L,1) find(qq_phi_entangled_shifted(:,2)>=L,1)]); 
    f_psi_pure_A_shifted = min([find(qq_psi_pure_shifted(:,1)<=-L,1) find(qq_psi_pure_shifted(:,2)<=-L,1)]);
    f_psi_pure_B_shifted = min([find(qq_psi_pure_shifted(:,1)>=L,1) find(qq_psi_pure_shifted(:,2)>=L,1)]); 
    f_psi_entangled_A_shifted = min([find(qq_psi_entangled_shifted(:,1)<=-L,1) find(qq_psi_entangled_shifted(:,2)<=-L,1)]);
    f_psi_entangled_B_shifted = min([find(qq_psi_entangled_shifted(:,1)>=L,1) find(qq_psi_entangled_shifted(:,2)>=L,1)]); 

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

txt = {strcat('theta\_ph=',string(theta_1)),strcat('theta\_el=',string(theta_2)),strcat('mu\_ph=',string(mu_1)),strcat('mu\_el=',string(mu_2)),strcat('sigma\_ph=',string(sigma_1)),strcat('sigma\_el=',string(sigma_2)),strcat('k\_ph=',string(k_1)),strcat('k\_el=',string(k_2)),strcat('omega=',string(omega))};

figure(1)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for pure product without boundary condition');
text(mu_2+.5,.5,txt);
hold off;

figure(2)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for entangled pair without boundary condition');
text(mu_2+.5,.5,txt);
hold off;

figure(3)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for pure product with boundary condition');
text(mu_2+.5,.5,txt);
hold off;

figure(4)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for entangled pair with boundary condition');
text(mu_2+.5,.5,txt);
hold off;

figure(5)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for pure product without boundary condition and shifted Bob');
text(mu_2+.5,.5,txt);
hold off;

figure(6)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for entangled pair without boundary condition and shifted Bob');
text(mu_2+.5,.5,txt);
hold off;

figure(7)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for pure product with boundary condition and shifted Bob');
text(mu_2+.5,.5,txt);
hold off;

figure(8)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for entangled pair with boundary condition and shifted Bob');
text(mu_2+.5,.5,txt);
hold off;

figure(9)
scatter(yy_phi_pure_A,yy_phi_pure_B);
xlabel('Arrival time at Alice','FontSize',20);
ylabel('Arrival time at Bob','FontSize',20);
title('Arrival times for pure product without boundary condition');
xlim([-L,L])
ylim([0 time])
text(mu_2+.5,.5,txt);

figure(10)
scatter(yy_phi_entangled_A,yy_phi_entangled_B);
xlabel('Arrival time at Alice','FontSize',20);
ylabel('Arrival time at Bob','FontSize',20);
title('Arrival times for entangled pair without boundary condition');
xlim([-L,L])
ylim([0 time])
text(mu_2+.5,.5,txt);

figure(11)
scatter(yy_psi_pure_A,yy_psi_pure_B);
xlabel('Arrival time at Alice','FontSize',20);
ylabel('Arrival time at Bob','FontSize',20);
title('Arrival times for pure product with boundary condition');
xlim([-L,L])
ylim([0 time])
text(mu_2+.5,.5,txt);

figure(12)
scatter(yy_psi_entangled_A,yy_psi_entangled_B);
xlabel('Arrival time at Alice','FontSize',20);
ylabel('Arrival time at Bob','FontSize',20);
title('Arrival times for entangled pair with boundary condition');
xlim([-L,L])
ylim([0 time])
text(mu_2+.5,.5,txt);

figure(13)
scatter(yy_phi_pure_A_shifted,yy_phi_pure_B_shifted);
xlabel('Arrival time at Alice','FontSize',20);
ylabel('Arrival time at Bob','FontSize',20);
title('Arrival times for pure product without boundary condition and shifted Bob');
xlim([-L,L])
ylim([0 time])
text(mu_2+.5,.5,txt);

figure(14)
scatter(yy_phi_entangled_A_shifted,yy_phi_entangled_B_shifted);
xlabel('Arrival time at Alice','FontSize',20);
ylabel('Arrival time at Bob','FontSize',20);
title('Arrival times for entangled pair without boundary condition and shifted Bob');
xlim([-L,L])
ylim([0 time])
text(mu_2+.5,.5,txt);

figure(15)
scatter(yy_psi_pure_A_shifted,yy_psi_pure_B_shifted);
xlabel('Arrival time at Alice','FontSize',20);
ylabel('Arrival time at Bob','FontSize',20);
title('Arrival times for pure product with boundary condition and shifted Bob');
xlim([-L,L])
ylim([0 time])
text(mu_2+.5,.5,txt);

figure(16)
scatter(yy_psi_entangled_A_shifted,yy_psi_entangled_B_shifted);
xlabel('Arrival time at Alice','FontSize',20);
ylabel('Arrival time at Bob','FontSize',20);
title('Arrival times for entangled pair with boundary condition and shifted Bob');
xlim([-L,L])
ylim([0 time])
text(mu_2+.5,.5,txt);

figure(17)
cdfplot(yy_phi_pure_A);
hold on
cdfplot(yy_phi_pure_B);
hold on
cdfplot(yy_phi_pure_A_shifted);
hold on
cdfplot(yy_phi_pure_B_shifted);
hold on
title('Cumulative distribution of arrival time for pure product without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
text(time-.5,.5,txt);

figure(18)
cdfplot(yy_phi_entangled_A);
hold on
cdfplot(yy_phi_entangled_B);
hold on
cdfplot(yy_phi_entangled_A_shifted);
hold on
cdfplot(yy_phi_entangled_B_shifted);
hold on
title('Cumulative distribution of arrival time for entangled pair without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
text(time-.5,.5,txt);

figure(19)
cdfplot(yy_psi_pure_A);
hold on
cdfplot(yy_psi_pure_B);
hold on
cdfplot(yy_psi_pure_A_shifted);
hold on
cdfplot(yy_psi_pure_B_shifted);
hold on
title('Cumulative distribution of arrival time for pure product with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
text(time-.5,.5,txt);

figure(20)
cdfplot(yy_psi_entangled_A);
hold on
cdfplot(yy_psi_entangled_B);
hold on
cdfplot(yy_psi_entangled_A_shifted);
hold on
cdfplot(yy_psi_entangled_B_shifted);
hold on
title('Cumulative distribution of arrival time for entangled pair with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
text(time-.5,.5,txt);

figure(21)
histogram(yy_phi_pure_A);
hold on
histogram(yy_phi_pure_B);
hold on
histogram(yy_phi_pure_A_shifted);
hold on
histogram(yy_phi_pure_B_shifted);
hold on
title('Histogram of arrival time for pure product without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
text(time-.5,.5,txt);

figure(22)
histogram(yy_phi_entangled_A);
hold on
histogram(yy_phi_entangled_B);
hold on
histogram(yy_phi_entangled_A_shifted);
hold on
histogram(yy_phi_entangled_B_shifted);
hold on
title('Histogram of arrival time for entangled pair without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
text(time-.5,.5,txt);

figure(23)
histogram(yy_psi_pure_A);
hold on
histogram(yy_psi_pure_B);
hold on
histogram(yy_psi_pure_A_shifted);
hold on
histogram(yy_psi_pure_B_shifted);
hold on
title('Histogram of arrival time for pure product with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
text(time-.5,.5,txt);

figure(24)
histogram(yy_psi_entangled_A);
hold on
histogram(yy_psi_entangled_B);
hold on
histogram(yy_psi_entangled_A_shifted);
hold on
histogram(yy_psi_entangled_B_shifted);
hold on
title('Histogram of arrival time for entangled pair with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
text(time-.5,.5,txt);

space = (-L:mesh:L);
times = (0:mesh:time);

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
    for s = 1:2/mesh*L+1
        phiMminL1 = phiMinus(times(t),-L,theta_1,sigma_1,k_1,mu_1,omega);
        phiMminL2 = phiMinus(times(t),-L,theta_2,sigma_2,k_2,mu_2,omega);
        phiMsp1 = phiMinus(times(t),space(s),theta_1,sigma_1,k_1,mu_1,omega);
        phiMsp2 = phiMinus(times(t),space(s),theta_2,sigma_2,k_2,mu_2,omega);
        phiPminL1 = phiPlus(times(t),-L,theta_1,sigma_1,k_1,mu_1,omega);
        phiPminL2 = phiPlus(times(t),-L,theta_2,sigma_2,k_2,mu_2,omega);
        phiPsp1 = phiPlus(times(t),space(s),theta_1,sigma_1,k_1,mu_1,omega);
        phiPsp2 = phiPlus(times(t),space(s),theta_2,sigma_2,k_2,mu_2,omega);
        phiMposL1 = phiMinus(times(t),L,theta_1,sigma_1,k_1,mu_1,omega);
        phiMposL2 = phiMinus(times(t),L,theta_2,sigma_2,k_2,mu_2,omega);
        phiPposL1 = phiPlus(times(t),L,theta_1,sigma_1,k_1,mu_1,omega);
        phiPposL2 = phiPlus(times(t),L,theta_2,sigma_2,k_2,mu_2,omega);

        psiMminL1 = psiMinus(times(t),-L,theta_1,sigma_1,k_1,mu_1,omega,L);
        psiMminL2 = psiMinus(times(t),-L,theta_2,sigma_2,k_2,mu_2,omega,L);
        psiMsp1 = psiMinus(times(t),space(s),theta_1,sigma_1,k_1,mu_1,omega,L);
        psiMsp2 = psiMinus(times(t),space(s),theta_2,sigma_2,k_2,mu_2,omega,L);
        psiPminL1 = psiPlus(times(t),-L,theta_1,sigma_1,k_1,mu_1,omega,L);
        psiPminL2 = psiPlus(times(t),-L,theta_2,sigma_2,k_2,mu_2,omega,L);
        psiPsp1 = psiPlus(times(t),space(s),theta_1,sigma_1,k_1,mu_1,omega,L);
        psiPsp2 = psiPlus(times(t),space(s),theta_2,sigma_2,k_2,mu_2,omega,L);
        psiMposL1 = psiMinus(times(t),L,theta_1,sigma_1,k_1,mu_1,omega,L);
        psiMposL2 = psiMinus(times(t),L,theta_2,sigma_2,k_2,mu_2,omega,L);
        psiPposL1 = psiPlus(times(t),L,theta_1,sigma_1,k_1,mu_1,omega,L);
        psiPposL2 = psiPlus(times(t),L,theta_2,sigma_2,k_2,mu_2,omega,L);

        phiMminL11 = phiMinus(times(t),-L1,theta_1,sigma_1,k_1,mu_3,omega);
        phiMminL12 = phiMinus(times(t),-L1,theta_2,sigma_2,k_2,mu_4,omega);
        phiPminL11 = phiPlus(times(t),-L1,theta_1,sigma_1,k_1,mu_3,omega);
        phiPminL12 = phiPlus(times(t),-L1,theta_2,sigma_2,k_2,mu_4,omega);
        phiMposL11 = phiMinus(times(t),L1,theta_1,sigma_1,k_1,mu_3,omega);
        phiMposL12 = phiMinus(times(t),L1,theta_2,sigma_2,k_2,mu_4,omega);
        phiPposL11 = phiPlus(times(t),L1,theta_1,sigma_1,k_1,mu_3,omega);
        phiPposL12 = phiPlus(times(t),L1,theta_2,sigma_2,k_2,mu_4,omega);
        
        psiMminL11 = psiMinus(times(t),-L1,theta_1,sigma_1,k_1,mu_3,omega,L1);
        psiMminL12 = psiMinus(times(t),-L1,theta_2,sigma_2,k_2,mu_4,omega,L1);
        psiMspL11 = psiMinus(times(t),space(s),theta_1,sigma_1,k_1,mu_3,omega,L1);
        psiMspL12 = psiMinus(times(t),space(s),theta_2,sigma_2,k_2,mu_4,omega,L1);
        psiPminL11 = psiPlus(times(t),-L1,theta_1,sigma_1,k_1,mu_3,omega,L1);
        psiPminL12 = psiPlus(times(t),-L1,theta_2,sigma_2,k_2,mu_4,omega,L1);
        psiPspL11 = psiPlus(times(t),space(s),theta_1,sigma_1,k_1,mu_3,omega,L1);
        psiPspL12 = psiPlus(times(t),space(s),theta_2,sigma_2,k_2,mu_4,omega,L1);
        psiMposL11 = psiMinus(times(t),L1,theta_1,sigma_1,k_1,mu_3,omega,L1);
        psiMposL12 = psiMinus(times(t),L1,theta_2,sigma_2,k_2,mu_4,omega,L1);
        psiPposL11 = psiPlus(times(t),L1,theta_1,sigma_1,k_1,mu_3,omega,L1);
        psiPposL12 = psiPlus(times(t),L1,theta_2,sigma_2,k_2,mu_4,omega,L1);

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

        fun_phi_pure_A(t,s) = j10(arr_phi_j10_A(5),arr_phi_j10_A(6),arr_phi_j10_A(7),arr_phi_j10_A(8))+j01(arr_phi_j01_A(5),arr_phi_j01_A(6),arr_phi_j01_A(7),arr_phi_j01_A(8));
        fun_phi_pure_B(t,s) = j10(arr_phi_j10_B(5),arr_phi_j10_B(6),arr_phi_j10_B(7),arr_phi_j10_B(8))+j01(arr_phi_j01_B(5),arr_phi_j01_B(6),arr_phi_j01_B(7),arr_phi_j01_B(8));
        fun_phi_entangled_A(t,s) = j10(arr_phi_j10_A(1),arr_phi_j10_A(2),arr_phi_j10_A(3),arr_phi_j10_A(4))+j01(arr_phi_j01_A(1),arr_phi_j01_A(2),arr_phi_j01_A(3),arr_phi_j01_A(4));
        fun_phi_entangled_B(t,s) = j10(arr_phi_j10_B(1),arr_phi_j10_B(2),arr_phi_j10_B(3),arr_phi_j10_B(4))+j01(arr_phi_j01_B(1),arr_phi_j01_B(2),arr_phi_j01_B(3),arr_phi_j01_B(4));
        fun_psi_pure_A(t,s) = j10(arr_psi_j10_A(5),arr_psi_j10_A(6),arr_psi_j10_A(7),arr_psi_j10_A(8))+j01(arr_psi_j01_A(5),arr_psi_j01_A(6),arr_psi_j01_A(7),arr_psi_j01_A(8));
        fun_psi_pure_B(t,s) = j10(arr_psi_j10_B(5),arr_psi_j10_B(6),arr_psi_j10_B(7),arr_psi_j10_B(8))+j01(arr_psi_j01_B(5),arr_psi_j01_B(6),arr_psi_j01_B(7),arr_psi_j01_B(8));
        fun_psi_entangled_A(t,s) = j10(arr_psi_j10_A(1),arr_psi_j10_A(2),arr_psi_j10_A(3),arr_psi_j10_A(4))+j01(arr_psi_j01_A(1),arr_psi_j01_A(2),arr_psi_j01_A(3),arr_psi_j01_A(4));
        fun_psi_entangled_B(t,s) = j10(arr_psi_j10_B(1),arr_psi_j10_B(2),arr_psi_j10_B(3),arr_psi_j10_B(4))+j01(arr_psi_j01_B(1),arr_psi_j01_B(2),arr_psi_j01_B(3),arr_psi_j01_B(4));
        
        fun_phi_pure_A_shifted(t,s) = j10(arr_phi_shifted_j10_A(5),arr_phi_shifted_j10_A(6),arr_phi_shifted_j10_A(7),arr_phi_shifted_j10_A(8))+j01(arr_phi_shifted_j01_A(5),arr_phi_shifted_j01_A(6),arr_phi_shifted_j01_A(7),arr_phi_shifted_j01_A(8));
        fun_phi_pure_B_shifted(t,s) = j10(arr_phi_shifted_j10_B(5),arr_phi_shifted_j10_B(6),arr_phi_shifted_j10_B(7),arr_phi_shifted_j10_B(8))+j01(arr_phi_shifted_j01_B(5),arr_phi_shifted_j01_B(6),arr_phi_shifted_j01_B(7),arr_phi_shifted_j01_B(8));
        fun_phi_entangled_A_shifted(t,s) = j10(arr_phi_shifted_j10_A(1),arr_phi_shifted_j10_A(2),arr_phi_shifted_j10_A(3),arr_phi_shifted_j10_A(4))+j01(arr_phi_shifted_j01_A(1),arr_phi_shifted_j01_A(2),arr_phi_shifted_j01_A(3),arr_phi_shifted_j01_A(4));
        fun_phi_entangled_B_shifted(t,s) = j10(arr_phi_shifted_j10_B(1),arr_phi_shifted_j10_B(2),arr_phi_shifted_j10_B(3),arr_phi_shifted_j10_B(4))+j01(arr_phi_shifted_j01_B(1),arr_phi_shifted_j01_B(2),arr_phi_shifted_j01_B(3),arr_phi_shifted_j01_B(4));
        fun_psi_pure_A_shifted(t,s) = j10(arr_psi_shifted_j10_A(5),arr_psi_shifted_j10_A(6),arr_psi_shifted_j10_A(7),arr_psi_shifted_j10_A(8))+j01(arr_psi_shifted_j01_A(5),arr_psi_shifted_j01_A(6),arr_psi_shifted_j01_A(7),arr_psi_shifted_j01_A(8));
        fun_psi_pure_B_shifted(t,s) = j10(arr_psi_shifted_j10_B(5),arr_psi_shifted_j10_B(6),arr_psi_shifted_j10_B(7),arr_psi_shifted_j10_B(8))+j01(arr_psi_shifted_j01_B(5),arr_psi_shifted_j01_B(6),arr_psi_shifted_j01_B(7),arr_psi_shifted_j01_B(8));
        fun_psi_entangled_A_shifted(t,s) = j10(arr_psi_shifted_j10_A(1),arr_psi_shifted_j10_A(2),arr_psi_shifted_j10_A(3),arr_psi_shifted_j10_A(4))+j01(arr_psi_shifted_j01_A(1),arr_psi_shifted_j01_A(2),arr_psi_shifted_j01_A(3),arr_psi_shifted_j01_A(4));
        fun_psi_entangled_B_shifted(t,s) = j10(arr_psi_shifted_j10_B(1),arr_psi_shifted_j10_B(2),arr_psi_shifted_j10_B(3),arr_psi_shifted_j10_B(4))+j01(arr_psi_shifted_j01_B(1),arr_psi_shifted_j01_B(2),arr_psi_shifted_j01_B(3),arr_psi_shifted_j01_B(4));
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

figure(25)
plot(times,mu_phi_pure_A);
hold on
plot(times,mu_phi_pure_B);
hold on
plot(times,mu_phi_pure_A_shifted);
hold on
plot(times,mu_phi_pure_B_shifted);
hold on
title('Probability density function of arrival times for pure product without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
text(time-.5,.5,txt);

figure(26)
plot(times,mu_phi_entangled_A);
hold on
plot(times,mu_phi_entangled_B);
hold on
plot(times,mu_phi_entangled_A_shifted);
hold on
plot(times,mu_phi_entangled_B_shifted);
hold on
title('Probability density function of arrival times for entangled pair without boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
text(time-.5,.5,txt);

figure(27)
plot(times,mu_phi_pure_A);
hold on
plot(times,mu_phi_pure_B);
hold on
plot(times,mu_phi_pure_A_shifted);
hold on
plot(times,mu_phi_pure_B_shifted);
hold on
title('Probability density function of arrival times for pure product with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
text(time-.5,.5,txt);

figure(28)
plot(times,mu_phi_entangled_A);
hold on
plot(times,mu_phi_entangled_B);
hold on
plot(times,mu_phi_entangled_A_shifted);
hold on
plot(times,mu_phi_entangled_B_shifted);
hold on
title('Probability density function of arrival times for entangled pair with boundary condition');
legend('Alice','Bob','Alice shifted','Bob shifted','Location','southeast','fontsize',12);
xlim([0 time])
text(time-.5,.5,txt);

save('two_body_variables.mat')