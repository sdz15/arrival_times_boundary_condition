function val = psiPlus(t,s,theta,sigma,k,mu,omega,L)

if t>=L-s
    f = -phiPlus(t-L+s,L,theta,sigma,k,mu,omega);
else
    f = 0;
end

n = 0;
chisum = 0;

if t > 0
    while (t>(2*n*L))
        chisum = chisum + chi(t,s,theta,sigma,k,mu,omega,L,1,n);
        n = n+1;
    end
end

val = phiPlus(t,s,theta,sigma,k,mu,omega)+f+chisum; 