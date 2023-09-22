function val = phiMinus(t,s,theta,sigma,k,mu,omega)
    x1 = psiMinusR(s-t,theta,sigma,k,mu);
    x2 = -(omega/2)*integral(@(sig) besselfun(radical(t,s,sig),omega,1).*(t+s-sig).*psiMinusR(sig,theta,sigma,k,mu),s-t,s+t);
    x3 = -(1i*omega/2)*integral(@(sig) besselfun(radical(t,s,sig),omega,0).*psiPlusR(sig,theta,sigma,k,mu),s-t,s+t);
    
    val = x1+x2+x3;
end