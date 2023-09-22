function val = besselfun(x,omega,index)
    if (norm(x)<=1e-6)
        val = (omega^index)/(factorial(index)*2^index);
        return
    end
    
    val = besselj(index,omega*x)./(x.^index);
end