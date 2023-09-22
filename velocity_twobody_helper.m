function val = velocity_twobody_helper(psiM_11,psiM_12,psiM_21,psiM_22,psiP_11,psiP_12,psiP_21,psiP_22) % PASS IN FUNCTIONS WITHOUT SQUARE OR ABSOLUTE VALUE
    a = psiM_11.*psiM_22;
    b = psiM_11.*psiP_22;
    c = psiP_11.*psiM_22;
    d = psiP_11.*psiP_22;
    e = psiM_12.*psiM_21;
    f = psiM_12.*psiP_21;
    g = psiP_12.*psiM_21;
    h = psiP_12.*psiP_21;
    
    x1 = abs(a-e).^2;
    x2 = abs(b-f).^2;
    x3 = abs(c-g).^2;
    x4 = abs(d-h).^2;
    
    x5 = abs(a).^2;
    x6 = abs(b).^2;
    x7 = abs(c).^2;
    x8 = abs(d).^2;

    val = [x1 x2 x3 x4 x5 x6 x7 x8];
    % FIRST FOUR VALUES ARE ENTANGLED, NEXT FOUR ARE PURE
end