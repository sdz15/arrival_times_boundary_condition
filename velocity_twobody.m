function val = velocity_twobody(psiM_11,psiM_12,psiM_21,psiM_22,psiP_11,psiP_12,psiP_21,psiP_22,entangled) % PASS IN FUNCTIONS WITHOUT SQUARE OR ABSOLUTE VALUE
    arr = velocity_twobody_helper(psiM_11,psiM_12,psiM_21,psiM_22,psiP_11,psiP_12,psiP_21,psiP_22);
    
    [x1,x2,x3,x4]= deal(0);

    if entangled == true
        x1 = arr(1);
        x2 = arr(2);
        x3 = arr(3);
        x4 = arr(4);
    else
        x1 = arr(5);
        x2 = arr(6);
        x3 = arr(7);
        x4 = arr(8);
    end

    j_00 = j00(x1,x2,x3,x4);

    if (norm(j_00)>=1e-6)
        j_10 = j10(x1,x2,x3,x4);
        j_01 = j01(x1,x2,x3,x4);
        val = [j_10; j_01]./j_00;
        return
    end
    val = [0; 0];
end