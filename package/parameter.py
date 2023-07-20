

#######各温度のparameter########

def param86K():
    M1s    = 203.78362			
    M2s    = 184.6378
    gamma1 = 1.76 * 10**7 #rad/(s・Oe)
    gamma2 = 1.76 * 10**7 #rad/(s・Oe)
    lamb   = 1550 #G/(emu/cm^3)
    H1u    = 2000 #Oe
    H2u    = 2000 #Oe
    Ho     = 0 #Oe
    alpha1p =  0.09
    alpha1m = -0.09 
    alpha2p = -0.09
    alpha2m = +0.09 
    param86Kp = [alpha1p,alpha2p,M1s,M2s,gamma1,gamma2,lamb,H1u,H2u,Ho]
    param86Km = [alpha1m,alpha2m,M1s,M2s,gamma1,gamma2,lamb,H1u,H2u,Ho]
    return param86Kp,param86Km

def param66K():
    M1s    = 231.0191	
    M2s    = 185.22546
    gamma1 = 1.76 * 10**7 #rad/(s・Oe)
    gamma2 = 1.76 * 10**7 #rad/(s・Oe)
    lamb   = 1550 #G/(emu/cm^3)
    H1u    = 2000 #Oe
    H2u    = 2000 #Oe
    Ho     = 0 #Oe
    alpha1p = 0.11
    alpha2p = -0.11
    alpha1m = -0.06 
    alpha2m = +0.06  
    param86Kp = [alpha1p,alpha2p,M1s,M2s,gamma1,gamma2,lamb,H1u,H2u,Ho]
    param86Km = [alpha1m,alpha2m,M1s,M2s,gamma1,gamma2,lamb,H1u,H2u,Ho]
    return param86Kp,param86Km



