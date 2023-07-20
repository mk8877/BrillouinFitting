from numpy import *

#####磁気感受率の関数を定義#####
def Chi1xx(omega,alpha1,alpha2,M1s,M2s,gamma1,gamma2,lamb,H1u,H2u,Ho):
    H1 = Ho + H1u + lamb*M2s +2500- 4*pi*(M1s-M2s)
    H2 = Ho - H2u - lamb*M1s +2500- 4*pi*(M1s-M2s)
    d_Re = (-(1+alpha2**2)*gamma1*H1 + (1+alpha1*alpha2)*gamma2*lamb*M2s + 2*alpha1*alpha2*gamma2*H2)*omega**2 + gamma1*gamma2**2*(H1*H2**2 + H1*H2*lamb*M2s + H2*lamb**2*M1s*M2s + lamb**3*M1s*M2s**2)
    d_Im = (1+alpha2**2)*alpha1*omega**3 + (2*gamma1*gamma2 *H1*H2*alpha2 + gamma1*gamma2 *H1*lamb*M2s*alpha2 + gamma1*gamma2*lamb**2*M1s*M2s*alpha2 - gamma2**2*H2**2*alpha1 - gamma2**2*H2*lamb*M2s*alpha1)*omega
    n = ((omega+gamma1*H1+1j*alpha1*omega)*(omega+gamma2*H2-1j*alpha2*omega) + gamma1*gamma2*lamb**2*M1s*M2s) * ((omega-gamma1*H1-1j*alpha1*omega)*(omega-gamma2*H2+1j*alpha2*omega) + gamma1*gamma2*lamb**2*M1s*M2s)
    return (d_Re - 1j*d_Im)/n*gamma1*M1s

def Chi1xy(omega,alpha1,alpha2,M1s,M2s,gamma1,gamma2,lamb,H1u,H2u,Ho):
    H1 = Ho + H1u + lamb*M2s +2500- 4*pi*(M1s-M2s)
    H2 = Ho - H2u - lamb*M1s +2500- 4*pi*(M1s-M2s)
    d_Re = ((alpha1-alpha2)*gamma2*lamb*M2s - 2*alpha2*gamma2*H2)*omega**2
    d_Im = (1+alpha2**2)*omega**3 - (gamma2**2 *H2**2 + gamma2**2 *H2*lamb*M2s + gamma1*gamma2*H1*lamb*M2s - gamma1*gamma2*lamb**2*M1s*M2s)*omega
    n = ((omega+gamma1*H1+1j*alpha1*omega)*(omega+gamma2*H2-1j*alpha2*omega) + gamma1*gamma2*lamb**2*M1s*M2s) * ((omega-gamma1*H1-1j*alpha1*omega)*(omega-gamma2*H2+1j*alpha2*omega) + gamma1*gamma2*lamb**2*M1s*M2s)
    return (d_Re + 1j*d_Im)/n*gamma1*M1s

def Chi2xx(omega,alpha1,alpha2,M1s,M2s,gamma1,gamma2,lamb,H1u,H2u,Ho):
    H1 = Ho + H1u + lamb*M2s +2500- 4*pi*(M1s-M2s)
    H2 = Ho - H2u - lamb*M1s +2500- 4*pi*(M1s-M2s)
    d_Re = ((1+alpha1**2)*gamma2*H2 + (1+alpha1*alpha2)*gamma1*lamb*M1s - 2*alpha1*alpha2*gamma1*H1)*omega**2 + gamma1**2*gamma2*(-H1**2*H2 + H1*H2*lamb*M1s - H1*lamb**2*M1s*M2s + lamb**3*M1s**2*M2s)
    d_Im = (1+alpha1**2)*alpha2*omega**3 + (2*gamma1*gamma2 *H1*H2*alpha1 - gamma1*gamma2 *H2*lamb*M1s*alpha1 + gamma1*gamma2*lamb**2*M1s*M2s*alpha1 - gamma1**2*H1**2*alpha2 + gamma1**2*H1*lamb*M1s*alpha2)*omega
    n = ((omega+gamma1*H1+1j*alpha1*omega)*(omega+gamma2*H2-1j*alpha2*omega) + gamma1*gamma2*lamb**2*M1s*M2s) * ((omega-gamma1*H1-1j*alpha1*omega)*(omega-gamma2*H2+1j*alpha2*omega) + gamma1*gamma2*lamb**2*M1s*M2s)
    return (d_Re - 1j*d_Im)/n*gamma2*M2s

def Chi2xy(omega,alpha1,alpha2,M1s,M2s,gamma1,gamma2,lamb,H1u,H2u,Ho):
    H1 = Ho + H1u + lamb*M2s +2500- 4*pi*(M1s-M2s)
    H2 = Ho - H2u - lamb*M1s +2500- 4*pi*(M1s-M2s)
    d_Re = ((alpha1-alpha2)*gamma1*lamb*M1s - 2*alpha1*gamma1*H1)*omega**2
    d_Im = (1+alpha1**2)*omega**3 - (gamma1**2 *H1**2 - gamma1**2 *H1*lamb*M1s - gamma1*gamma2*H2*lamb*M1s - gamma1*gamma2*lamb**2*M1s*M2s)*omega
    n = ((omega+gamma1*H1+1j*alpha1*omega)*(omega+gamma2*H2-1j*alpha2*omega) + gamma1*gamma2*lamb**2*M1s*M2s) * ((omega-gamma1*H1-1j*alpha1*omega)*(omega-gamma2*H2+1j*alpha2*omega) + gamma1*gamma2*lamb**2*M1s*M2s)
    return (d_Re - 1j*d_Im)/n*gamma2*M2s