
import numpy as np, gzip, os
p=r'E:/AAAwork/LLP/code/SensCalc-main/spectra/SM particles/DoubleDistr_LHC_h.txt'
d=np.loadtxt(p)
th=np.unique(d[:,0]); E=np.unique(d[:,1])
V=d[:,2].reshape(len(th),len(E))
# integral over theta,E with trapezoid
def trapz_2d(th,E,V):
    out=0.0
    for i in range(len(E)-1):
        dE=E[i+1]-E[i]
        col=0.5*(V[:,i]+V[:,i+1])
        out+=np.trapezoid(col,th)*dE
    return out
tot=trapz_2d(th,E,V)
print('integral f dtheta dE = %.6f'%tot)
eta=-np.log(np.tan(th/2))
# dN/deta = int f dE * dtheta/deta
g=np.array([np.trapezoid(V[i,:],E) for i in range(len(th))])  # int dE f(theta,E)
dthdeta=2*np.exp(-eta)/(1+np.exp(-2*eta))
dNdeta=g*np.abs(dthdeta)
# integrate over eta
print('int dNdeta deta = %.4f'%np.trapezoid(dNdeta,eta))
for lo,hi in [(2,5),(2,4),(3,5),(2,3)]:
    m=(eta>=lo)&(eta<=hi)
    print('fraction eta in [%g,%g] = %.4e'%(lo,hi,np.trapezoid(dNdeta[m],eta[m])))
# show dN/deta at eta=2,3,4,5
for e0 in [1,2,2.5,3,4,5,6]:
    i=np.argmin(np.abs(eta-e0))
    print('  eta=%.2f: dN/deta=%.4e'%(eta[i],dNdeta[i]))
print()
print('E grid weights near top: E=%.1f int f dtheta=%.3e'%(E[-1], np.trapezoid(V[:,-1],th)))
print('E=%.1f: %.3e, E=%.1f: %.3e'%(E[-2],np.trapezoid(V[:,-2],th),E[30],np.trapezoid(V[:,30],th)))
