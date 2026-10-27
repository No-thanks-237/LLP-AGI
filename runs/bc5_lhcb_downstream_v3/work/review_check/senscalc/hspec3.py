
import numpy as np
p=r'E:/AAAwork/LLP/code/SensCalc-main/spectra/SM particles/DoubleDistr_LHC_h.txt'
d=np.loadtxt(p)
th=np.unique(d[:,0]); E=np.unique(d[:,1]); V=d[:,2].reshape(len(th),len(E))
def integ2(th,E,V,weight=None):
    out=0.0
    for i in range(len(E)-1):
        dE=E[i+1]-E[i]
        col=0.5*(V[:,i]+V[:,i+1])
        if weight is not None: col=col*weight
        out+=np.trapezoid(col, th)*dE
    return out
print('(a) int f dth dE          = %.6g'%integ2(th,E,V))
print('(b) int f 2pi sinth dth dE= %.6g'%integ2(th,E,V,2*np.pi*np.sin(th)))
lt=np.log10(th); lE=np.log10(E)
print('(c) int f dlogth dlogE    = %.6g'%integ2(lt,lE,V))
print('(d) int f 2pi sinth dth dE with dE as GeV  = %.6g'%integ2(th,E,V,2*np.pi*np.sin(th)))
# E spectrum: dN/dE = int f dtheta
gE=np.array([np.trapezoid(V[:,j],th) for j in range(len(E))])
print('E integral(a) = %.4g'%np.trapezoid(gE,E))
print('E grid 20 values:', ['%.1f'%x for x in E[:20]])
print('dN/dE at E=200,500,1000,2000:', ['%.3e'%gE[np.argmin(np.abs(E-x))] for x in [200,500,1000,2000]])
# theta: fraction of int f dE carried below theta=0.264 rad (eta>2)
gT=np.array([np.trapezoid(V[i,:],E) for i in range(len(th))])
c=0
for i in range(len(th)-1):
    if th[i]<0.264:
        c+=np.trapezoid(gT[i:i+2],th[i:i+2])
print('fraction of (a) with theta<0.264 rad (eta>2) = %.4f'%(c/integ2(th,E,V)))
