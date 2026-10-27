import numpy as np
d=r'E:/AAAwork/LLP-AGI/runs/bc5_lhcb_downstream_v3/work/physics_review_20261025/'
a=np.loadtxt(d+'reference_2023_Higgs.txt')
th=np.unique(a[:,0]); E=np.unique(a[:,1]); V=a[:,2].reshape(len(th),len(E))
print('2023: ntheta=%d nE=%d ; E %.1f-%.1f ; theta %.2e..%.4f'%(len(th),len(E),E[0],E[-1],th[0],th[-1]))
eta=-np.log(np.tan(th/2)); o=np.argsort(eta); eta=eta[o]; V=V[o,:]
g=np.array([np.trapezoid(V[i,:],E) for i in range(len(th))])
dn=g*np.abs(2*np.exp(-eta)/(1+np.exp(-2*eta))); dn=dn/np.trapezoid(dn,eta)
for e0 in [0,1,2,3,4,5]:
    i=np.argmin(np.abs(eta-e0)); print('  2023 dN/deta(%.0f)=%.4f'%(e0,dn[i]))
for lo,hi in [(2,5)]:
    m=(eta>=lo)&(eta<=hi); print('  2023 frac[%g,%g]=%.4f'%(lo,hi,np.trapezoid(dn[m],eta[m])))
print()
c=np.genfromtxt(d+'reference_2023_ctau.dat', invalid_raise=False)
ok=np.isfinite(c[:,1]); ms=c[ok,0]; ct=c[ok,1]
print('ctau table: %d finite rows, m %.4f-%.4f'%(len(ms),ms[0],ms[-1]))
steps=np.round(np.diff(ms),6); u=np.unique(steps)
print(' unique steps n=%d:'%len(u), u[:15])
for mm in [5.082255,13.023095,33.371213,53.419632,62.49]:
    i=np.argmin(np.abs(ms-mm)); G=1.973269804593025e-16/ct[i]
    print('  m=%.4f ctau=%.6e m  Gamma=%.6e GeV'%(ms[i],ct[i],G))