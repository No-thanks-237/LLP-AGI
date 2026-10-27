"""Read-only review of v3.4. Output is confined to this review directory."""
from pathlib import Path
import hashlib,json,math,random,sys,tarfile,time
import numpy as np
HERE=Path(__file__).resolve().parent
RUN=HERE.parent.parent
sys.path.insert(0,str(RUN/'work/analysis'))
import bc5_lib as L
import bc5_production as P
import run_scan as R
import run_sensitivity as RS
import scalar_widths as W
out={'model_tag':R.MODEL_TAG,'date':'2026-10-25','hashes':{}}
sources=['work/analysis/bc5_lib.py','work/analysis/bc5_production.py','work/analysis/run_scan.py',
         'work/analysis/run_sensitivity.py','work/widths/scalar_widths.py',
         'work/widths/reference_br_model.py','work/magnus/samples/normalisation.json',
         'outputs/data/signal_scan.json','outputs/data/signal_scan_refbr.json']
for p in sources:
    out['hashes'][p]=hashlib.sha256((RUN/p).read_bytes()).hexdigest()
def save():
    (HERE/'recheck_results.json').write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')

# Neutral kaon daughters must share a vertex. Existing code groups each kaon's two charged
# daughters consecutively, so within each adjacent pair the two non-null offsets must agree.
rng=random.Random(25101)
same=[]; examples=[]
for _ in range(4000):
    p4=(math.sqrt(1.5**2+8**2+100**2),8,0,100)
    d=L.s_decay_products(1.5,'KK',rng,s_p4=p4)
    if not d or d[0][6] is None:continue
    for i in range(0,len(d),2):
        if i+1>=len(d):continue
        a=np.array(d[i][6]); b=np.array(d[i+1][6])
        same.append(np.linalg.norm(a-b)<1e-12)
        if len(examples)<3:examples.append({'first_offset_m':a.tolist(),'second_offset_m':b.tolist()})
out['neutral_K_shared_vertex']={'pairs':len(same),'shared_vertex_fraction':float(np.mean(same)),
    'examples':examples}
print('neutral K',out['neutral_K_shared_vertex'],flush=True)

# Four-body kinematics are checked both as a distribution and for rejected generation attempts.
rng=random.Random(25102)
n=20000
events=[L.nbody(1.0,[0.0]*4,rng) for _ in range(n)]
a=np.array([x for x in events if x is not None])
out['massless_four_body']={'attempts':n,'none':sum(x is None for x in events),
    'mean_E':a[:,:,1].mean(axis=0).tolist(),
    'se_E':(a[:,:,1].std(axis=0)/math.sqrt(len(a))).tolist(),
    'max_E_residual':float(abs(a[:,:,1].sum(axis=1)-1).max()),
    'max_p_residual':float(abs(a[:,:,2:5].sum(axis=1)).max())}
out['massive_four_body_failures']={}
for m in [0.5584,0.56,0.565616,0.58,0.6,0.661654,0.8,1.0,1.5]:
    t=time.time();n=1000
    ds=[L.nbody(m,[W.M_PI_CHARGED]*4,rng) for _ in range(n)]
    good=np.array([x for x in ds if x is not None])
    out['massive_four_body_failures'][str(m)]={'attempts':n,'none':sum(x is None for x in ds),
       'failure_fraction':sum(x is None for x in ds)/n,
       'mean_E_good':good[:,:,1].mean(axis=0).tolist() if len(good) else [],
       'multi_BR':L.channel_branching(m)[0].get('multi',0),'seconds':time.time()-t}
print('phase space',out['massless_four_body'],out['massive_four_body_failures'],flush=True)

# Independent quadrature of the massive three-body marginal.
n=15000
ds=[L.nbody(1.0,[W.M_PI_CHARGED]*3,rng) for _ in range(n)]
a=np.array([x for x in ds if x is not None])
s23=1+W.M_PI_CHARGED**2-2*a[:,0,1]
s=np.linspace(4*W.M_PI_CHARGED**2,(1-W.M_PI_CHARGED)**2,100001)
lam1=(1-(W.M_PI_CHARGED+np.sqrt(s))**2)*(1-(W.M_PI_CHARGED-np.sqrt(s))**2)
lam2=s*(s-4*W.M_PI_CHARGED**2)
w=np.sqrt(np.maximum(lam1*lam2,0))/s
out['massive_three_body']={'mean_s23':float(s23.mean()),'se':float(s23.std()/math.sqrt(len(s23))),
    'independent_integral':float(np.trapezoid(s*w,s)/np.trapezoid(w,s)),
    'none':sum(x is None for x in ds)}

# Integral of the linear response over the exponential, compared to high-precision mpmath.
import mpmath as mp
mp.mp.dps=60
out['lifetime']={'constant_response':[],'linear_response':[]}
for z0 in [0,1.01,1.249,2.4,2.6]:
    for alpha in [.1,1,10,100]:
        got=L._p_window([alpha],[z0],np.ones((1,13)),L.ZGRID)[0]
        exact=0 if z0>=2.5 else math.exp(-alpha*(max(1,z0)-z0))-math.exp(-alpha*(2.5-z0))
        out['lifetime']['constant_response'].append({'z0':z0,'alpha':alpha,'got':float(got),'exact':exact})
for alpha in np.logspace(-12,3,80):
    aa=mp.mpf(str(alpha)); lo=mp.mpf(1); hi=mp.mpf('2.5')
    xx=aa*(hi-lo)
    exact=mp.exp(-aa*lo)*(1-(1+xx)*mp.exp(-xx))/xx
    got=L._p_window([alpha],[0],np.linspace(0,1,13)[None,:],L.ZGRID)[0]
    if exact>mp.mpf('1e-290'):
        out['lifetime']['linear_response'].append({'alpha':float(alpha),'relative_error':float(mp.mpf(got)/exact-1)})
print('massive 3body',out['massive_three_body'],'lifetime max linear relative error',
      max(abs(x['relative_error']) for x in out['lifetime']['linear_response']),flush=True)

scan=json.loads((RUN/'outputs/data/signal_scan.json').read_text())
out['scan']={k:v for k,v in scan.items() if k not in ['signal_per_fb_per_channel','theta2','nsamp_by_mass','masses']}
out['scan']['mass_count']=len(scan['masses'])
out['MC']={}
for lum in [25,300]:
    cr=RS.crossing_errors(scan,2.3,lum)
    los=[v[2] for v in cr.values() if np.isfinite(v[2])]
    his=[v[3] for v in cr.values() if np.isfinite(v[3])]
    out['MC'][str(lum)]={'lo_median':float(np.median(los)),'lo_max':float(max(los)),
       'lo_over5':sum(v>.05 for v in los),'lo_n':len(los),
       'hi_median':float(np.median(his)),'hi_max':float(max(his)),
       'hi_over5':sum(v>.05 for v in his),'hi_n':len(his)}
print('MC',out['MC'],flush=True)
save()

# On-disk Higgs input, not a claimed change in a report.
d=np.loadtxt(RUN/'work/magnus/samples/higgs.txt.gz',comments='#')
mask=abs(d[:,10])==62
a=d[mask]
eta=np.arcsinh(a[:,4]/np.maximum(np.hypot(a[:,2],a[:,3]),1e-30))
ev,counts=np.unique(a[:,0],return_counts=True)
out['higgs']={'rows':len(d),'selected_rows':len(a),'selected_unique_events':len(ev),
    'max_selected_per_event':int(counts.max()),'max_abs_eta':float(abs(eta).max()),
    'all_selected_within_eta6':bool(np.all(abs(eta)<=6.00001)),
    'mean_pT':float(np.hypot(a[:,2],a[:,3]).mean()),
    'all_record_unique_events':len(np.unique(d[:,0]))}
archive=RUN/'work/magnus/bundles/higgs3.tar.gz'
with tarfile.open(archive) as t:
    src=t.extractfile('parents.cc').read().decode()
    out['higgs']['archived_eta_filter_lines']=[(i+1,l.strip()) for i,l in enumerate(src.splitlines()) if 'eta()' in l]
print('Higgs',out['higgs'],flush=True)
save()

# Library support at every actually open channel of the production grid.
L.load_hadron_lib(str(RUN/'work/magnus/hadron_lib'))
out['library']={}
for ch,c in L.HADRON_LIB.chan.items():
    mapch='tau' if ch=='tautau' else ch
    checks=[]
    for m in scan['masses']:
        br=L.channel_branching(m)[0].get(mapch,0)
        if br<=0:continue
        node=float(c['mass_nodes'][np.argmin(abs(c['mass_nodes']-m))])
        checks.append({'requested':m,'node':node,'rel_mismatch':abs(node-m)/max(node,m)})
    # verify the node's stable charged subset has a physically possible neutral complement
    starts=c['start'];energy=np.add.reduceat(c['e'],starts)
    px=np.add.reduceat(c['px'],starts);py=np.add.reduceat(c['py'],starts);pz=np.add.reduceat(c['pz'],starts)
    nodes=c['mass'][starts]
    impossible=energy+np.sqrt(px*px+py*py+pz*pz)>nodes+1e-4
    # rotational invariance through the public decay sampler (which rotates the raw library)
    r=random.Random(25103);testmass=62.49
    cos2=[]
    for _ in range(3000):
        for v in L.s_decay_products(testmass,mapch,r):
            pp=sum(x*x for x in v[2:5])
            if pp>0:cos2.append(v[4]**2/pp)
    out['library'][ch]={'nodes':len(c['mass_nodes']),'attempts':int(c['n_attempt']),
        'open_grid_nodes':checks,'max_mismatch':max(x['rel_mismatch'] for x in checks),
        'impossible_charged_subset_fraction_at_own_node':float(np.mean(impossible)),
        'rotated_cos2_mean_62p49':float(np.mean(cos2))}
print('libraries', {k:{p:v[p] for p in ['nodes','max_mismatch','impossible_charged_subset_fraction_at_own_node','rotated_cos2_mean_62p49']} for k,v in out['library'].items()},flush=True)
save()

# Conditional quartic density sampled vs a direct integration of the paper's integrand.
r=random.Random(25104)
out['quartic_q2']={}
for channel,spec in [('K',P.M_KP),('K*892',.892),('K0*700',.700),('K1_1270',1.270)]:
    ms=1.0;mx=P.M_BP;n=10000
    draws=[P.sample_dipion_mass(ms,mx,spec,channel,r) for _ in range(n)]
    ss=np.array([v*v for v in draws if v is not None])
    q=np.linspace(4*ms**2,(mx-spec)**2,20001)
    den=np.array([P._q2_weight(ms,mx,spec,channel,v) for v in q])
    expected=np.trapezoid(q*den,q)/np.trapezoid(den,q)
    out['quartic_q2'][channel]={'n':n,'none':sum(v is None for v in draws),'mean':float(ss.mean()),
        'se':float(ss.std()/math.sqrt(len(ss))),'quadrature':float(expected)}

out['widths']={str(m):W.widths(m) for m in [.3,1,1.5,2,2.000001,3,5,10,30,62.49]}
out['B_formula_checks']={}
for m in [0,1,3,4.699999,4.7,4.78]:
    g2=(P.M_B_QUARK/P.V)**2*P.XI_SB**2
    inc=g2*(P.M_BP**2-m*m)**2/(32*math.pi*P.M_BP**3*P.GAM_BP)
    fk=.33/(1-m*m/37.5)
    lam=(1-(P.M_KP-m)**2/P.M_BP**2)*(1-(P.M_KP+m)**2/P.M_BP**2)
    me2=.25*(P.M_BP**2-P.M_KP**2)**2/(P.M_B_QUARK-P.M_S_QUARK)**2*fk*fk
    exc=g2*me2*math.sqrt(max(lam,0))/(16*math.pi*P.M_BP*P.GAM_BP)
    got=sum(P.br_B_to_S(m).values())
    out['B_formula_checks'][str(m)]={'actual':got,'inclusive':inc,'exclusive':exc,
        'actual_over_expected':got/(inc if m<4.7 else exc)}
print('quartic',out['quartic_q2'],'B formulas',out['B_formula_checks'],flush=True)
save()

assert all(hashlib.sha256((RUN/p).read_bytes()).hexdigest()==h for p,h in out['hashes'].items())
print('DONE; source and production artifacts unchanged.',flush=True)
