"""Reproducible hand-written synthetic HepMC reader/geometry fixtures.

No event generator, physics model, or tested C++ transport is imported.
The endpoint oracle uses independent z-variable RK4 and straight-line geometry.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

SEED = 314159  # identifier; fixtures contain no random draws
EPS = 1e-5     # m, safely above ROOT Float_t rounding at these boundaries
KB = 0.299792458
MUON_MASS = 0.1056583755


def particle(x=.5, y=.2, z=2., px=0., py=0., pz=20., charge=1, neutral=False):
    mass = 0. if neutral else MUON_MASS
    return dict(pid=22 if neutral else (-13 if charge > 0 else 13), status=1,
                charge=0 if neutral else charge, mass=mass,
                vertex=[x, y, z, .17], momentum=[px, py, pz, math.sqrt(px*px+py*py+pz*pz+mass*mass)])


def oracle(p):
    x,y,z,ct = p['vertex']; px,py,pz,e = p['momentum']
    q = p['charge']; norm = math.sqrt(px*px+py*py+pz*pz); path = 0.
    if q == 0 or pz <= 0 or z > 7.7:
        return None
    def straight(zend):
        nonlocal x,y,z,ct,path
        dz=zend-z; ds=dz*norm/pz
        x+=dz*px/pz; y+=dz*py/pz; ct+=ds*e/norm; path+=ds; z=zend
    if z < 3.: straight(3.)
    if z < 7.:
        # This invariant decides whether the forward branch reaches the exit.
        target_px=px-KB*q*(7.-z)
        if px*px+pz*pz-target_px*target_px <= 0: return None
        n=2000; h=(7.-z)/n
        state=[x,y,px,pz,path,ct]
        def deriv(v):
            xx,yy,ppx,ppz,ll,tt=v
            return [ppx/ppz,py/ppz,-KB*q,KB*q*ppx/ppz,norm/ppz,e/ppz]
        def add(a,b,f): return [aa+f*bb for aa,bb in zip(a,b)]
        for _ in range(n):
            a=deriv(state); b=deriv(add(state,a,h/2)); c=deriv(add(state,b,h/2)); d=deriv(add(state,c,h))
            state=[v+h*(aa+2*bb+2*cc+dd)/6 for v,aa,bb,cc,dd in zip(state,a,b,c,d)]
        x,y,px,pz,path,ct=state; z=7.
    straight(7.7); front=[x,y,z,ct]
    straight(9.4); end=[x,y,z,ct]
    if any(abs(v[0])>3.24 or abs(v[1])>2.415 for v in (front,end)): return None
    dx=end[0]-front[0];dy=end[1]-front[1];den=dx*dx+dy*dy
    t=max(0.,min(1.,-(front[0]*dx+front[1]*dy)/den)) if den else 0.
    if math.hypot(front[0]+t*dx,front[1]+t*dy)<.09: return None
    return dict(front=front,end=end,momentum=[px,py,pz,e],path=path)


def cases():
    out=[]
    def add(name,p,accepted): out.append((name,[[p]],[accepted]))
    for boundary in (3.,7.,7.7,9.4):
        for label,d in [('minus',-EPS),('equal',0.),('plus',EPS)]:
            add(f'birth_z_{boundary}_{label}',particle(z=boundary+d),boundary+d<=7.7)
    for axis,limit in [('x',3.24),('y',2.415)]:
        for sign in (-1,1):
            for label,d in [('inside',-EPS),('equal',0.),('outside',EPS)]:
                kw=dict(z=7.7,x=.5,y=.2);kw[axis]=sign*(limit+d)
                add(f'rect_{axis}_{sign}_{label}',particle(**kw),d<=0)
    for label,d in [('inside',-EPS),('equal',0.),('outside',EPS)]:
        add(f'hole_radius_{label}',particle(z=7.7,x=.09+d,y=0),d>=0)
    add('hole_cross_between_good_endpoints',particle(z=7.7,x=.17,y=0,px=-2,pz=10),False)
    # Front inside, back just outside/inside: avoids an ill-conditioned exact slope.
    for label,d in [('inside',-EPS),('outside',EPS)]:
        add(f'back_rect_{label}',particle(z=7.7,x=.5,y=.2,px=(3.24+d-.5)/1.7*10,pz=10),d<0)
    add('positive_charge',particle(charge=1),True)
    add('negative_charge',particle(charge=-1),True)
    add('turn_before_field_exit',particle(z=3,pz=.1),False)
    add('backward',particle(pz=-20),False)
    add('neutral_only_zero_tracks',particle(neutral=True),False)
    add('LLP_birth_not_decay',particle(x=.6,y=.3,z=2.3,px=.3,py=.1),True)
    out.append(('different_parents_same_event',[[particle(x=.5,y=.2,z=2.)],
                                             [particle(x=.7,y=-.2,z=2.4,charge=-1)]],[True,True]))
    out.append(('opposite_pair_same_parent',[[particle(x=.5,y=.2),particle(x=.5,y=.2,charge=-1)]],[True,True]))
    return out


def construct():
    result=[]
    for number,(name,groups,expected_acceptance) in enumerate(cases(),1):
        parents=[]; daughters=[]; cursor=1+len(groups); accept_cursor=0
        for parent_index,group in enumerate(groups):
            mom=[sum(p['momentum'][k] for p in group) for k in range(4)]
            parent=dict(pid=25,status=2,charge=0,mass=math.sqrt(max(0,mom[3]**2-sum(v*v for v in mom[:3]))),
                        vertex=[0.,0.,0.,0.],momentum=mom,m1=0,m2=-1,d1=cursor,d2=cursor+len(group)-1,
                        end_vertex=-(parent_index+2),birth_vertex=-1,accepted=False)
            parents.append(parent)
            for p in group:
                p=dict(p,m1=parent_index+1,m2=-1,d1=-1,d2=-1,end_vertex=0,birth_vertex=-(parent_index+2))
                pred=oracle(p); expected=expected_acceptance[accept_cursor];accept_cursor+=1
                if bool(pred)!=expected: raise AssertionError(f'oracle disagrees with geometry design: {name}')
                p['accepted']=expected;p['transport']=pred;daughters.append(p);cursor+=1
        total=[sum(p['momentum'][k] for p in parents) for k in range(4)]
        # An explicit incoming source particle supplies HepMC3's implicit root
        # vertex 0, which the bundled reader expects in addition to E's vertices.
        source=dict(pid=25,status=4,charge=0,mass=math.sqrt(max(0,total[3]**2-sum(v*v for v in total[:3]))),
                    vertex=[0.,0.,0.,0.],momentum=total,m1=-1,m2=-1,d1=1,d2=len(parents),
                    end_vertex=-1,birth_vertex=0,accepted=False)
        weight=(-.5 if number%7==0 else 1+number/100)
        result.append(dict(number=number,name=name,weight=weight,weights=[weight,2*weight],
                           groups=groups,particles=[source]+parents+daughters,parent_count=len(groups)))
    return result


def fmt(v): return format(v,'.17g')


def hepmc(events,version,momentum_unit,position_unit):
    pf=1000. if momentum_unit=='MEV' else 1.;vf=100. if position_unit=='CM' else 1000.
    lines=['HepMC::Version '+('2.06.11' if version==2 else '3.02.05'),
           'HepMC::'+('IO_GenEvent' if version==2 else 'Asciiv3')+'-START_EVENT_LISTING']
    for ev in events:
        nvert=1+ev['parent_count'];parts=ev['particles']
        if version==2:
            lines.append(f"E {ev['number']} 0 0 0 0 0 -1 {nvert} 0 0 0 2 "+' '.join(map(fmt,ev['weights'])))
        else: lines.append(f"E {ev['number']} {nvert} {len(parts)}")
        lines.append(f'U {momentum_unit} {position_unit}')
        if version==2: lines += ['C 0 0','F 0 0 0 0 0 0 0']
        else:
            lines += ['W '+' '.join(map(fmt,ev['weights'])), 'A 0 mpi 0','A 0 signal_process_id 0',
                      'A 0 event_scale 0','A 0 alphaQCD 0','A 0 alphaQED 0',
                      'A 0 GenCrossSection 0 0','A 0 GenPdfInfo 0 0 0 0 0 0 0']
            p=parts[0];mom=' '.join(fmt(v*pf) for v in p['momentum']);mass=fmt(p['mass']*pf)
            lines.append(f"P 1 0 {p['pid']} {mom} {mass} {p['status']}")
        for vertex in range(-1,-nvert-1,-1):
            indices=[i for i,p in enumerate(parts) if p['birth_vertex']==vertex]
            pos=parts[indices[0]]['vertex']
            coords=' '.join(fmt(v*vf) for v in pos)
            if version==2:
                lines.append(f'V {vertex} 0 {coords} {int(vertex==-1)} {len(indices)} 0')
                if vertex==-1:
                    p=parts[0];mom=' '.join(fmt(v*pf) for v in p['momentum']);mass=fmt(p['mass']*pf)
                    lines.append(f"P 1 {p['pid']} {mom} {mass} {p['status']} 0 0 -1 0")
            else:
                incoming='1' if vertex==-1 else str(-vertex)
                lines.append(f'V {vertex} 0 [{incoming}] @ {coords}')
            for i in indices:
                p=parts[i];mom=' '.join(fmt(v*pf) for v in p['momentum']);mass=fmt(p['mass']*pf)
                if version==2: lines.append(f"P {i+1} {p['pid']} {mom} {mass} {p['status']} 0 0 {p['end_vertex']} 0")
                else: lines.append(f"P {i+1} {vertex} {p['pid']} {mom} {mass} {p['status']}")
    lines.append('HepMC::'+('IO_GenEvent' if version==2 else 'Asciiv3')+'-END_EVENT_LISTING')
    return '\n'.join(lines)+'\n'


def write(output):
    output.mkdir(parents=True,exist_ok=True)
    if any(output.iterdir()): raise FileExistsError('fixture output directory must be empty')
    events=construct();rows=[]
    for ev in events:
        rows.append(f"E {ev['number']} {ev['name']} {fmt(ev['weight'])} {len(ev['particles'])} "+str(sum(p['accepted'] for p in ev['particles'])))
        for i,p in enumerate(ev['particles']):
            t=p.get('transport') or dict(front=[0]*4,end=[0]*4,momentum=[0]*4,path=0)
            nums=[i,p['pid'],p['status'],p['charge'],p['m1'],p['m2'],p['d1'],p['d2'],
                  *p['vertex'],*p['momentum'],int(p['accepted']),*t['front'],*t['end'],*t['momentum'],t['path']]
            rows.append('P '+' '.join(map(fmt,nums)))
    (output/'expected.tsv').write_text('\n'.join(rows)+'\n',encoding='utf-8')
    files=[]
    for version in (2,3):
        for mu in ('GEV','MEV'):
            for lu in ('MM','CM'):
                name=f'boundary_hepmc{version}_{mu.lower()}_{lu.lower()}.hepmc'
                (output/name).write_text(hepmc(events,version,mu,lu),encoding='ascii');files.append(name)
    spec=dict(kind='synthetic_detector_software_fixture',seed=SEED,random_draws=0,epsilon_m=EPS,
              event_count=len(events),oracle='independent z-variable RK4 and designed straight-line boundary cases',
              generator_models=False,sensitivity=False,events=events,files=files,
              hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir()})
    (output/'expected.json').write_text(json.dumps(spec,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(event_count=len(events),file_count=len(files),seed=SEED,
                         accepted_tracks=sum(p['accepted'] for ev in events for p in ev['particles']),
                         zero_track_events=sum(not any(p['accepted'] for p in ev['particles']) for ev in events))))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True,type=Path)
    write(parser.parse_args().output)
