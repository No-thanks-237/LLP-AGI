"""Audit charge-dependent endpoints, multicomponent contours, and vertex quadrature."""
from pathlib import Path
import json,math,sys
import numpy as np
from scipy.integrate import quad
HERE=Path(__file__).resolve().parent
RUN=HERE.parents[2]
sys.path.insert(0,str(RUN/'work/analysis'))
import models as M
import event_engine as E
import response_table as RT
import final_analysis as FA
import twohdm_lib as L


def main():
    L.load_hadron_lib(str(RUN/'work/magnus/hadron_lib'))
    model=M.TwoHDMModel(verbose=False)
    P=model.P
    out={'endpoints_GeV':{'charged_B_to_charged_K':M.M_BP-M.M_KP,
                         'neutral_B_to_neutral_K':M.M_B0-M.M_K0L,
                         'neutral_B_code_endpoint_using_charged_K':M.M_B0-M.M_KP}}
    mh=4.784
    tab=RT.load_table(RT.cache_path(model,mh))
    src=json.load(open(RUN/'work/analysis/scan_2hdm/m_4.784.json'))
    nominal_weight=model.history_weights
    def corrected_weight(m,t,ta):
        w=nominal_weight(m,t,ta)
        if m>=M.M_B0-M.M_K0L and 'B:mix:B0' in w:w['B:mix:B0']=0.
        return w
    rows=[]
    for lumi in (25.,300.):
        c=src['curves'][str(lumi)]['crossings'][-1]
        tb=c['param_refined']
        old=E.evaluate(tab,model,mh,tb,lumi)
        model.history_weights=corrected_weight
        new=E.evaluate(tab,model,mh,tb,lumi)
        cr=E.refine_crossing(tab,model,mh,lumi,tb*.5,tb*1.1)
        model.history_weights=nominal_weight
        rows.append({'lumi_fb':lumi,'tanb_original':tb,'tanb_removing_closed_B0':cr['param_refined'],
                     'signal_original':old['s_total'],'signal_removing_closed_B0':new['s_total'],
                     'relative_boundary_change':cr['param_refined']/tb-1})
    out['m4p784_closed_B0_experiment']=rows
    out['inclusive_rate_near_charged_endpoint']=[{'mass':m,'Bplus':P.br_B_to_XsH(m,100,'B+')['total'],
                                                'Bzero':P.br_B_to_XsH(m,100,'B0')['total']}
               for m in (4.784,4.78566-1e-8,4.785663-1e-9,4.785663+1e-9,4.785973+1e-9)]
    multi={}
    for mode in ('paper_mb','exclusive'):
        s=FA.load_scan(RUN/('work/analysis/scan_2hdm_'+mode)/'signal_scan.json')
        rr=FA.boundary_rows(s)
        cases=[]
        for mass,meta in s['mass_files'].items():
            for lumi,cs in meta['crossings'].items():
                if len(cs)>2:
                    row=next(r for r in rr if r['mass_GeV']==float(mass) and r['lumi_fb']==float(lumi))
                    cases.append({'m':float(mass),'lumi':float(lumi),'crossings':cs,'boundary_csv_row':row})
        multi[mode]=cases
    out['multiple_crossings']=multi
    # Constant acceptance has an independent exact convolution CDF.
    vals=[]
    for alpha in (.001,.1,1,10,30):
        lam=1/alpha
        calc=float(L.p_window_vertex(1.,np.array([1.]),np.ones((1,13)),1/alpha,np.array([lam]))[0])
        exact=(1+alpha)*math.exp(-alpha)-(1+2.5*alpha)*math.exp(-2.5*alpha)
        vals.append({'alpha':alpha,'lambda':lam,'relative_difference':calc/exact-1,'exact':exact})
    out['vertex_equal_scales_removable_singularity']=vals
    # Parent-subsampling + an overall yield factor is only linear-order unbiased for
    # event OR: exactly two parents, each q, retaining exactly one of the two.
    q=.1
    out['parent_subsampling_OR_toy']={'q':q,'true':1-(1-q)**2,'estimator':2*q,
                                     'relative_bias':2*q/(1-(1-q)**2)-1}
    nominal=FA.load_scan(RUN/'outputs/data/signal_scan.json')
    points=[p for p in nominal['points'] if p['mH']<4.786 and abs(p['tanb']-5)<1e-10]
    rel=[p['ds']['25']/p['s']['25'] for p in points]
    out['nominal_tanb5_mc']={'min':min(rel),'median':float(np.median(rel)),'max':max(rel),
                            'worst_mass':points[int(np.argmax(rel))]['mH']}
    (HERE/'physics_edges_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
