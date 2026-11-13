"""Separate kinetic moments, raw HDG fields and mixed source errors.

Compare the saved GSIS step 150 with CIS step 35. The GSIS stopping test
was not met: these are snapshot differences, not converged GSIS errors.
Compute both cell-average L2 and exact-in-space P4 L2 from saved nodal data.
"""
from pathlib import Path
import csv
import json

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'results/reference_cavity_residual_platform'


def reference_p4_integration():
    # Same equispaced reference nodes/order as Basis_Function.f90.
    nodes=np.array([(i/4,j/4) for j in range(5) for i in range(5-j)])
    powers=[(degree-j,j) for degree in range(5) for j in range(degree+1)]
    vandermonde=np.array([[x**i*y**j for i,j in powers] for x,y in nodes])
    coefficients=np.linalg.inv(vandermonde)
    x,w=np.polynomial.legendre.leggauss(6);x=(x+1)/2;w=w/2
    points=np.array([(u,(1-u)*t) for u in x for t in x])
    weights=np.array([wu*wt*(1-u) for u,wu in zip(x,w) for wt in w])
    monomials=np.array([[xx**i*yy**j for i,j in powers] for xx,yy in points])
    shape=np.einsum('pi,ij->pj',monomials,coefficients)
    assert np.max(abs(shape.sum(axis=1)-1))<1e-12
    mass=np.einsum('pi,p,pj->ij',shape,weights,shape)
    integral=np.einsum('pi,p->i',shape,weights)
    assert abs(mass.sum()-.5)<1e-12 and np.linalg.eigvalsh(mass).min()>0
    return mass,integral


def main():
    cells={mode:np.loadtxt(FOLDER/mode/'final_cells.dat') for mode in ['cis','gsis']}
    nodes={mode:np.loadtxt(FOLDER/mode/'final_nodes.dat').reshape(-1,15,22) for mode in ['cis','gsis']}
    area=cells['cis'][:,1];mass,integral=reference_p4_integration()
    assert np.allclose(cells['cis'][:,:4],cells['gsis'][:,:4])
    for mode in ['cis','gsis']:
        assert np.max(abs(2*area[:,None]*integral-nodes[mode][:,:,2]))<1e-13
        for node_slice,cell_slice in [(slice(3,10),slice(4,11)),(slice(10,16),slice(11,17)),(slice(16,22),slice(17,23))]:
            recovered=2*np.einsum('i,tik->tk',integral,nodes[mode][:,:,node_slice])
            assert np.max(abs(recovered-cells[mode][:,cell_slice]))<1e-11
    norms={
        'cell_average_L2':lambda z:np.sqrt(np.einsum('t,tk,tk->',area,z,z)),
        'full_P4_L2':lambda z:np.sqrt(np.einsum('t,tik,ij,tjk->',2*area,z,mass,z)),
    }
    fields={
        'cell_average_L2':dict(cis=cells['cis'][:,4:10],kinetic=cells['gsis'][:,4:10],
                               mixed_source=cells['gsis'][:,11:17],raw_hdg=cells['gsis'][:,17:23]),
        'full_P4_L2':dict(cis=nodes['cis'][:,:,3:9],kinetic=nodes['gsis'][:,:,3:9],
                          mixed_source=nodes['gsis'][:,:,10:16],raw_hdg=nodes['gsis'][:,:,16:22]),
    }
    rows=[]
    for discretization,group in fields.items():
        norm=norms[discretization]
        for baseline in ['cis','kinetic']:
            reference=group[baseline]
            for name in ['kinetic','mixed_source','raw_hdg']:
                value=group[name]
                for quantity,ids,offset in [('velocity',[1,2],0),('density_perturbation',[0],0),
                                          ('complete_density',[0],1),('temperature_perturbation',[3],0),
                                          ('complete_temperature',[3],1),('heat_flux',[4,5],0)]:
                    difference=value[...,ids]-reference[...,ids]
                    absolute=norm(difference);denominator=norm(reference[...,ids]+offset)
                    rows.append(dict(discretization=discretization,baseline=baseline,quantity=quantity,
                                     gsis_field=name,absolute_L2=float(absolute),reference_L2=float(denominator),
                                     relative_percent=float(100*absolute/denominator)))
    g=fields['cell_average_L2'];correction=g['raw_hdg']-g['kinetic'];feedback=g['mixed_source']-g['kinetic']
    mixing=np.sum(area[:,None]*correction*feedback)/np.sum(area[:,None]*correction**2)
    gn=fields['full_P4_L2']
    assert np.max(abs(gn['mixed_source']-gn['kinetic']-mixing*(gn['raw_hdg']-gn['kinetic'])))<1e-11
    metadata=dict(cis_step=35,gsis_step=150,gsis_converged=False,reference='Converged CIS at the same mesh, model and velocity grid, not an exact physical reference.',
                  density_and_temperature='Both perturbation-normalized and complete-field-normalized errors are provided; absolute differences are identical.',
                  spatial_integration='Full P4 errors use 6x6 Gauss-Legendre Duffy integration, exact for the squared degree-4 interpolant. Cell means are checked against saved reference outputs.',
                  mixing_coefficient=float(mixing),errors=rows)
    with (FOLDER/'macro_field_errors.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (FOLDER/'macro_field_errors.json').write_text(json.dumps(metadata,indent=2)+'\n')
    for row in rows:
        if row['discretization']=='full_P4_L2' and row['baseline']=='cis':
            print(row['gsis_field'],row['quantity'],f"{row['relative_percent']:.9g}%",flush=True)
    print('PASS P4 integration, recovered cell means and nodal source mixing identity. Mixing coefficient:',mixing)


if __name__=='__main__':main()
