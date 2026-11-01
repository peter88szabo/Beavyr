"""Transition properties for Beavyr. No environment activation is performed here.

State tracking uses Behemoth's biorthogonal occupied/virtual overlap construction
(src/hamiltonian/tddft/state_overlap.rs). TDDFT uses normalized X pseudo-CI
wavefunctions for tracking; near-degenerate/ambiguous matches are rejected.
"""
import os
import numpy as np


def spin_overlap(previous, current, cross):
    op, vp = previous
    oc, vc = current
    soo = op.T @ cross @ oc
    u, sigma, vt = np.linalg.svd(soo)
    if sigma.min() < 1.e-8:
        raise ValueError("Singular occupied-space overlap; cannot track excited state")
    v = vt.T
    sov = u.T @ (op.T @ cross @ vc)
    svo = (vp.T @ cross @ oc) @ v
    p = svo / sigma[None, :]
    q = sov / sigma[:, None]
    r = vp.T @ cross @ vc - p @ sov
    return u, v, sigma, p, q, r


def state_overlaps(reference, current, cross, root):
    blocks = [spin_overlap(a, b, cross) for a, b in zip(reference['orbitals'], current['orbitals'])]
    ref = reference['vectors'][root]
    refnorm = np.sqrt(sum(np.sum(x*x) for x in ref))
    result = []
    for state in current['vectors']:
        norm = np.sqrt(sum(np.sum(x*x) for x in state))
        if min(norm, refnorm) < 1.e-12:
            raise ValueError("Empty excitation vector")
        aa, bb, cc, det = [], [], [], []
        for xp, xc, block in zip(ref, state, blocks):
            u, v, sigma, p, q, r = block
            xp = u.T @ xp / refnorm
            xc = v.T @ xc / norm
            aa.append(np.sum(xp * p.T))
            bb.append(np.sum(xc * q))
            cc.append(np.sum((xp @ r) * xc / sigma[:, None]))
            det.append(np.prod(sigma) * np.linalg.det(u) * np.linalg.det(v))
        overlap = det[0]**2 * (2*aa[0]*bb[0]+cc[0]) if len(blocks)==1 else np.prod(det)*(sum(aa)*sum(bb)+sum(cc))
        result.append(overlap)
    return np.array(result)


def match(reference, current, cross, root):
    overlaps = state_overlaps(reference, current, cross, root)
    order = np.argsort(-np.abs(overlaps))
    best = order[0]
    if abs(overlaps[best]) < 0.7 or (len(order)>1 and abs(overlaps[best])-abs(overlaps[order[1]]) < 0.15):
        raise ValueError(f"Ambiguous excited-state match: overlaps={overlaps}. Reduce displacement or increase roots; derivatives were not accepted.")
    return current['dipoles'][best] * np.sign(overlaps[best]), int(best), float(overlaps[best])


def main():
    module = os.environ['BEAVYR_PROPERTY_BACKEND']
    root = int(os.environ['BEAVYR_ROOT'])-1
    roots = int(os.environ['BEAVYR_ROOTS'])
    charge, multi = int(os.environ['BEAVYR_CHARGE']), int(os.environ['BEAVYR_MULTIPLICITY'])
    basis, functional = os.environ['BEAVYR_BASIS'], os.environ['BEAVYR_FUNCTIONAL']
    tda = os.environ['BEAVYR_TDA']=='true'
    hf = os.environ['BEAVYR_HF']=='true'
    ht = os.environ['BEAVYR_HT']=='true'
    step = float(os.environ['BEAVYR_STEP'])
    cores = min(4, max(1, int(os.environ['BEAVYR_CORES'])))
    memory = max(256, int(os.environ['BEAVYR_MEMORY']))
    rows = open('geometry.bohr').read().splitlines()
    atoms = [line.split()[0] for line in rows]
    geometry = np.array([[float(x) for x in line.split()[1:]] for line in rows])
    if module=='pyscf':
        from pyscf import gto, scf, dft, lib
        lib.num_threads(cores)
        def evaluate(coords):
            mol=gto.M(atom=list(zip(atoms,coords)),unit='Bohr',basis=basis,charge=charge,spin=multi-1,verbose=3,max_memory=memory)
            mf=(scf.RHF(mol) if multi==1 else scf.UHF(mol)) if hf else (dft.RKS(mol) if multi==1 else dft.UKS(mol))
            if not hf: mf.xc=functional
            mf.conv_tol=1.e-11
            mf.max_cycle=150
            mf.kernel()
            if not mf.converged: raise ValueError('SCF did not converge')
            td=mf.TDA() if tda else mf.TDHF() if hf else mf.TDDFT()
            td.nstates=roots
            td.conv_tol=1.e-9
            td.max_cycle=150
            td.kernel()
            if len(td.e)<roots or not np.all(td.converged): raise ValueError('Excited-state calculation did not converge for every requested root')
            coefficients=[mf.mo_coeff] if multi==1 else mf.mo_coeff
            occupations=[mf.mo_occ] if multi==1 else mf.mo_occ
            orbitals=[(c[:,occ>0],c[:,occ==0]) for c,occ in zip(coefficients,occupations)]
            vectors=[[np.asarray(xy[0])] if multi==1 else [np.asarray(x) for x in xy[0]] for xy in td.xy]
            return dict(orbitals=orbitals,vectors=vectors,dipoles=np.asarray(td.transition_dipole()),mol=mol)
        def cross(a,b): return gto.intor_cross('int1e_ovlp',a['mol'],b['mol'])
    else:
        import psi4
        from psi4.driver.procrouting.response.scf_response import tdscf_excitations
        psi4.set_num_threads(cores)
        psi4.set_memory(f'{memory} MB')
        psi4.core.set_output_file('psi4-properties.log',False)
        def evaluate(coords):
            psi4.core.clean()
            text=f'{charge} {multi}\n'+'\n'.join(f'{a} {x:.14f} {y:.14f} {z:.14f}' for a,(x,y,z) in zip(atoms,coords))+'\nunits bohr\nsymmetry c1\nno_reorient\nno_com\n'
            mol=psi4.geometry(text)
            psi4.set_options({'basis':basis,'save_jk':True,'reference':('rhf' if multi==1 else 'uhf') if hf else ('rks' if multi==1 else 'uks'),'e_convergence':1.e-11,'d_convergence':1.e-10,'maxiter':150})
            _,wfn=psi4.energy('hf' if hf else functional,molecule=mol,return_wfn=True)
            states=tdscf_excitations(wfn,states=roots,tda=tda,r_convergence=1.e-8,maxiter=150)
            if len(states)<roots: raise ValueError('Insufficient converged excited roots')
            orbitals=[(np.array(wfn.Ca_subset('AO','OCC')),np.array(wfn.Ca_subset('AO','VIR')))]
            if multi!=1: orbitals.append((np.array(wfn.Cb_subset('AO','OCC')),np.array(wfn.Cb_subset('AO','VIR'))))
            spins=['ALPHA'] if multi==1 else ['ALPHA','BETA']
            # Psi4 R=X+Y, L=X-Y. Reconstruct X for Behemoth's pseudo-CI overlap.
            vectors=[[0.5*(np.array(s['RIGHT EIGENVECTOR '+spin])+np.array(s['LEFT EIGENVECTOR '+spin])) for spin in spins] for s in states]
            return dict(orbitals=orbitals,vectors=vectors,dipoles=np.array([s['ELECTRIC DIPOLE TRANSITION MOMENT (LEN)'] for s in states]),basis=wfn.basisset())
        def cross(a,b): return np.array(psi4.core.MintsHelper(a['basis']).ao_overlap(a['basis'],b['basis']))
    reference=evaluate(geometry)
    derivatives=[]
    with open('tracking.csv','w') as log:
        log.write('coordinate,displacement_bohr,matched_root,signed_overlap\n')
        if ht:
            for i in range(geometry.size):
                values=[]
                for sign in [1,-1]:
                    displaced=geometry.copy()
                    displaced.flat[i]+=sign*step
                    current=evaluate(displaced)
                    dipole,selected,overlap=match(reference,current,cross(reference,current),root)
                    log.write(f'{i+1},{sign*step},{selected+1},{overlap:.12f}\n');log.flush()
                    values.append(dipole)
                derivatives.append((values[0]-values[1])/(2*step))
                print(f'Dipole derivative {i+1}/{geometry.size}',flush=True)
    if not np.isfinite(reference['dipoles'][root]).all() or not np.isfinite(derivatives).all(): raise ValueError('Non-finite transition properties')
    with open('transition.properties','w') as f:
        f.write('BEAVYR_TRANSITION_V1\n'+str(len(atoms))+'\n')
        for a,xyz in zip(atoms,geometry):f.write(a+' '+' '.join(f'{v:.14e}' for v in xyz)+'\n')
        for vector in [reference['dipoles'][root],*derivatives]: f.write(' '.join(f'{v:.14e}' for v in vector)+'\n')
    print('BEAVYR_TRANSITION_SUCCESS',flush=True)

if __name__=='__main__':main()
