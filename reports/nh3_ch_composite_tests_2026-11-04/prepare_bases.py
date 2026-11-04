#!/usr/bin/env python3
"""Convert installed Psi4 jun bases to job-local CFOUR/MRCC format."""
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import shutil
from pyscf.gto.basis import parse_gaussian

ROOT=Path(__file__).resolve().parent/'bases'
ROOT.mkdir(exist_ok=True)
source=Path('/home/peter/micromamba/envs/psi4env/share/psi4/basis')
metadata=[]
custom=[]
for zeta in ('t','q'):
    label=f'jun-cc-pV{zeta.upper()}Z'
    path=source/f'jun-cc-pv{zeta}z.gbs'
    shutil.copy2(path,ROOT/path.name)
    counts={}
    for element in ('H','C','N'):
        shells=parse_gaussian.load(str(path),element)
        groups=defaultdict(list)
        for shell in shells:
            l=shell[0]
            assert isinstance(l,int) and all(len(p)==len(shell[1]) for p in shell[1:])
            for col in range(1,len(shell[1])):
                groups[l].append({p[0]:p[col] for p in shell[1:]})
        angular=sorted(groups)
        exponents={l:sorted({e for c in groups[l] for e in c},reverse=True) for l in angular}
        block=[f'{element}:{label}',f'Converted without coefficient changes from Psi4 {path.name}', '',
               str(len(angular)), ' '.join(map(str,angular)),
               ' '.join(str(len(groups[l])) for l in angular),
               ' '.join(str(len(exponents[l])) for l in angular), '']
        for l in angular:
            block.append(' '.join(f'{e:.16E}' for e in exponents[l]))
            block.append('')
            for e in exponents[l]:
                block.append(' '.join(f'{c.get(e,0):.16E}' for c in groups[l]))
            block.append('')
        # Reparse the numerical CFOUR payload and check every coefficient.
        tokens=' '.join(block[3:]).split(); k=0
        n=int(tokens[k]);k+=1
        ls=list(map(int,tokens[k:k+n]));k+=n
        nc=list(map(int,tokens[k:k+n]));k+=n
        ne=list(map(int,tokens[k:k+n]));k+=n
        assert ls==angular
        for l,cnum,enum in zip(ls,nc,ne):
            es=list(map(float,tokens[k:k+enum]));k+=enum
            assert es==exponents[l]
            for e in es:
                values=list(map(float,tokens[k:k+cnum]));k+=cnum
                assert values==[c.get(e,0) for c in groups[l]]
        assert k==len(tokens)
        counts[element]=sum((2*l+1)*len(groups[l]) for l in angular)
        custom.append('\n'.join(block)+'\n')
    expected={'H':14,'C':39,'N':39} if zeta=='t' else {'H':30,'C':71,'N':71}
    assert counts==expected,(label,counts)
    metadata.append(dict(label=label,source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),spherical_basis_functions=counts))
(ROOT/'jun_GENBAS').write_text('\n'.join(custom))
standard='\n'.join((Path('/home/peter/Programs/MRCC/binary/BASIS')/e).read_text() for e in ('H','C','N'))
(ROOT/'GENBAS').write_text(standard+'\n'+'\n'.join(custom))
(ROOT/'sources.json').write_text(json.dumps(metadata,indent=2)+'\n')
print(json.dumps(metadata,indent=2))
