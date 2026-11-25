"""Null simulation of duplicated instruments (reinfection exposure, 2,188-case EC outcome).
Usage: python3 sim_duplicates.py Task1_harmonised_small.csv
Input: harmonised exposure-outcome data (columns SNP, beta.exposure, se.outcome).
Under the null, each SNP's outcome effect ~ N(0, se.outcome^2); IVW with multiplicative
random effects (SE floored at the fixed-effect SE); two-sided normal test at 0.05.
A duplicated row re-uses its lead SNP's outcome draw (upper bound for the effect of duplication).
Scenarios (400,000 replicates each, seed 2026):
  A  8 independent loci
  B  13 SNPs, duplicate copies removed
  C  16 SNPs = 13 + 3 duplicate copies (as in the unpruned 16-SNP analysis)
  D  13 SNPs + 3 extra copies of each of the three highest-weight loci
     (weight = beta.exposure^2 / se.outcome^2, among the 8 loci); 22 rows
"""
import csv, sys, numpy as np
from scipy import stats
rows=list(csv.DictReader(open(sys.argv[1])))
bx={r['SNP']:float(r['beta.exposure']) for r in rows}
sy={r['SNP']:float(r['se.outcome']) for r in rows}
base=list(bx)
copies16={'rs7829447':'rs11996757','rs7460431':'rs7463208','rs7464228':'rs7463208'}
loci8=["rs2486961","rs10208306","rs28530774","rs11996757","rs6999003","rs821932","rs7238797","rs28505079"]
top3=sorted(loci8,key=lambda s:-(bx[s]/sy[s])**2)[:3]
design=lambda L: np.array([base.index(s) for s in L])
scen={'C 16 SNPs (3 copies)':design(base+[copies16[c] for c in copies16]),'B 13 SNPs':design(base),
 'A 8 loci':design(loci8),'D 13 + 3 copies of top-3 loci':design(base+[s for s in top3 for _ in range(3)])}
BX=np.array([bx[s] for s in base]); SY=np.array([sy[s] for s in base])
rng=np.random.default_rng(2026)
def run(D,N=400000,chunk=50000):
    k=len(D); x=BX[D]; se=SY[D]; rej=0
    for _ in range(N//chunk):
        y=(rng.standard_normal((chunk,len(base)))*SY)[:,D]
        b=y/x; w=1/(se/np.abs(x))**2; sw=w.sum()
        e=(w*b).sum(1)/sw; Q=(w*(b-e[:,None])**2).sum(1)
        t=e/np.sqrt(np.maximum(1,Q/(k-1))/sw)
        rej+=(2*stats.norm.sf(abs(t))<0.05).sum()
    return rej/N
print('top-3 loci:',top3)
for n,D in scen.items(): print(n,len(D),'rows, false-positive rate %.5f'%run(D))
