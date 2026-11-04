"""Certified targeted moments for adaptive first-block Krylov enrichment.

For a new exact dyadic polynomial seed a this builds the column
    C3 a = E^* M^3 E a
and, optionally, the new row/column of the C4 Gram matrix
    a_i^* E^* M^4 E a_j = <M^2 E a_i, M^2 E a_j>.
Only one parity is processed per invocation.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from collections import defaultdict
from functools import lru_cache
import argparse,json,time,gc
from flint import arb,arb_poly,ctx,arb_mat

POWERS=[(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11),(13,13),(16,2),(17,17),(19,19)]
N=864; DEG=1727

def rational(q): q=Q(q); return arb(q.numerator)/q.denominator
def lograt(q): return rational(q).log()
def endpoint(key): side,q=key; return arb(side)-lograt(q)*arb(2)/3

def load_seed(path):
    d=json.loads(Path(path).read_text()); bits=int(d['bits'])
    vals=[arb(int(x))/(1<<bits) for x in d['raw_Legendre_coefficients']]
    if len(vals)!=N: raise ValueError((path,len(vals)))
    return vals

def save_vec(values,path,meta,bits=700):
    rows=[];maxw=0
    for i,v in enumerate(values):
        scaled=v*(1<<bits); lo=int(scaled.lower().floor().unique_fmpz()); hi=int(scaled.upper().ceil().unique_fmpz())
        if hi<lo or hi-lo>4: raise ArithmeticError((Path(path).name,i,hi-lo,v.str(20)))
        rows.append([str(lo),str(hi)]);maxw=max(maxw,hi-lo)
    Path(path).write_text(json.dumps(meta|{'status':'PASS_INCREMENT_INTERVALS','bits':bits,'max_width_units':maxw,'rows':rows},separators=(',',':'))+'\n')

def load_mat(path):
    d=json.loads(Path(path).read_text()); den=1<<(d['bits']+1);M=arb_mat(d['dimension'],d['dimension'])
    for i,row in enumerate(d['rows']):
        for j,(lo,hi) in enumerate(row): M[i,j]=arb(int(lo)+int(hi))/den+arb(0,arb(int(hi)-int(lo))/den)
    return M,d

def save_mat(M,path,meta,bits=700):
    rows=[];maxw=0
    for i in range(M.nrows()):
        row=[]
        for j in range(M.ncols()):
            v=M[i,j]*(1<<bits);lo=int(v.lower().floor().unique_fmpz());hi=int(v.upper().ceil().unique_fmpz())
            if hi<lo or hi-lo>6: raise ArithmeticError((i,j,hi-lo,M[i,j].str(20)))
            row.append([str(lo),str(hi)]);maxw=max(maxw,hi-lo)
        rows.append(row)
    Path(path).write_text(json.dumps(meta|{'status':'PASS_C4_GRAM_INTERVALS','bits':bits,'dimension':M.nrows(),'max_width_units':maxw,'rows':rows},separators=(',',':'))+'\n')

def legendre_polynomial(coeff,p):
    # Match the 31 October operation order exactly.  This matters for a
    # bit-for-bit Arb interval regression even when the enclosed real
    # polynomial is mathematically identical.
    x=arb_poly([0,1]);prev=arb_poly([1]);cur=x
    f=(prev*coeff[0]) if p==0 else (cur*coeff[0])
    for deg in range(2,DEG+1):
        prev,cur=cur,((2*deg-1)*x*cur-(deg-1)*prev)*(arb(1)/deg)
        if deg%2==p: f += cur*coeff[deg//2]
    expected=1726+p
    if f.degree()!=expected: raise ArithmeticError((f.degree(),expected))
    return f

def paths(k):
    shifts=[(Q(1),arb(24)/5)]
    for n,p in POWERS:
        w=-arb(p).log()/arb(n).sqrt();shifts.extend([(Q(n),w),(Q(1,n),w)])
    groups=defaultdict(dict);active=0
    for path in product(shifts,repeat=k):
        q=Q(1);partial=[q];weight=arb(1)
        for rr,w in path: q*=rr;partial.append(q);weight*=w
        low=min(partial);high=max(partial);length=arb(3)-lograt(high/low)
        if length<0: continue
        if not length>0: raise ArithmeticError('zero support path')
        active+=1
        for key,ww in [((1,high),weight),((-1,low),-weight)]: groups[q][key]=groups[q].get(key,arb(0))+ww
    return groups,active

def build_c3(seed_path,p,out,precision=8192):
    start=time.monotonic();ctx.prec=precision;coeff=load_seed(seed_path);f=legendre_polynomial(coeff,p);groups,active=paths(3);selected=sorted(q for q in groups if q>=1)
    @lru_cache(maxsize=128)
    def primitive(key):
        x=endpoint(key);v=x;vals=[]
        for j in range(2*DEG+1):vals.append(v/(j+1));v*=x
        return arb_poly(vals)
    moments=[arb(0) for _ in range(DEG+1)]
    for ii,q in enumerate(selected):
        boundary=arb_poly([])
        for key,w in groups[q].items():boundary+=primitive(key)*w
        disp=lograt(q)*arb(2)/3;shifted=f(arb_poly([disp,1]));degree=shifted.degree();conv=arb_poly(list(reversed(shifted.coeffs())))*boundary
        fac=arb(3)/2*(1 if q==1 else 2)
        for j in range(p,DEG+1,2): moments[j]+=fac*conv[degree+j]
        if ii%16==0: print('C3',ii+1,'/',len(selected),'seconds',round(time.monotonic()-start,1),flush=True)
    vals=[];x=arb_poly([0,1]);prev=arb_poly([1]);cur=x
    for deg in range(DEG+1):
        if deg==0:P=prev
        elif deg==1:P=cur
        else: prev,cur=cur,((2*deg-1)*x*cur-(deg-1)*prev)*(arb(1)/deg);P=cur
        if deg%2==p: vals.append(sum((P[j]*moments[j] for j in range(p,deg+1,2)),arb(0)))
    if len(vals)!=N: raise ArithmeticError(len(vals))
    save_vec(vals,out,{'power':3,'parity':p,'precision_bits':ctx.prec,'nonempty_paths':active,'selected_displacements':len(selected),'seed':str(seed_path),'elapsed_seconds':round(time.monotonic()-start,3)})
    primitive.cache_clear();gc.collect();print('C3 PASS',out,flush=True)

def build_c4_row(seed_paths,p,old_gram,out,precision=8192):
    start=time.monotonic();ctx.prec=precision;seeds=[load_seed(x) for x in seed_paths];fs=[legendre_polynomial(x,p) for x in seeds];r=len(fs);new=r-1
    groups,active=paths(2);events={}
    for q in groups:
        for key in groups[q]: events.setdefault(key,[arb_poly([]) for _ in range(r)])
    keys=sorted(events,key=lambda k:float(endpoint(k)))
    if endpoint(keys[0])!=-1 or endpoint(keys[-1])!=1: raise ArithmeticError('endpoint coverage')
    allq=sorted(q for q in groups if q>=1)
    for ii,q in enumerate(allq):
        disp=lograt(q)*arb(2)/3
        positives=[f(arb_poly([disp,1])) for f in fs]
        for qq in ([q] if q==1 else [q,1/q]):
            for sidx in range(r):
                shifted=positives[sidx] if qq==q else arb_poly([v*((-1)**(p+j)) for j,v in enumerate(positives[sidx].coeffs())])
                if qq not in groups: continue
                for key,w in groups[qq].items(): events[key][sidx]-=shifted*w
        if ii%16==0: print('C4 shifts',ii+1,'/',len(allq),'seconds',round(time.monotonic()-start,1),flush=True)
    running=[arb_poly([]) for _ in range(r)];cross=[arb(0) for _ in range(r)]
    for i,key in enumerate(keys[:-1]):
        for sidx in range(r): running[sidx]+=events[key][sidx]
        left=endpoint(key);right=endpoint(keys[i+1])
        if right<=0: continue
        if left<0:left=arb(0)
        gn=running[new]
        for j in range(r):
            prim=(running[j]*gn).integral();cross[j]+=3*(prim(right)-prim(left))
        if i%16==0: print('C4 pieces',i+1,'/',len(keys)-1,'seconds',round(time.monotonic()-start,1),flush=True)
    M=arb_mat(r,r)
    if r>1:
        if old_gram is None: raise SystemExit('need --old-c4 for r>1')
        old,_=load_mat(old_gram)
        if old.nrows()!=r-1: raise SystemExit(('old gram dimension',old.nrows(),r-1))
        for i in range(r-1):
            for j in range(r-1): M[i,j]=old[i,j]
    for j in range(r):M[j,new]=M[new,j]=cross[j]
    save_mat(M,out,{'parity':p,'precision_bits':ctx.prec,'nonempty_M2_paths':active,'seeds':[str(x) for x in seed_paths],'elapsed_seconds':round(time.monotonic()-start,3)})
    print('C4 ROW PASS',out,[x.str(18) for x in cross],flush=True)

def main():
    pa=argparse.ArgumentParser();pa.add_argument('--parity',choices=['even','odd'],required=True);pa.add_argument('--seed',action='append',required=True);pa.add_argument('--old-c4',type=Path);pa.add_argument('--c3-out',type=Path,required=True);pa.add_argument('--c4-out',type=Path);pa.add_argument('--stage',choices=['c3','c4','all'],default='all');pa.add_argument('--precision',type=int,default=8192)
    a=pa.parse_args();p=int(a.parity=='odd')
    if a.stage in ('c3','all'):build_c3(a.seed[-1],p,a.c3_out,a.precision)
    if a.stage in ('c4','all'):
        if a.c4_out is None:raise SystemExit('--c4-out required')
        build_c4_row(a.seed,p,a.old_c4,a.c4_out,a.precision)
if __name__=='__main__':main()
