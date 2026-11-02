"""Independent piecewise-polynomial quadrature check of 24 A/B entries.
Eight Gauss nodes integrate every polynomial used here exactly (degree <=10).
"""
from pathlib import Path
from fractions import Fraction
import json
from flint import arb,ctx
from physical_even_trial import prime_powers
from certify_spectrum import load_source
OUT=Path(__file__).resolve().parent
def leg(k,x):
    if k==0:return arb(1)
    prev=arb(1);cur=x
    for j in range(2,k+1):prev,cur=cur,((2*j-1)*x*cur-(j-1)*prev)/j
    return cur
def main():
    ctx.prec=2048;L=arb(3)/2
    shifts=[(arb(0),arb(24)/5)]
    for n,p in prime_powers(Fraction(3,2)):
        d=arb(n).log();w=-arb(p).log()/arb(n).sqrt();shifts.extend([(d,w),(-d,w)])
    gauss=[arb.legendre_p_root(8,k,weight=True) for k in range(8)]
    checks=[]
    for p in [0,1]:
        A=load_source('A',p);B=load_source('B',p)
        for i,j in [(p,p),(p,p+2),(p,p+4),(p+2,p+2),(p+2,p+4),(p+4,p+4)]:
            va=arb(0);vb=arb(0)
            for ia,(da,wa) in enumerate(shifts):
                for db,wb in shifts:
                    left=max(-L,-L-da,-L-db);right=min(L,L-da,L-db)
                    if right<left:continue
                    assert right>left
                    center=(left+right)/2;half=(right-left)/2
                    value=sum((weight*leg(i,(center+half*node+da)/L)*leg(j,(center+half*node+db)/L) for node,weight in gauss),arb(0))*half
                    vb+=wa*wb*value
                    if ia==0:va+=wb*value
            for name,value,M in [('A',va,A),('B',vb,B)]:
                error=(value-M[i//2,j//2].mid()).abs_upper()
                assert error<arb(2)**-511,(name,i,j,error)
                checks.append({'matrix':name,'degrees':[i,j],'error':error.str(20)})
    assert len(checks)==24
    result={'status':'PASS_24_INDEPENDENT_LOW_DEGREE_ENTRIES','method':'8-node Gauss rule exact for polynomial degree <=15; actual degree <=10; every triple support intersection retained','checks':checks}
    (OUT/'low_source_crosscheck.json').write_text(json.dumps(result,indent=2)+'\n');print(result['status'])
if __name__=='__main__':main()
