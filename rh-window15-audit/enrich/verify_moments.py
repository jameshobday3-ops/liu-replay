"""Cross-check the new moment route against certified A/B and scalar products."""
import json,argparse,time
from pathlib import Path
from flint import arb,arb_mat,ctx
from inputs import baseline_spectrum as base
from build_targeted_moments import load_intervals,proposals
OUT=Path(__file__).resolve().parent
def main():
    pa=argparse.ArgumentParser();pa.add_argument('--baseline',type=Path,required=True);pa.add_argument('--column2-only',action='store_true');args=pa.parse_args()
    ctx.prec=1024;base.OUT=args.baseline;coeff=proposals();result=[];start=time.monotonic()
    for p in [0,1]:
        A=base.load_source('A',p);B=base.load_source('B',p);a=arb_mat([[x] for x in coeff[p]])
        Aa=A*a;Ba=B*a;b=arb_mat([[Aa[i,0]*(4*i+2*p+1)/3] for i in range(864)])
        new=load_intervals(OUT/f'moment2_column_p{p}.json')
        assert len(new)==864 and all(Ba[i,0].contains(new[i]) for i in range(864))
        assert not Ba[0,0].contains(new[0]+1)
        norm2=(a.transpose()*Ba)[0,0];resnorm=norm2-(Aa.transpose()*b)[0,0];assert resnorm>0
        row={'parity':p,'all_second_moment_column_entries_contained':True,'negative_control':'REJECTED','residual_norm_squared':resnorm.str(60)}
        if not args.column2_only:
            third=load_intervals(OUT/f'moment3_column_p{p}.json');scalars=load_intervals(OUT/f'scalar_moments_p{p}.json')
            third_scalar=sum((coeff[p][i]*third[i] for i in range(864)),arb(0))
            assert norm2.contains(scalars[0]);assert third_scalar.overlaps(scalars[1]);assert not third_scalar.overlaps(scalars[1]+1)
            row.update({'direct_second_moment_scalar_contained':True,'third_moment_scalar_independent_routes_overlap':True,'third_moment_scalar':third_scalar.str(60),'fourth_moment_scalar':scalars[2].str(60)})
        result.append(row);print(row,flush=True)
    receipt={'status':'PASS_SECOND_MOMENT_COLUMN_CHECK' if args.column2_only else 'PASS_TARGETED_MOMENT_CROSSCHECKS','results':result,'elapsed_seconds':round(time.monotonic()-start,3)}
    (OUT/('column2_validation.json' if args.column2_only else 'moment_validation.json')).write_text(json.dumps(receipt,indent=2)+'\n')
if __name__=='__main__':main()
