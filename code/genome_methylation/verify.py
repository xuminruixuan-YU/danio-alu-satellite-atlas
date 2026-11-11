"""Recalculate deposited summary arithmetic and check frozen plotting inputs.

This is not a BAM replay or a fresh full per-CpG decoding audit.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd

BASE=Path(__file__).resolve().parents[2];DATA=BASE/'data/genome_methylation'
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);a=parser.parse_args()
    if a.output.exists() and any(a.output.iterdir()):raise ValueError('Choose a new or empty output directory')
    a.output.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((DATA/'parameters.json').read_text());assert cfg['reference']=='GCF_049306965.1' if 'reference' in cfg else cfg['assembly']=='GCF_049306965.1'
    genomes=pd.read_csv(DATA/'genome_categories.tsv',sep='\t')
    chromosomes=pd.read_csv(DATA/'chromosomes.tsv',sep='\t')
    catalogue=pd.read_csv(DATA/'catalog_by_chromosome.tsv',sep='\t')
    for filename in ['genome_categories.tsv','chromosomes.tsv','catalog_by_chromosome.tsv','large_regions.tsv','S4_13_fixed_200kb_backgrounds.tsv']:
        t=pd.read_csv(DATA/filename,sep='\t');ok=t.effective_CpGs>0
        assert np.array_equal(t.reference_CpGs,t.effective_CpGs+t.missing_CpGs)
        assert np.array_equal(t.valid_observations,t.C_observations+t.m_observations+t.h_observations)
        assert np.allclose(t.loc[ok,'CpG_equal_5mC'],t.loc[ok,'CpG_fraction_sum']/t.loc[ok,'effective_CpGs'])
        assert np.allclose(t.loc[ok,'observation_weighted_5mC'],t.loc[ok,'m_observations']/t.loc[ok,'valid_observations'])
        assert t.loc[~ok,'CpG_equal_5mC'].isna().all()
    checks=[]
    for item in cfg['items']:
        run=item['run'];c=chromosomes[chromosomes.run==run].set_index('reference')
        assert len(c)==26 and int(c.nuclear.sum())==25
        w=pd.read_csv(DATA/f'{run}.windows.tsv.gz',sep='\t',low_memory=False)
        old=w[w.window_bp!=200000];new=w[w.window_bp==200000]
        assert np.array_equal(old.display,(old.effective_CpGs>=10)&(old.independent_reads>=3))
        support=(new.distinct_reads_lower_bound>=3)|(new.distinct_reads_exact>=3)
        assert np.array_equal(new.display,(new.effective_CpGs>=10)&support)
        assert w.loc[w.effective_CpGs==0,'CpG_equal_5mC'].isna().all()
        for (ref,scale),q in w.groupby(['reference','window_bp']):
            q=q.sort_values('start0');r=c.loc[ref]
            assert q.start0.iloc[0]==0 and np.array_equal(q.start0.to_numpy()[1:],q.end0.to_numpy()[:-1])
            for field in ['reference_CpGs','effective_CpGs','valid_observations','C_observations','m_observations','h_observations']:
                assert int(q[field].sum())==int(r[field]),(run,ref,scale,field)
            s=float((q.CpG_equal_5mC*q.effective_CpGs).sum())
            assert np.isclose(s/r.effective_CpGs,r.CpG_equal_5mC,rtol=1e-10)
        # Annotation categories partition each reference sequence, including missing positions.
        for ref,r in c.iterrows():
            cat=catalogue[(catalogue.run==run)&(catalogue.reference==ref)]
            for field in ['reference_CpGs','effective_CpGs','valid_observations']:
                assert int(cat[field].sum())==int(r[field])
        checks.append(dict(run=run,reference_records=26,windows=len(w),passed=True))
    overview=[]
    for scope,g in genomes.groupby('scope'):
        valid=g[g.effective_CpGs>0]
        if valid.empty:continue
        overview.append(dict(scope=scope,individuals=len(valid),CpG_equal_mean=valid.CpG_equal_5mC.mean(),observation_weighted_mean=valid.observation_weighted_5mC.mean()))
    pd.DataFrame(overview).to_csv(a.output/'recomputed_individual_equal.tsv',sep='\t',index=False)
    result={'passed':True,'scope':'deposited summary arithmetic, non-overlapping windows and plot gates; no new raw projection',
            'samples':checks,'missing_not_imputed_zero':True,'200kb_read_bounds_not_mislabelled_exact':True}
    (a.output/'checked.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
