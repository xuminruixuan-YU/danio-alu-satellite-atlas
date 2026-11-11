"""Reference-CpG projection, equivalent gates to the existing regional pipeline.

Vectorization changes computation, not the selection rule. Every returned parent
is retained; exact coordinate checks concern only returned AS>=0.8*primary hits.
"""
import csv,gzip,itertools,json,os,re,sys,time
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pysam
from plot_common import BASE, save, now
D=Path(os.environ.get('DANIO_PROJECTION_OUTPUT','genome_projection_output')).resolve()
ROOT=BASE
from legacy_locus_cpg_coordinates import compile_placement,qualify_cpg
from legacy_mod_probability_calls import bounds,call

FIELDS=['C','m','h','unknown','quantization_tie','low_confidence','below_Q10']

def guard():
    s=os.statvfs(D);free=s.f_bavail*s.f_frsize/2**30
    mem={k:int(v.split()[0]) for k,v in (l.split(':',1) for l in Path('/proc/meminfo').read_text().splitlines())}
    avail=mem['MemAvailable']/2**20
    if free<300 or avail<160:raise RuntimeError(f'Resource reserve reached: disk {free:.1f} GiB; available RAM {avail:.1f} GiB')

def rec(a):
    return dict(read=a.query_name,flag=int(a.flag),reference=a.reference_name or '',start0=a.reference_start,
                end0=a.reference_end or -1,mapq=int(a.mapping_quality),AS=a.get_tag('AS') if a.has_tag('AS') else None,CIGAR=a.cigarstring or '*')

def parent_rows(path,fmt):
    if fmt=='SAM':
        with open(path) as f:
            for line in f:
                if line.startswith('@'):continue
                v=line.rstrip('\n').split('\t')
                yield dict(read=v[0],flag=int(v[1]),reference='' if v[2]=='*' else v[2],start0=int(v[3])-1,
                           AS=next((int(z[5:]) for z in v[11:] if z.startswith('AS:i:')),None),CIGAR=v[5])
    else:
        with gzip.open(path,'rt') as f:
            for v in csv.DictReader(f,delimiter='\t'):
                v['flag']=int(v['flag']);v['start0']=int(v['start0']);v['AS']=int(v['AS']) if v['AS'] else None
                yield v

def records_with_parents(part,limit=None):
    groups=itertools.groupby(parent_rows(part['parents'],part['parents_format']),lambda x:x['read'])
    with pysam.AlignmentFile(part['bam'],'rb',threads=2) as bam:
        for n,a in enumerate(bam):
            if limit is not None and n>=limit:return
            assert not a.is_secondary and not a.is_supplementary
            name,pp=next(groups);assert name==a.query_name,(part['name'],name,a.query_name)
            yield a,list(pp)
    assert next(groups,None) is None

def mapped_positions(comp,q):
    out=np.full(len(q),-1,dtype=np.int64)
    if not comp.matches:return out
    a=np.asarray(comp.matches,dtype=np.int64)
    ix=np.searchsorted(a[:,0],q,side='right')-1;valid=ix>=0
    j=np.maximum(ix,0);valid &= q<a[j,1]
    out[valid]=a[j[valid],2]+q[valid]-a[j[valid],0]
    return out

def classify(m,h,q):
    z=np.full(len(m),3,dtype=np.uint8)
    implicit=(m==-1)&(h==-1);z[implicit]=0
    ok=(m>=0)&(h>=0)&(m<=255)&(h<=255)&(m+h<=256)
    cm=np.maximum(0,254-m-h);ch=np.maximum(0,256-m-h)
    wm=ok&(m>=h+1)&(m>=ch)
    wh=ok&(h>=m+1)&(h>=ch)
    wc=ok&(cm>=m+1)&(cm>=h+1)
    wins=wm.astype(int)+wh.astype(int)+wc.astype(int)
    z[ok & (wins!=1)]=4;z[ok & (wins==1)]=5
    z[wm&(m>=205)]=1;z[wh&(h>=205)]=2;z[wc&(cm>=205)]=0
    z[q<10]=6
    return z

def project(a,pp,pos_by_ref,diag,validate=False):
    diag['input_reads']+=1
    if a.is_unmapped:diag['first_failure_unmapped']+=1;return None
    if a.mapping_quality<20:diag['first_failure_MAPQ']+=1;return None
    if any(int(x['flag'])&2048 or 'H' in x['CIGAR'] for x in pp):diag['first_failure_supplementary_or_H']+=1;return None
    sec=sum(bool(x['flag']&256) and not x['flag']&2048 for x in pp)
    if sec>=100:diag['first_failure_secondary_cap']+=1;return None
    pri=[x for x in pp if not x['flag']&(256|2048)];assert len(pri)==1
    score=pri[0]['AS']
    if score is None or score<=0:diag['first_failure_nonpositive_AS']+=1;return None
    assert pri[0]['reference']==a.reference_name and pri[0]['start0']==a.reference_start and pri[0]['CIGAR']==a.cigarstring
    high=[x for x in pp if x['AS'] is not None and x['AS']>=.8*score]
    diag['returned_parent_records']+=len(pp);diag['high_scoring_returned_parents']+=len(high)
    diag['read_has_one_high_score' if len(high)==1 else 'read_has_multiple_high_scores']+=1
    seq=a.get_forward_sequence();quality=a.get_forward_qualities()
    if not a.has_tag('MM') or not a.has_tag('ML') or quality is None:diag['first_failure_missing_tags_or_quality']+=1;return None
    forms=[s.split(',')[0] for s in a.get_tag('MM').split(';') if s]
    if not ('C+m.' in forms and 'C+h.' in forms and (not a.has_tag('dx') or a.get_tag('dx')!=1)):
        diag['first_failure_unsupported_encoding_or_dx1']+=1;return None
    mods=a.modified_bases_forward
    if mods is None:raise ValueError('MM/ML decode failed '+a.query_name)
    b=np.frombuffer(seq.encode(),dtype=np.uint8);p=np.flatnonzero((b[:-1]==67)&(b[1:]==71))
    if not len(p):diag['eligible_no_read_CpG']+=1;return None
    compiled=[compile_placement(x) for x in high];primary=next(c for c,x in zip(compiled,high) if not x['flag']&(256|2048))
    qc=len(seq)-1-p if a.is_reverse else p;qg=qc-1 if a.is_reverse else qc+1
    pc=mapped_positions(primary,qc);pg=mapped_positions(primary,qg);c=np.minimum(pc,pg)
    valid=(pc>=0)&(pg>=0)&(np.abs(pc-pg)==1)
    reference_pos=pos_by_ref[a.reference_name]
    ix=np.searchsorted(reference_pos,c)
    if len(reference_pos):valid &= (ix<len(reference_pos))&(reference_pos[np.minimum(ix,len(reference_pos)-1)]==c)
    else:valid &= False
    p=p[valid];qc=qc[valid];qg=qg[valid];pc=pc[valid];pg=pg[valid];c=c[valid];ix=ix[valid]
    diag['read_CG_matching_reference_CG_before_coordinate']+=len(p)
    pre=ix.copy();stable=np.ones(len(p),dtype=bool)
    for comp in compiled:
        if comp.unusable_reason or comp.reference!=a.reference_name or comp.reverse!=primary.reverse or comp.query_length!=len(seq):stable[:]=False;break
        stable &= (mapped_positions(comp,qc)==pc)&(mapped_positions(comp,qg)==pg)
    if validate:
        old=np.array([qualify_cpg(int(t),a.reference_name,int(r),primary,compiled)=='position_concordant' for t,r in zip(p,c)],dtype=bool)
        assert np.array_equal(old,stable),(a.query_name,'coordinate equivalence failed')
    diag['coordinate_conflict_or_incomplete']+=int((~stable).sum())
    p=p[stable];ix=ix[stable];c=c[stable]
    diag['position_concordant_observations']+=len(p)
    mm=dict(mods.get(('C',0,'m'),[]));hh=dict(mods.get(('C',0,'h'),[]))
    m=np.fromiter((mm.get(int(x),-1) for x in p),dtype=np.int16,count=len(p));h=np.fromiter((hh.get(int(x),-1) for x in p),dtype=np.int16,count=len(p))
    q=np.frombuffer(quality,dtype=np.uint8);qmin=np.minimum(q[p],q[p+1]);states=classify(m,h,qmin)
    if validate:
        truth=np.array([FIELDS.index(call(bounds(mm.get(int(x)),hh.get(int(x)),int(x) not in mm and int(x) not in hh),.8)) if quality[x]>=10 and quality[x+1]>=10 else 6 for x in p],dtype=np.uint8)
        assert np.array_equal(states,truth),(a.query_name,'modification equivalence failed')
    assert len(np.unique(ix))==len(ix),('duplicate read-reference-CpG',a.query_name)
    diag['eligible_reads']+=1;diag['valid_observations']+=int((states<3).sum())
    return a.reference_name,ix,states,pre,c

class CounterSet:
    def __init__(self,config,run):
        self.config=config;self.run=run;self.diag=Counter();self.pos={};self.counts={};self.pre={};self.reads={}
        self.scales=config['window_bp'];self.n_by_ref=Counter();self.part=Counter()
        self.arrays=json.loads((D/'inputs/array_windows.json').read_text())['matrices']
        self.features=defaultdict(list);self.feature_reads=Counter();self.catalog_reads=Counter();self.labels={}
        for x in self.arrays:
            a=x['array'];self.features[a['seqid']].append((a['array_id'],int(a['start0']),int(a['end0'])))
        for x in config['reference_records']:
            r=x['reference'];pos=np.load(D/'inputs/reference_CpGs'/(r+'.npy')).astype(np.int64,copy=False);self.pos[r]=pos
            self.counts[r]=np.zeros((len(pos),7),dtype=np.uint32);self.pre[r]=np.zeros(len(pos),dtype=np.uint32)
            self.labels[r]=np.zeros(len(pos),dtype=np.uint8)
            self.reads[r]={w:np.zeros((x['length']+w-1)//w,dtype=np.uint32) for w in self.scales}
        with (D/'inputs/fragments.tsv').open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                if row.get('species','TU' if row.get('fragment','').startswith('TU.') else '')!='TU':continue
                ref=row['seqid'];lo=int(row['start0']);hi=int(row['end0']);family=row['family']
                assert family in ['CL1','CL2'];bit=1 if family=='CL1' else 2
                a,b=np.searchsorted(self.pos[ref],[lo,hi-1]);self.labels[ref][a:b] |= bit
        f=config['focus']['array'];self.focus=f
        self.focus_path=D/'per_run'/run/'focus_valid_observations.tsv.gz'
        resume=(D/'per_run'/run/'reuse_checkpoint/receipt.json').exists()
        if resume:
            checkpoint=json.loads((D/'per_run'/run/'reuse_checkpoint/receipt.json').read_text())
            if 'focus_file_bytes' in checkpoint:
                with self.focus_path.open('r+b') as f:f.truncate(checkpoint['focus_file_bytes'])
        self.focus_out=gzip.open(self.focus_path,'at' if resume else 'wt',compresslevel=1)
        if not resume:self.focus_out.write('read\treference\tposition0\tcall\tMAPQ\n')

    def add(self,a,parents,part,validate=False):
        out=project(a,parents,self.pos,self.diag,validate)
        self.part[part+'_reads']+=1
        if out is None:return
        ref,ix,state,pre,c=out
        self.pre[ref][pre]+=1
        self.counts[ref][ix,state]+=1
        valid=state<3;positions=c[valid];calls=state[valid]
        if len(positions):
            self.n_by_ref[ref]+=1;self.part[part+'_contributing_reads']+=1
            for label in np.unique(self.labels[ref][ix[valid]]):self.catalog_reads[f'{ref}|{int(label)}']+=1
            for w in self.scales:self.reads[ref][w][np.unique(positions//w)]+=1
            for key,lo,hi in self.features[ref]:
                if np.any((positions>=lo)&(positions+1<hi)):self.feature_reads[key]+=1
            f=self.focus
            if ref==f['seqid']:
                mask=(positions>=int(f['plot_start0']))&(positions+1<int(f['plot_end0']))
                for p,s in zip(positions[mask],calls[mask]):self.focus_out.write(f'{a.query_name}\t{ref}\t{p}\t{FIELDS[s]}\t{a.mapping_quality}\n')
        self.part[part+'_valid_observations']+=len(positions)

    def checkpoint(self,stage,extra=None):
        save(D/'per_run'/self.run/'progress.json',dict(at_HKT=now(),stage=stage,diagnostics=dict(self.diag),parts=dict(self.part),**(extra or {})))

    def cache_reuse(self):
        dest=D/'per_run'/self.run/'reuse_checkpoint';dest.mkdir(exist_ok=True)
        self.focus_out.close()
        for ref in self.pos:
            np.savez_compressed(dest/(ref+'.npz'),counts=self.counts[ref],pre_coordinate=self.pre[ref],
                                **{f'reads_{w}':v for w,v in self.reads[ref].items()})
        save(dest/'receipt.json',dict(passed=True,at_HKT=now(),diagnostics=dict(self.diag),parts=dict(self.part),
             contributing_read_IDs_by_reference=dict(self.n_by_ref),large_region_distinct_read_IDs=dict(self.feature_reads),catalog_distinct_reads=dict(self.catalog_reads),focus_file_bytes=self.focus_path.stat().st_size))
        self.focus_out=gzip.open(self.focus_path,'at',compresslevel=1)

    def load_reuse(self):
        dest=D/'per_run'/self.run/'reuse_checkpoint';state=json.loads((dest/'receipt.json').read_text())
        self.diag=Counter(state['diagnostics']);self.part=Counter(state['parts'])
        self.n_by_ref=Counter(state['contributing_read_IDs_by_reference']);self.feature_reads=Counter(state['large_region_distinct_read_IDs'])
        self.catalog_reads=Counter(state['catalog_distinct_reads'])
        for ref in self.pos:
            z=np.load(dest/(ref+'.npz'));self.counts[ref]=z['counts'];self.pre[ref]=z['pre_coordinate']
            for w in self.scales:self.reads[ref][w]=z[f'reads_{w}']

    def output(self):
        self.focus_out.close();dest=D/'per_run'/self.run/'counts';dest.mkdir(exist_ok=True)
        for ref in self.pos:
            np.savez_compressed(dest/(ref+'.npz'),counts=self.counts[ref],pre_coordinate=self.pre[ref],
                                **{f'reads_{w}':v for w,v in self.reads[ref].items()})
        save(D/'per_run'/self.run/'projection_receipt.json',dict(passed=True,at_HKT=now(),run=self.run,
             diagnostics=dict(self.diag),parts=dict(self.part),contributing_read_IDs_by_reference=dict(self.n_by_ref),
             large_region_distinct_read_IDs=dict(self.feature_reads),catalog_distinct_reads=dict(self.catalog_reads),scope='current reference projection filters; returned candidates only'))

def selfcheck(config):
    m,h=np.meshgrid(np.arange(256,dtype=np.int16),np.arange(256,dtype=np.int16));m=m.ravel();h=h.ravel()
    actual=classify(m,h,np.full(len(m),30,dtype=np.uint8))
    expected=np.array([FIELDS.index(call(bounds(int(a),int(b)),.8)) for a,b in zip(m,h)],dtype=np.uint8)
    assert np.array_equal(actual,expected)
    for a,b in [(-1,-1),(-1,200),(200,-1)]:
        aa=np.array([a],dtype=np.int16);bb=np.array([b],dtype=np.int16)
        assert classify(aa,bb,np.array([30]))[0]==FIELDS.index(call(bounds(None if a<0 else a,None if b<0 else b,a<0 and b<0),.8))
    item=next(x for x in config['items'] if x['run']=='ERR16782316');positions={r['reference']:np.load(D/'inputs/reference_CpGs'/(r['reference']+'.npy')).astype(np.int64,copy=False) for r in config['reference_records']}
    diag=Counter();n=0
    for a,pp in records_with_parents(item['parts'][0],limit=1000):project(a,pp,positions,diag,True);n+=1
    save(D/'checks/vector_equivalence.json',dict(passed=True,all_ML_byte_pairs=65536,implicit_unknown_cases=3,
         real_input_records=n,real_coordinate_and_state_checks=dict(diag),reference='legacy regional coordinate and probability functions',analysis_parameters_changed=False))

if __name__=='__main__':selfcheck(json.loads((D/'inputs/config.json').read_text()))
