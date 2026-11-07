"""Unit-cost global edit distances; preserve sampling and keep old metrics diagnostic."""
from functools import lru_cache
import math
import edlib

FIELDS=['edit_bases','max_length_bp','substitution_bases','insertion_bases','deletion_bases','indel_events','both_callable_columns','terminal_indel_bases','legacy_substitutions','legacy_gap_bases','legacy_gap_events','legacy_both_columns']

def columns(x,y,l,h):
    s=i=d=e=b=0;prev=''
    for u,v in zip(x[l:h],y[l:h]):
        if u!='-' and v!='-':b+=1;s+=u!=v;kind=''
        elif u=='-':i+=1;kind='I'
        else:d+=1;kind='D'
        if kind and kind!=prev:e+=1
        prev=kind
    return s,i,d,e,b

@lru_cache(maxsize=65536)
def one(a,b):
    a,b=sorted((a,b));assert a and b and set(a+b)<=set('ACGT')
    z=edlib.align(a,b,mode='NW',task='path');p=edlib.getNiceAlignment(z,a,b)
    x,y=list(p['query_aligned']),list(p['target_aligned'])
    anchors=[j for j,(u,v) in enumerate(zip(x,y)) if u!='-' and v!='-']
    assert anchors
    lo,hi=anchors[0],anchors[-1]+1
    # Same equivalent internal-gap convention as the original comparison; raw ED unchanged.
    for g,other in [(x,y),(y,x)]:
        j=lo+1
        while j<hi:
            if g[j]!='-':j+=1;continue
            k=j+1
            while k<hi and g[k]=='-':k+=1
            original_end=k
            while j>lo+1 and g[j-1] in 'ACGT' and g[j-1]==other[j-1]==other[k-1] and all(c in 'ACGT' for c in other[j:k]):
                g[k-1]=g[j-1];g[j-1]='-';j-=1;k-=1
            j=original_end
    assert ''.join(x).replace('-','')==a and ''.join(y).replace('-','')==b
    s,i,d,e,bc=columns(x,y,0,len(x));ls,li,ld,le,lb=columns(x,y,lo,hi)
    assert s+i+d==z['editDistance'];assert bc+i+d==len(x)
    return (z['editDistance'],max(len(a),len(b)),s,i,d,e,bc,len(x)-(hi-lo),ls,li+ld,le,lb)

class Totals:
    def __init__(self):self.values=[0]*len(FIELDS);self.n=0;self.distances=[];self.no_terminal=[];self.affected=0
    def add(self,v):
        self.n+=1;self.values=[x+y for x,y in zip(self.values,v)];self.distances.append(v[0]/v[1]);self.affected+=v[7]>0
        if v[7]==0:self.no_terminal.append(v)
    def row(self):
        out=dict(zip(FIELDS,self.values));ed,L,s,i,d,e,bc,term,ls,lg,le,lb=self.values
        nl=sum(v[1] for v in self.no_terminal);ne=sum(v[0] for v in self.no_terminal)
        out.update(pairs=self.n,normalized_edit_distance=ed/L if L else None,pair_equal_distance=math.fsum(self.distances)/self.n if self.n else None,event_sensitivity=(s+e)/(bc+e) if bc+e else None,legacy_per_base_difference=(ls+lg)/(lb+lg) if lb+lg else None,legacy_per_event_difference=(ls+le)/(lb+le) if lb+le else None,terminal_affected_pairs=self.affected,no_terminal_pairs=len(self.no_terminal),no_terminal_distance=ne/nl if nl else None)
        return out

def total(pairs):
    t=Totals()
    for a,b in pairs:t.add(one(a,b))
    return t.row()
