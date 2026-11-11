"""Model probability intervals; no exact-zero imputation of implicitly low calls."""
def bounds(m,h,implicit=False):
    if implicit:
        assert m is None and h is None
        return {'m':(0.,12/256),'h':(0.,12/256),'C':(232/256,1.)}
    if m is None or h is None or not 0<=m<=255 or not 0<=h<=255 or m+h>256:return None
    return {'m':(m/256,(m+1)/256),'h':(h/256,(h+1)/256),'C':(max(0.,1-(m+h+2)/256),max(0.,1-(m+h)/256))}

def call(b,threshold=.8):
    if b is None:return 'unknown'
    candidates=[k for k,(lo,hi) in b.items() if all(lo>=b[j][1] for j in b if j!=k)]
    if len(candidates)!=1:return 'quantization_tie'
    k=candidates[0]
    return k if b[k][0]>=threshold else 'low_confidence'

def selfcheck():
    assert call(bounds(None,None,True),.8)=='C'
    assert call(bounds(None,None,True),.95)=='low_confidence'
    assert call(bounds(255,0),.95)=='m' and call(bounds(0,255),.95)=='h'
    assert call(bounds(200,12),.8)=='low_confidence'
    assert call(bounds(205,0),.8)=='m' and call(bounds(204,0),.8)=='low_confidence'
    assert call(bounds(128,128),.5)=='quantization_tie'
    assert call(bounds(None,128))=='unknown'
    assert bounds(255,255) is None
    return {'synthetic_encoding_checks':10,'scope':'implementation only; not biological validation'}
