#!/usr/bin/env python3
"""Polybrutus D68: exact Euclidean labelled metric and chosen ideal weighted graph.
Conservative perpendicular unit extension of the immutable Polybrutus D67 v0.1.0.
No physical spacetime dimension or experimentally measured resonance is claimed.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter
from hashlib import sha256
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
Q=s.Rational
core=[(0,0,1),(Q(1,2),0,Q(1,2)),(-Q(1,4),Q(3,4),Q(1,2)),
      (-Q(1,4),-Q(3,4),Q(1,2)),(0,0,-1),(Q(1,2),0,-Q(1,2)),
      (-Q(1,2),0,-Q(1,2))]
n=73
pairs=list(combinations(range(n),2))
G=s.zeros(n)
for i in range(7):
    for j in range(7):
        G[i,j]=3*core[i][0]*core[j][0]+core[i][1]*core[j][1]+core[i][2]*core[j][2]
for i in range(7,12):
    for j in range(7,12):
        G[i,j]=1 if i==j else -Q(1,4)
for i in range(12,n):
    G[i,i]=1
assert all(G[i,i]==1 for i in range(n))
assert DomainMatrix.from_Matrix(G).to_field().rank()==68
D=2*s.ones(n)-2*G
for i in range(n):
    D[i,i]=0
assert len(pairs)==2628 and all(D[i,j]>0 for i,j in pairs)

parent_path=ROOT/'prior'/'d67'/'results'/'certificate_d67.json'
parent_archive=ROOT/'prior'/'d67'/'frozen_scientific_tag_v0.1.0.zip'
assert sha256(parent_archive.read_bytes()).hexdigest()=='0fa69a4d16dc016172e126d4fc426e73b733a4d456c2b8702c5455f92484675c'
from zipfile import ZipFile
with ZipFile(parent_archive) as _frozen:
    assert _frozen.testzip() is None
    assert _frozen.read('results/certificate_d67.json')==parent_path.read_bytes()
parent_bytes=parent_path.read_bytes()
assert sha256(parent_bytes).hexdigest()=='f639385724f75f5624484bce2df5cd48292ced2184bed55ba09db2e1e30ed629'
parent=json.loads(parent_bytes)
assert parent['status']=='PASS' and parent['dimension']==67 and parent['vertices']==72
former={(int(e['i'])-1,int(e['j'])-1):s.Rational(e['d2'])
        for e in parent['labelled_squared_distances']}
assert len(former)==2556
assert all(D[i,j]==v for (i,j),v in former.items())
assert all(D[i,72]==2 for i in range(72))
hist=Counter(4*D[i,j] for i,j in pairs)
expected=Counter({4:6,7:2,8:2597,9:3,10:10,12:6,13:2,16:2})
assert hist==expected and sum(hist.values())==2628
assert sum(D[i,j] for i,j in pairs)==Q(21071,4)

H=s.eye(n)-s.ones(n)/n
B=-H*D*H/2
assert B==H*G*H
assert DomainMatrix.from_Matrix(B).to_field().rank()==68
# The first and fifth canonical points are antipodal; origin lies in the affine hull.
assert G[0,4]==-1

L=s.zeros(n)
for i,j in pairs:
    w=1/D[i,j]
    L[i,i]+=w;L[j,j]+=w;L[i,j]-=w;L[j,i]-=w
assert L*s.ones(n,1)==s.zeros(n,1)
assert DomainMatrix.from_Matrix(L-s.Rational(73,2)*s.eye(n)).to_field().rank()==11
assert DomainMatrix.from_Matrix(L-Q(36)*s.eye(n)).to_field().rank()==69
old=s.zeros(72)
for i,j in combinations(range(72),2):
    w=1/D[i,j]
    old[i,i]+=w;old[j,j]+=w;old[i,j]-=w;old[j,i]-=w
assert L[:72,:72]-old==s.eye(72)/2
assert all(L[i,72]==-Q(1,2) for i in range(72))
assert L[72,72]==Q(36)

# Four adversarial negative controls.
damaged=D.copy();damaged[0,72]=damaged[72,0]=Q(9,4)
assert Counter(4*damaged[i,j] for i,j in pairs)!=hist
duplicate=G.copy();duplicate[72,:]=G[71,:];duplicate[:,72]=G[:,71];duplicate[72,72]=1
assert DomainMatrix.from_Matrix(duplicate).to_field().rank()==67
tilted_dot_to_old_core=Q(3,10)
assert Q(3,10)**2+Q(91,100)==1
assert tilted_dot_to_old_core!=0 and 2-2*tilted_dot_to_old_core!=D[0,72]
swapped=D.copy();swapped.row_swap(0,7);swapped.col_swap(0,7)
assert swapped!=D and Counter(4*swapped[i,j] for i,j in pairs)==hist

report={
  'schema':'POLYBRUTUS-D68-EXACT-CERTIFICATE-v1',
  'status':'PASS','dimension':68,'vertices':73,'pairs':2628,
  'affine_rank':68,'centered_gram_rank':68,'ambient_storage_coordinates':69,
  'radius':1,'parent':'Topbrutus/polybrutus-d67@v0.1.0',
  'parent_tag_commit':'dad596e38305066e4cb1207bb8fc947779eb445e',
  'parent_doi':'10.5281/zenodo.24432937',
  'parent_certificate_sha256':sha256(parent_bytes).hexdigest(),
  'parent_frozen_tag_archive_sha256':sha256(parent_archive.read_bytes()).hexdigest(),
  'parent_frozen_tag_embedded_as_verified_zip':True,
  'parent_distances_preserved':2556,'new_distances':72,
  'new_squared_distances':['2']*72,
  'histogram_k_4d2':{str(k):int(v) for k,v in sorted(hist.items())},
  'sum_squared_distances':'21071/4',
  'gram_psd_certificate':'B=-(H D H)/2=H G H; G=XX^T exact',
  'negative_controls':'4/4 PASS',
  'oscillator_model':'Chosen ideal equal-mass complete graph with w_ij=1/d_ij^2; no physical measurement',
  'eigenvalue_73_over_2_multiplicity':62,
  'eigenvalue_36_multiplicity':4,
  'spectral_recurrence':'chi73(z)=z*(z-73/2)*chi72(z-1/2)/(z-1/2)',
  'spectral_closed_form':'chi73(z)=z*(z-73/2)^62*(z-36)^4*f6(z-33), f6(z)=chiP7(z)/z',
  'novel_universal_theorem':'NOT CLAIMED','physical_discovery':'NOT CLAIMED',
  'human_peer_review':'NOT_DONE','sympy':s.__version__,
  'labelled_squared_distances':[{'i':i+1,'j':j+1,'d2':str(D[i,j]),'k':str(4*D[i,j])} for i,j in pairs]
}
out=ROOT/'results'/'certificate_d68.json'
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')
print('PRIMARY D68 PASS: 73 vertices, 2628 labelled pairs, affine rank 68; 2556/2556 previous pairs preserved')
print('IDEAL GRAPH PASS: eigenvalue 73/2 x62; eigenvalue 36 x4; 4/4 negative controls')
print('CERTIFICATE_SHA256='+sha256(out.read_bytes()).hexdigest())
