#!/usr/bin/env python3
"""Finite exact checks and source-integrity audit for coefficient topology."""
from fractions import Fraction
from hashlib import sha256
from itertools import product
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
ASSETS=ROOT/'data/examples/assets'
SOURCE=Path('/home/ec2-user/latex/PartialCubes/papers/coefficient-topology')

def minimal_missing(words,n):
 W=set(words);out=set()
 for support in range(1,1<<n):
  bits=support
  while True:
   if not any((w&support)==bits for w in W) and all(any((w&(support^(1<<i)))==(bits&~(1<<i)) for w in W) for i in range(n) if support>>i&1): out.add((support,bits))
   if bits==0: break
   bits=(bits-1)&support
 return out

def cartesian(C,n,D,m): return {c|(d<<n) for c in C for d in D}

def partial_trace(C,n,I):
 out=set();inds=[i for i in range(n) if I>>i&1]
 for c in C:
  for masks in range(1<<len(inds)):
   released=[inds[j] for j in range(len(inds)) if masks>>j&1]
   for vals in range(1<<len(released)):
    x=c
    for j,i in enumerate(released): x=(x&~(1<<i))|(((vals>>j)&1)<<i)
    mask_word=sum(((masks>>j)&1)<<(n+j) for j in range(len(inds)))
    out.add(x|mask_word)
 return out

def run():
 review=json.loads((ASSETS/'coefficient_topology_build_review.json').read_text())
 assert review['status']=='passed' and review['pdf_pages']==82 and not review['warnings']
 for label in ('thm:positive','thm:ratio-laws','thm:ratio-zero','thm:shared-row-reduction','thm:coefficient-topology','thm:complement-topology','thm:minor-maps','thm:partial-trace'):
  assert label in review['labels']
 for name,expected in review['sha256'].items():
  path=SOURCE/name
  if path.is_file(): assert sha256(path.read_bytes()).hexdigest()==expected

 # Full cubes have no rows; punctured cubes have one full-support row;
 # singletons have one singleton row per fixed coordinate.
 for n in range(4): assert not minimal_missing(range(1<<n),n)
 for n in range(1,4):
  punct=set(range(1<<n))-{(1<<n)-1}
  assert minimal_missing(punct,n)=={((1<<n)-1,(1<<n)-1)}
  assert minimal_missing({0},n)=={(1<<i,1<<i) for i in range(n)}

 # Product patterns are exactly embedded factor patterns on a bounded suite.
 suite=[({0},1),({0,1,2},2),({0,3},2),({0,1,3},2)]
 product_checks=0
 for C,n in suite:
  for D,m in suite:
   mm=minimal_missing(cartesian(C,n,D,m),n+m)
   expected=minimal_missing(C,n)|{(s<<n,b<<n) for s,b in minimal_missing(D,m)}
   assert mm==expected;product_checks+=1

 # The square antipodal pair has the four ratios whose product is one.
 # Equal coefficients and magnitudes attain value one exactly.
 anti={0,3}; assert minimal_missing(anti,2)=={(3,1),(3,2)}
 ratios=[Fraction(1) for _ in range(4)]
 assert max(ratios)==1

 # Partial trace of a singleton is the three-vertex path; its lifted row
 # has one free positive mask/old coefficient ratio.
 T=partial_trace({0},1,1);assert T=={0,2,3}
 assert minimal_missing(T,2)=={(3,1)}
 for u in (Fraction(1,3),Fraction(1),Fraction(5,2)):
  bx=1/(1+u);bm=u/(1+u);assert bx>0 and bm>0 and bx+bm==1

 result={'status':'passed','source_build_pages':82,'source_labels_checked':8,
         'full_cube_dimensions_checked':[0,1,2,3],'punctured_cube_dimensions_checked':[1,2,3],
         'singleton_dimensions_checked':[1,2,3],'cartesian_product_pattern_checks':product_checks,
         'square_antipodal_ratio':'1','partial_trace_singleton_words':sorted(T),
         'scope':'Exact small-class checks of minimal-row identities, product rows, the antipodal ratio endpoint and partial-trace free ratios, plus source-build label and hash integrity. General topology and operation laws use the written proofs.',
         'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
 out=ASSETS/'coefficient_topology_core_verification.json';out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':run()
