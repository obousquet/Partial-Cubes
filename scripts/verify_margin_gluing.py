#!/usr/bin/env python3
"""Replay exact certificates and the obstruction for margin section gluing."""
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import importlib.util
import json

ROOT = Path(__file__).resolve().parents[1]
ASSET = ROOT / "data/examples/assets/margin_incompatible_extensions_verification.json"
spec = importlib.util.spec_from_file_location("cover", ROOT / "scripts/verify_cover_recursive.py")
cover = importlib.util.module_from_spec(spec); assert spec.loader; spec.loader.exec_module(cover)

def minimal_missing(words, n):
    words=set(words); result=set()
    for support in range(1,1<<n):
        sub=support
        while True:
            bits=sub
            if not any((w&support)==bits for w in words):
                if all(any((w&(support^(1<<i)))==(bits&~(1<<i)) for w in words)
                       for i in range(n) if support>>i&1):
                    result.add((support,bits))
            if sub==0: break
            sub=(sub-1)&support
    return result

def sign(word,i): return 1 if word>>i&1 else -1

def check_packet(packet):
    n=packet['n']; words=tuple(packet['words']); t=Fraction(packet['margin'])
    witnesses={int(c):tuple(map(Fraction,x)) for c,x in packet['witnesses'].items()}
    rows={(r['support'],r['bits']):{int(i):Fraction(x) for i,x in r['coefficients'].items()} for r in packet['rows']}
    assert set(rows)==minimal_missing(words,n)
    assert set(witnesses)==set(words)
    for c,x in witnesses.items():
        for i in range(n): assert t <= sign(c,i)*x[i] <= 1
    for (support,bits),row in rows.items():
        assert set(row)=={i for i in range(n) if support>>i&1}
        assert all(v>=0 for v in row.values()) and sum(row.values())==1
        for c,x in witnesses.items():
            assert sum(a*sign(bits,i)*x[i] for i,a in row.items()) <= -t
    return rows

def ample(words,n):
    W=set(words)
    return all(cover.shattered(W,n,s)==cover.strongly_shattered(W,n,s) for s in range(1<<n))

def run():
    report=json.loads(ASSET.read_text())
    assert report['status']=='passed' and report['arithmetic']=='fractions.Fraction; no optimizer'
    rows_a=check_packet(report['certificates']['A']); rows_b=check_packet(report['certificates']['B'])
    A=tuple(report['certificates']['A']['words']); B=tuple(report['certificates']['B']['words'])
    X=tuple(report['union']['words']); D=tuple(report['section_words_in_five_coordinates'])
    assert (len(A),len(B),len(D),len(X))==(25,25,17,33)
    for W,n in ((A,6),(B,6),(D,5),(X,7)):
        assert ample(W,n) and cover.is_partial_cube(W,n)
    assert all(c^(1<<i) in A for c in A for i in (0,1,2,3) if c>>i&1)
    assert all(c^(1<<i) in B for c in B for i in (0,1,2,3) if c>>i&1)
    assert all(c^(1<<i) in X for c in X for i in (0,1,2,3,6) if c>>i&1)
    def ratio(rows):
        p1=(49,49);p2=(50,50)
        return (rows[p2][5]/rows[p2][4])/(rows[p1][5]/rows[p1][4])
    assert ratio(rows_a)==Fraction(1,3) and ratio(rows_b)==3
    for pattern in ((49,49),(50,50),(56,8),(112,64)):
        assert pattern in minimal_missing(X,7)
    assert all(c in X for c in (42,25,97,82))
    result={'status':'passed','piece_sizes':[25,25],'section_size':17,'union_size':33,
            'all_ample_partial_cubes':True,'private_downward_closure':True,
            'piece_certificate_margins':['1/240','1/240'],'section_cycle_ratios':['1/3','3'],
            'union_margin_zero':'strict ratio cycle alpha_2 < beta_minus < alpha_1 < beta_plus < alpha_2',
            'scope':'Exact Fraction replay of both positive certificates, all minimal rows, ampleness, partial-cube metrics, private downward closure, and the four-pattern ratio obstruction.',
            'source_asset_sha256':sha256(ASSET.read_bytes()).hexdigest(),
            'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    out=ASSET.with_name('margin_gluing_verification.json');out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': run()
