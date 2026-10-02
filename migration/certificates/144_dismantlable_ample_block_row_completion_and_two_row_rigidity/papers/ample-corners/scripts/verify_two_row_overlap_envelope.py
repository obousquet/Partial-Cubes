"""Recover the envelope independently, rebuild the exact CNF and replay DRAT."""
import gzip, hashlib, json, subprocess, tempfile
from pathlib import Path
from pysat.formula import CNF
from explore_seed267_subfamilies import build, ROOT
from verify_full_base_defect import shadow

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    out=ROOT/'papers/ample-corners/evidence/two_row_overlap_envelope'
    r=json.loads((out/'report.json').read_text())
    source=ROOT/'papers/dismantlable-ample/evidence/order267_obstruction/record.json'
    assert digest(source)==r['source_sha256']
    X=set(json.loads(source.read_text())['family'])
    # Reduce coordinate zero by retaining edges, then contract coordinate one.
    C={x>>2 for x in X if x^1 in X}
    E={2*x+b for x in C for b in (0,1)}
    # Relabel the fresh least-significant coordinate to position ten.
    E={(x>>1)|((x&1)<<10) for x in E}
    assert C==set(r['base']) and E==set(r['envelope'])
    assert len(C)==len(shadow(C,10))==48 and len(E)==96
    pool,clauses,p=build(E,11)
    assert pool.top==r['variables'] and len(clauses)==r['clauses']
    assert r['status']=='UNSAT-DRAT-verified'
    with tempfile.TemporaryDirectory(dir='/home/ec2-user/damp-validation') as tmp:
        cnf=Path(tmp)/'search.cnf';CNF(from_clauses=clauses).to_file(str(cnf))
        assert digest(cnf)==digest(out/'search.cnf')==r['cnf_sha256']
        proof=Path(tmp)/'search.drat'
        with gzip.open(out/'search.drat.gz','rb') as f:proof.write_bytes(f.read())
        assert digest(proof)==r['proof_sha256']
        with (out/'replay.log').open('w') as f:
            cp=subprocess.run(['/home/ec2-user/damp-proof-archive/checkers/drat-trim',str(cnf),str(proof),'-t','45'],stdout=f,stderr=subprocess.STDOUT,timeout=55)
        assert cp.returncode==0 and 's VERIFIED' in (out/'replay.log').read_text()
    result=dict(status='PASS',base_order=48,envelope_order=96,formula_reconstructed=True,proof_replayed=True,scope='Every ample subfamily of this labelled envelope; no claim about all 96-word ample classes or other bases.')
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
