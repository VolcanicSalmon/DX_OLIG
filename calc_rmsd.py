import os
import sys
import glob
import mdtraj as md
import numpy as np
import pandas as pd
import re
score_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
out_csv   = sys.argv[2] if len(sys.argv) > 2 else 'rmsd.csv'
wt_aas = {int(k): v for k, v in json.loads(sys.argv[3]).items()}
wt_refs = {}
for pos, aa in wt_aas.items():
    ref_pdb = os.path.join(score_dir, f"wt_{pos}_{aa}_wt_0001.pdb")
    if os.path.exists(ref_pdb):
        wt_refs[pos] = md.load(ref_pdb)
        print(f"Loaded WT ref pos {pos} ({aa}): {ref_pdb}")
    else:
        print(f"{ref_pdb} not found")

rows = []
for pdb in sorted(glob.glob(os.path.join(score_dir, '*_wt_0001.pdb'))):
    basename = os.path.basename(pdb)
    m = re.match(r'(.+?)_(\d+)_([A-Z]+)_wt_0001\.pdb', basename)
    if not m:
        continue
    prefix = m.group(1)
    pos = int(m.group(2))
    aa = m.group(3)

    # skip silent mut
    if pos in wt_aas and aa == wt_aas[pos]:
        continue
    if pos not in wt_refs:
        print(f"{basename}: no WT ref for position {pos}, skipping")
        continue
    try:
        ref = wt_refs[pos]
        t = md.load(pdb)
        ca_ref = ref.topology.select('name CA')
        ca_mut = t.topology.select('name CA')
        if len(ca_ref) != len(ca_mut):
            print(f"{basename}: mismatching")
            continue
        t.superpose(ref, atom_indices=ca_mut, ref_atom_indices=ca_ref)
        rmsd = md.rmsd(t, ref, atom_indices=ca_mut, ref_atom_indices=ca_ref)[0] * 10
        rows.append({'description': basename.replace('.pdb', ''), 'ca_rmsd': rmsd})
        print(f"{basename}: CA RMSD = {rmsd:.3f} Å")
    except Exception as e:
        print(f"{basename}: not found ({e})")
rmsd_df = pd.DataFrame(rows)
rmsd_df.to_csv(out_csv, index=False)
print(f"\nSaved: {out_csv} ({len(rmsd_df)} rows)")
