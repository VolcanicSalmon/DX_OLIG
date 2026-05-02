import os
import sys
import glob
import json
import pandas as pd
import re

score_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
out_csv   = sys.argv[2] if len(sys.argv) > 2 else 'ddg.csv'
wt_aas = {int(k): v for k, v in json.loads(sys.argv[3]).items()}

def parse_sc(sc_path):
    rows = []
    header = None
    with open(sc_path) as f:
        for line in f:
            if line.startswith('SCORE:') and ('total_score' in line or 'score' in line) and header is None:
                header = line.strip().split()[1:]
            elif line.startswith('SCORE:') and header:
                vals = line.strip().split()[1:]
                if len(vals) == len(header):
                    rows.append(dict(zip(header, vals)))
    return rows


wt_scores = {}
for pos, aa in wt_aas.items():
    sc_file = os.path.join(score_dir, f"wt_{pos}_{aa}_score.sc")
    if os.path.exists(sc_file):
        rows = parse_sc(sc_file)
        if rows:
            row = rows[0]
            wt_scores[pos] = float(row.get('total_score', row.get('score', 0)))
            print(f"WT pos {pos} ({aa}): {wt_scores[pos]:.3f}")
    else:
        print(f"WARNING: {sc_file} not found")

# load all mutant
rows = []
for sc in sorted(glob.glob(os.path.join(score_dir, '*.sc'))):
    basename = os.path.basename(sc)
    # skip wt_pos* files
    if re.match(r'wt_pos\d+', basename):
        continue
    for row in parse_sc(sc):
        desc = row.get('description', '')
        m = re.match(r'(.+?)_(\d+)_([A-Z]+)_', desc)
        if m:
            row['prefix'] = m.group(1)
            row['position'] = int(m.group(2))
            row['mutant_aa'] = m.group(3)
        else:
            row['prefix'] = desc
            row['position'] = None
            row['mutant_aa'] = None
        rows.append(row)

df = pd.DataFrame(rows)
if 'score' in df.columns and 'total_score' not in df.columns:
    df.rename(columns={'score': 'total_score'}, inplace=True)

score_cols = ['total_score', 'cart_bonded', 'fa_atr', 'fa_dun', 'fa_elec',
              'fa_intra_rep', 'fa_intra_sol_xover4', 'fa_rep', 'fa_sol',
              'hbond_bb_sc', 'hbond_lr_bb', 'hbond_sc', 'hbond_sr_bb',
              'lk_ball_wtd', 'omega', 'p_aa_pp', 'rama_prepro', 'ref', 'yhh_planarity']
for col in score_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

df['position']  = pd.to_numeric(df['position'], errors='coerce').astype('Int64')
df['ddG_total'] = df.apply(
    lambda r: r['total_score'] - wt_scores[r['position']]
    if pd.notna(r['position']) and r['position'] in wt_scores else None,
    axis=1
)
df['ddG_total'] = pd.to_numeric(df['ddG_total'], errors='coerce')
df['is_wt'] = df.apply(
    lambda r: r['mutant_aa'] == wt_aas.get(r['position'], None)
    if pd.notna(r['position']) else False,
    axis=1
)

df.to_csv(out_csv, index=False)
print(f"\nSaved: {out_csv} ({len(df)} rows)")
print(f"Positions covered: {sorted(df['position'].dropna().unique())}")
print(f"Silent mutations (ddG should be ~0):")
print(df[df['is_wt']][['description', 'position', 'mutant_aa', 'total_score', 'ddG_total']].to_string())
