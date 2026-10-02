# Purpose: Retained external functional result summaries
# Manuscript/panels: Supplementary Table 7 linked data
# Reviewer comments: R2 Comment 4
# Execution requires the inputs and working-directory layout described in README.md; "02_Reviewer_Analyses/13_External_functional_summary.py"
"""Prepare auditable supplementary extracts without discarding complete outputs."""
from pathlib import Path
import csv, math

out = Path('33_External_Functional_Analysis_20260929')
res = out / 'results'
def read(name):
    with (res / name).open() as f:
        return list(csv.DictReader(f))
def write(name, rows):
    if not rows:
        return
    with (res / name).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

go = {'GO:0006260','GO:0098813','GO:0000819','GO:0045229','GO:0043062',
      'GO:0030198','GO:0002443','GO:0002449','GO:0002250'}
# Theory-led comparison list, not chosen according to significant results.
kegg = {'hsa04110','hsa03030','hsa04512','hsa04510','hsa04062','hsa04060',
        'hsa04660','hsa04662','hsa04650','hsa04612'}
gsea = read('GSEA_all_complete.csv')
write('GSEA_nine_discovery_GO_all_groups.csv', [r for r in gsea if r['pathway'] in go])
write('GSEA_related_KEGG_all_groups.csv', [r for r in gsea if r['pathway'] in kegg])
ora = read('ORA_top300_up_complete.csv')
write('ORA_top300_nine_GO_and_related_KEGG.csv', [r for r in ora if r['ID'] in go | kegg])
top = []
for ds in sorted({r['Dataset'] for r in gsea}):
    for scope in sorted({r['Scope'] for r in gsea}):
        for group in ['C1','C2','C3']:
            for coll in ['GO_BP','KEGG']:
                for direction in [1,-1]:
                    rr = [r for r in gsea if r['Dataset']==ds and r['Scope']==scope
                          and r['Group']==group and r['Collection']==coll
                          and r['padj']!='NA' and float(r['NES'])*direction > 0]
                    rr.sort(key=lambda r:(float(r['padj']),-abs(float(r['NES'])),r['pathway']))
                    top.extend(rr[:10])
write('GSEA_top10_each_direction.csv',top)
print('Complete outputs retained; fixed-program and top-10 summary extracts saved.')
