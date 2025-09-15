#!/usr/bin/env python3
import sys, re, pysam
from collections import defaultdict, Counter

# usage: label_te_umis.py <stellarscope-updated.bam> <te.gtf> > umi_labels.tsv

bam_path, gtf_path = sys.argv[1], sys.argv[2]

# --- build a minimal TE interval index: chrom -> list of (start,end,te_id)
idx = {}
with open(gtf_path) as fh:
    for ln in fh:
        if not ln or ln[0] == '#': continue
        chrom, _, _, beg, end, _, _, _, attrs = ln.rstrip('\n').split('\t')
        m = dict(re.findall(r'(\S+)\s+"([^"]*)"', attrs))
        te_id = m.get('transcript_id') or m.get('gene_id') or m.get('gene_name') or m.get('ID') or m.get('Name')
        if not te_id: continue
        idx.setdefault(chrom, []).append((int(beg)-1, int(end), te_id))

# quick overlap helper (linear scan is fine for typical TE density; swap for intervaltree if huge)
def hit_te(chrom, s, e):
    for a,b,t in idx.get(chrom, ()):
        if a < e and s < b:  # overlaps
            yield t

# --- walk BAM and collect per-(CB,UB) evidence
bam = pysam.AlignmentFile(bam_path, 'rb')
umi_has_gene = defaultdict(bool)
umi_te_hits  = defaultdict(Counter)

for a in bam:
    if a.is_unmapped or a.is_secondary or a.is_supplementary: continue
    try:
        cb = a.get_tag('CB'); ub = a.get_tag('UB')
    except KeyError:
        continue
    if a.has_tag('GX') and a.get_tag('GX') not in ('', '-', None):  # gene present
        umi_has_gene[(cb,ub)] = True

    chrom = bam.get_reference_name(a.reference_id)
    for s,e in a.get_blocks():  # aligned blocks
        for te in hit_te(chrom, s, e):
            umi_te_hits[(cb,ub)][te] += 1

# --- emit labels: one row per molecule that overlaps a TE
print("CB\tUB\tTE_id\tlabel")
for (cb,ub), cnt in umi_te_hits.items():
    if not cnt: continue
    te_id, _ = max(cnt.items(), key=lambda kv: (kv[1], kv[0]))  # majority vote TE
    lab = "TE+Gene" if umi_has_gene[(cb,ub)] else "TE-only"
    print(f"{cb}\t{ub}\t{te_id}\t{lab}")
