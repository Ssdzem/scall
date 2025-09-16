#!/usr/bin/env python3
import sys, re, pysam
from collections import defaultdict, Counter

# usage: label_umis_with_genes.py <stellarscope-updated.bam> <te.gtf> > umi_labels.tsv
bam_path, gtf_path = sys.argv[1], sys.argv[2]

# --- minimal TE interval index: chrom -> list of (start,end,te_id)
idx = {}
with open(gtf_path) as fh:
    for ln in fh:
        if not ln or ln[0] == '#': continue
        chrom,_,_,beg,end,_,_,_,attrs = ln.rstrip('\n').split('\t')
        m = dict(re.findall(r'(\S+)\s+"([^"]*)"', attrs))
        te_id = m.get('transcript_id') or m.get('gene_id') or m.get('gene_name') or m.get('ID') or m.get('Name')
        if not te_id: continue
        idx.setdefault(chrom, []).append((int(beg)-1, int(end), te_id))

def te_hits(chrom, s, e):
    for a,b,t in idx.get(chrom, ()):
        if a < e and s < b:  # overlap
            yield t

bam = pysam.AlignmentFile(bam_path, 'rb')

umi_te = defaultdict(Counter)   # (CB,UB) -> Counter({te_id: hits})
umi_gx = defaultdict(set)       # (CB,UB) -> set(gene_ids)
umi_gn = defaultdict(set)       # (CB,UB) -> set(gene_names)

for a in bam:  # sequential; no index required
    if a.is_unmapped or a.is_secondary or a.is_supplementary: continue
    try:
        cb = a.get_tag('CB'); ub = a.get_tag('UB')
    except KeyError:
        continue
    if a.has_tag('GX'):
        for gid in re.split(r'[;,]', str(a.get_tag('GX'))):
            gid = gid.strip()
            if gid: umi_gx[(cb,ub)].add(gid)
    if a.has_tag('GN'):
        for gname in re.split(r'[;,]', str(a.get_tag('GN'))):
            gname = gname.strip()
            if gname: umi_gn[(cb,ub)].add(gname)
    chrom = bam.get_reference_name(a.reference_id)
    for s,e in a.get_blocks():
        for te in te_hits(chrom, s, e):
            umi_te[(cb,ub)][te] += 1

print("CB\tUB\tTE_id\tgene_ids\tgene_names\tlabel")
seen = set(umi_te.keys()) | set(umi_gx.keys())
for key in seen:
    cb, ub = key
    te_id = ""
    if umi_te[key]:
        te_id, _ = max(umi_te[key].items(), key=lambda kv: (kv[1], kv[0]))  # majority TE
    gids = ";".join(sorted(umi_gx[key])) if umi_gx[key] else ""
    gns  = ";".join(sorted(umi_gn[key])) if umi_gn[key] else ""
    if te_id and gids:  label = "TE+Gene"
    elif te_id:         label = "TE-only"
    elif gids:          label = "Gene-only"
    else:               continue  # neither TE nor gene; skip
    print(f"{cb}\t{ub}\t{te_id}\t{gids}\t{gns}\t{label}")
