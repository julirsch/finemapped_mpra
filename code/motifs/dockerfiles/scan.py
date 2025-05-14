import sys
import pandas as pd
import twobitreader
import os
from tqdm import tqdm

window = 20

variants = pd.read_table(sys.argv[1], names="chr pos REF ALT".split(), header=None)

# load genome
hg38 = twobitreader.TwoBitFile("/motifs/hg38.2bit")
hg38cache = {ch: hg38[ch] for ch in variants.chr.unique()}


# score variants
def get_score(seqfile):
    os.system(f"/usr/local/bin/motif_liquidator /motifs/HOCOMOCOv11_core.meme {seqfile} -o motif_scores.out")
    return pd.read_table("motif_scores.out")


def subst(seq, refallele, altallele):
    assert seq[window - 1].upper() == refallele
    newseq = seq[:(window - 1)] + altallele + seq[window:]
    assert newseq[(window - 1)] == altallele[0]
    return newseq


def get_seqs(row):
    seq = hg38cache[row.chr][row.pos - window: row.pos + window].upper()
    return seq, subst(seq, row.REF, row.ALT)


with open("variant_seqs.fasta", "w") as seqfile:
    for ix, row in tqdm(variants.iterrows(), total=len(variants)):
        if len(row.REF) != 1:
            continue
        ref_seq, alt_seq = get_seqs(row)
        seqfile.write(f">{ix}_REF\n{ref_seq}\n")
        seqfile.write(f">{ix}_ALT\n{alt_seq}\n")

# extract and write score deltas
scores = get_score("variant_seqs.fasta")
scores.columns = "motif seqname start stop strand score pval qval matched".split()
scores['orig_row'] = scores['seqname'].str.split('_').str[0].astype(int)
scores['is_ref_alt'] = scores['seqname'].str.split('_').str[1]
scores_ref = scores.query("is_ref_alt == 'REF'").groupby("motif orig_row".split()).score.max()
scores_alt = scores.query("is_ref_alt == 'ALT'").groupby("motif orig_row".split()).score.max()
scores_ref.name = "refscore"
scores_alt.name = "altscore"
scores_ref_alt = pd.concat([scores_ref, scores_alt], axis=1).fillna(0)
scores_ref_alt["delta"] = scores_ref_alt.altscore - scores_ref_alt.refscore
for motif, grp in scores_ref_alt.reset_index().groupby("motif"):
    variants.loc[grp.orig_row, motif + "_delta"] = grp.delta.values
    variants.loc[grp.orig_row, motif + "_refscore"] = grp.refscore.values

variants.to_csv("variants_with_scores.txt", sep="\t", index=None)
