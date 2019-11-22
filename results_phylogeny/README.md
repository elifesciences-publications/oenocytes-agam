# Phylogenetic analysis of FA synthesis pathway

## Sequence retrieval

For each gene family, we define initial homologs from Pfam domain annotations, and retrieve the sequence for each individual domain:

* FAdesaturase: `FA_desaturase` (and we manually add `Anogam_AGAP003050-RA`)
* FAelongase: `ELO`
* FAreductase: `NAD_binding_4`
* FAsynthase: `ketoacyl-synt`
* FAdecarboxy: `p450`

## Alignment

For each gene family:

```bash
# alignment
mafft --localpair --reorder --maxiterate 1000  --thread 10 FAreductase.00.fasta > FAreductase.00.l.fasta
# trimming alignment with trimal
trimal -in FAreductase.00.l.fasta -out FAreductase.00.lt.fasta -automated1
```

## Phylogeny

For each gene family:

```bash
# phylogeny inference with best model testing and 1000 UFbootstraps
iqtree -s FAreductase.00.lt.fasta -m TEST -mset LG,WAG,JTT -nt 10 -bb 1000 -pre FAreductase.00.lt.iqt
```

Visualization with Figtree.

## Tools

* `iqtree` v1.6.10
* `mafft` v7.310
* `trimal` v1.4.rev22
