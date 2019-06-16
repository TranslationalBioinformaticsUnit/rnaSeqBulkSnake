#!/bin/env python
import pandas as pd

samples = str.split(snakemake.params.samples, ' ')

t = pd.read_table(snakemake.input[0], index_col=0, usecols=[0, 1], header=None, skiprows=4)

counts = [pd.read_table(f, index_col=0, usecols=[0, 1], header=None, skiprows=4)
          for f in snakemake.input]

# for t, sample in zip(counts, snakemake.params.samples):
for t, sample in zip(counts, samples):
    t.columns = [sample]


matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
matrix.to_csv(snakemake.output[0], sep="\t")
