
`matrices` contains alignment environment data needed by the SWPS3 aligner library. 

`default_aligner_params.json` is the default config file containing parameters for alignment/clustering. 

The values are explained in the context of merging two clusters as follows:

The representative sequences (S1, S2) are aligned, and produce a score. 
Based on the alignment, S1 and S2 also have a number of uncovered amino acids (potentially zero, if for example S1 aligns completely within S2).
If the alignment score is above `min_score`, the clusters will be fully merged, or partially merged. 
Whether to fully merge depends on `min_full_merge_score` and `max_aa_uncovered`. 

If S1 has < `max_aa_uncovered` amino acids uncovered and the alignment score was > `min_full_merge_score`, the clusters will be fully merged (S1 into S2).
Otherwise, if S2 has < `max_aa_uncovered` amino acids uncovered and the alignment score was > `min_full_merge_score`, the clusters will be fully merged (S2 into S1).
Otherwise, a partial merge takes place where all sequences of S2's cluster are aligned against S1 and assigned to S1's cluster if the alignment exceeds `min_score`. 

The current default values were empirically derived and shown to produce "good" results. 

