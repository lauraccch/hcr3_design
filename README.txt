Output:

- prelim_probes.fa: all probe sequences initially designed, before BLAST filtering and before probe number limitation
- probes.fa: blasted, filtered and number-limited final probe sequences
-> both of the fastas conatin the parts of the probes that habridize to the target (with NN in the middle)
-> when manually blating, use blastn & sequences should show in opposite direction than gene direction

TSV Files: The results of BLAST, against the specified database (if there are 2 TSV files, then the probes were blasted against the organism that contains the target, as well as another organism that will be present in the sample)