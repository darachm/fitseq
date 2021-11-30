singularity exec --bind ./PyFitSeq:/PyFitSeq ~/.singularity/fitseq-dev.sif \
    python3 /PyFitSeq/pyfitseq.py -i ppiseq_test_counts.csv \
        -t 0 1 2 3 4 -m 50 --min-step 0.0001 \
        -o test_out \
        2>&1 > test.err
