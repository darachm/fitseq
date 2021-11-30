# Run from the base repo directory, so bash testing/test.sh
singularity exec --bind PyFitSeq:/PyFitSeq ~/.singularity/fitseq-dev.sif \
    python3 /PyFitSeq/pyfitseq.py -i testing/data/ppiseq_test_counts.csv \
        -p 7 -t 0 1 2 3 4 -m 2 --min-step 0.0001 \
        --fitness_type m \
        -o testing/output/test_out \
        #2> testing/output/test.err
