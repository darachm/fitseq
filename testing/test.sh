# Run from the base repo directory, so bash testing/test.sh
singularity exec --bind PyFitSeq:/PyFitSeq ~/.singularity/fitseq-dev.sif \
    python3 /PyFitSeq/pyfitseq.py -i testing/data/ppiseq_test_counts_1000.csv \
        -p 8 -t 0 1 2 3 4 -m 100 --min-step 0.001 \
        --min-iter 1 \
        --max-chunk-size 100000 \
        --fitness_type m \
        --output-mean-fitness testing/output/test_means.csv \
        -o testing/output/test_out.csv
        #2> testing/output/test.err
