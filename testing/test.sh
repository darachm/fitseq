# Run from the base repo directory, so bash testing/test.sh
singularity exec --bind PyFitSeq:/PyFitSeq ~/.singularity/fitseq-dev.sif \
    python3 /PyFitSeq/pyfitseq.py -i testing/data/ppiseq_test_counts_100.csv \
<<<<<<< HEAD
        -p 8 -t 0 1 2 3 4 -m 100 --min-step 0.0000 \
        --min-iter 100 \
=======
        -p 1 -t 0 1 2 3 4 -m 50 --min-step 0.00000000 \
>>>>>>> 4e59a0e9d0b9553ee1698fb0269d7907554e6515
        --fitness_type m \
        -o testing/output/test_out \
        #2> testing/output/test.err
