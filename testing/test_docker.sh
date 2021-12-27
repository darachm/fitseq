# Run from the base repo directory, so bash testing/test_docker.sh
docker run \
    --mount type=bind,source=$(pwd)/testing,target=/testing \
    darachm/fitseq fitseq.py \
        -i testing/data/ppiseq_test_counts_1000.csv \
        -p 8 -t 0 1 2 3 4 \
        -m 20 --min-step 0.001 \
        --output-mean-fitness testing/output/test_means.csv \
        -o testing/output/test_out.csv
        #2> testing/output/test.err
