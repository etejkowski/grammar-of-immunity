#!/bin/sh
# Fetch official IMMREP23 challenge data (https://github.com/justin-barton/IMMREP23)
set -e
base=https://raw.githubusercontent.com/justin-barton/IMMREP23/main/data
for f in test.csv solutions.csv sample_submission.csv VDJdb_paired_chain.csv; do
    curl -sL -o "$(dirname "$0")/immrep23_$f" "$base/$f"
done
echo "IMMREP23 data fetched."
