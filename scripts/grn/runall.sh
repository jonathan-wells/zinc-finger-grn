#!/usr/bin/env bash

./extract_te_cres.sh
./extract_putative_tfs.sh
./extract_zfp_targets.sh
./cres_to_network.py
