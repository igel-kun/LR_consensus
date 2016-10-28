#!/bin/bash

read -r -a DISTS < <(bin/dists "$1" "$2" | tail -n +2)

echo "RF: ${DISTS[0]}"
echo "Tr: ${DISTS[1]}"
echo "RL: ${DISTS[2]}"
echo "nRF: ${DISTS[3]}"
echo "nTr: ${DISTS[4]}"
echo "nRL: ${DISTS[5]}"
