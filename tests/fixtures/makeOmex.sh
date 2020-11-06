#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

examples=(BIOMD0000000028 BIOMD0000000297 BIOMD0000000297.edited)
for example in "${examples[@]}"
do
    cd "${SCRIPT_DIR}/${example}/"
    tempfile=tmp.omex
    zip -r "${tempfile}" ex1/ ex2/ manifest.xml
    mv "${tempfile}" "../${example}.omex"
done
