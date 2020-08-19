#!/bin/bash
cd BIOMD0000000028/
zip -r test.omex ex1/ ex2/ manifest.xml
mv test.omex  ../BIOMD0000000028.omex