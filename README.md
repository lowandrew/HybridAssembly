# Hybrid Assembly

This repo contains a wrapper for Unicycler - the `hybrid_assembly.py` script runs trimming and
error correction on illumina reads and then renames things nicely and cleans up Unicycler files that we don't
really care about. More details to come...

## Installation

You'll need `unicycler` installed, as well as the `bbmap` suite - easiest way to get these is via conda (assuming you have
already added the `bioconda` channel):

`conda install unicycler bbmap`

Then clone this repo and use the `hybrid_assembly.py` script, which will tell you all you need to know about it if you
 give it the `-h` flag.