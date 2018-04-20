
# CountMappedReads - an in-house script for FPKM calculation

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](/LICENSE)
[![size](https://img.shields.io/github/size/webcaetano/craft/build/phaser-craft.min.js.svg)]()

CountMappedReads maps NGS reads on given reference sequences (by using bowtie2), then counts mapped reads and calculates FPKM.

## requirements
* bowtie2 (ver < 2.3.0)
* samtools
* Ruby (ver >=2.0)

## usage 
```
### CountMappedReads ver 1.0 ###

[description]
CountMappedReads maps NGS reads on given reference sequences (by using bowtie2), then counts mapped reads and calculates FPKM.

[usage]
$ CountMappedReads <reference fasta> <output dir> {-1 <pe1> -2 <pe2> | -U <up>} [options]

[arguments]
    - reference fasta      -- nucleotide fasta file for mapping/couting (e.g. genomes, contigs, genes)
    - output dir           -- output directory

[dependencies]
    - bowtie2 (ver < 2.3)
    - samtools
    - ruby (ver >=2.0)

[options]
  (general)
    -h, --help
    -v, --version

  (bowtie2)
    -1            [str] (required)    -- bowtie2 '-1' option (paired-end read 1)
    -2            [str] (required)    -- bowtie2 '-2' option (paired-end read 2)
    -U            [str] (required)    -- bowtie2 '-U' option (unpaired reads)
    --end-to-end        (default)     -- If specified, bowtie2 is run in 'end-to-end' alignment mode.
    --local                           -- If specified, bowtie2 is run in 'local' alignment mode.
    --score-min   [str]               -- bowtie2 '--score-min' option (default: 'L,0,-0.6' for 'end-to-end' mode, 'L,0,0.6' for 'local' mode)

  (computation)
    --threads     [int] (default: 1)  -- number of threads for run of bowtie2 and samtools
    --mem         [str] (default: 1G) -- maximum memory per thread for 'samtools sort'; suffix K/M/G is recognized (e.g., 800M)


[output files]
result/count-FPKM.tsv                 -- tab separated file (1: seq ID, 2: seq length, 3: FPKM, 4: mapped read count, 5: mapping type)

[note]
The '--score-min' option of bowtie2 controls a minimum acceptable level of mapping quality.
For example, with the 'end-to-end' mode, '--score-min L,0,-0.6' means an aligned region shows at least 90% identity, assuming that a read quality is high enough and a gap is not open.
Similarity, 'L,0,-0.3' means at least 95% identity, 'L,0,-0.9' means at least 85% identity, and 'L,0,-1.2' means at least 80% identity under the same assumption.
For deteils of the '--score-min' option and the alignment score, see bowtie2 manual page (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

If reference sequences are genes (rather than genomes), it might be appropriate to use bowtie2 'local' mode, considering a situation that only a part of a read is overlapping with a gene.
When calculation is performed by the 'local' mode, this tool performs a post-filtering process on the bowtie2 mapping, defined by two criteria (as follows) that cannot be adjusted by a user.
(1) Aligned length must be at least 50nt AND at least 0.5 * read length.
(2) Alignment score must be at least 1.2 * alignment length (i.e., suppose a read quality is high enough and a gap is not open, an aligned region shows at least 90% identity).
For the 'end-to-end' mode, no post-filtering process is performed (that is, the criteria above is not applied).
```
