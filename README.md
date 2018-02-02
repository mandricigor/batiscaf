# BATISCAF - BAd conTIg removal SCAFfolding


### General information

BATISCAF is a novel repeat aware scaffolding tool. 

The main steps of the algorithms are:

  1. Removal of "bad" contigs, i.e. contigs which are short and which are deemed to be repeated.
  2. Solving the trivial scaffolding problem on the set of reliable contigs
  3. Re-inserting the previously removed contigs into the scaffolds




### Software prerequisites


Before running BATISCAF make sure that the following software is installed on your Linux computer:

  * A short read aligner: [Bowtie](http://bowtie-bio.sourceforge.net/manual.shtml), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), [BWA](http://bio-bwa.sourceforge.net/), or another one
  * [IBM ILOG CPLEX >=12.7](https://www.ibm.com/analytics/data-science/prescriptive-analytics/cplex-optimizer) - optimizer (solving Integer Linear Programs and many other stuff)
  * [Python 2.7](https://www.python.org/download/releases/2.7/) with the following libraries:
    * [NetworkX >=2.0](https://networkx.github.io/)
    * [Biopython](http://biopython.org/)



### Running BATISCAF

  * First, get BATISCAF from [GitHub](https://github.com/mandricigor/batiscaf) by executing the following command:
  
  ```
  git clone https://github.com/mandricigor/batiscaf.git
  ```
  and change the local directory:
  ```
  cd batiscaf
  ```

  * Next, prepare alignment files in .SAM format for the paired-end reads. We recommend using Bowtie2:
  ```
  bowtie2-build -q -f $CONTIG_FILE $INDEX_FILE
  bowtie2 --quiet --no-hd --reorder -k 10 -q -p 10 -x $INDEX_FILE -U $READ1_FASTQ -S $SAM1_FILE
  bowtie2 --quiet --no-hd --reorder -k 10 -q -p 10 -x $INDEX_FILE -U $READ2_FASTQ -S $SAM2_FILE
  ```

  * Obtain the scaffolding graph in .graphml format:
  ```
  python scaffolding_graph.py -o OUTPUT_GRAPHML -c CONTIGS_FASTA -m1 MAPPINGS1 -m2 MAPPINGS2 -i INS_SIZE -p PAIR_MODE -s STD_DEV
  ```
  * Run BATISCAF:
  ```
  batiscaf.py --graphml SCAFFOLDING_GRAPH --fasta FASTA
  ```
  By default, the scaffolds will be written to the file output.scaffolds.batiscaf.fasta.



### Help

  * In order to get more help for the scaffolding_graph.py helper script, execute the command:
  ```
  python scaffolding_graph.py -h
  ```
  You should see the following output:
  ```
  usage: scaffolding_graph.py [-h] [-o OUTPUT_GRAPHML] [-c CONTIGS_FASTA] -m1
                            MAPPINGS1 -m2 MAPPINGS2 -i INS_SIZE -p PAIR_MODE
                            -s STD_DEV

BATISCAF scaffolding graph construction helper script. Produces a .graphml
file

optional arguments:
  -h, --help         show this help message and exit
  -o OUTPUT_GRAPHML  output graphml file
  -c CONTIGS_FASTA   fasta file with contigs
  -m1 MAPPINGS1      comma separated list of .sam files (first read in the
                     read pair)
  -m2 MAPPINGS2      comma separated list of .sam files (second read in the
                     read pair)
  -i INS_SIZE        insert sizes (comma separated values)
  -p PAIR_MODE       pair modes (fr - innie style -> <-, rf - outtie style <-
                     ->) (comma separated values)
  -s STD_DEV         libraries standard deviations (comma separated values)
  ```
  
  * In order to get more help for the batiscaf.py main script, execute the command:
  ```
  python batiscaf.py -h
  ```
  You should see the following output:
  ```
  usage: batiscaf.py [-h] --graphml SCAFFOLDING_GRAPH [--fasta FASTA]
                   [--filter_threshold FILTER_THRESHOLD] [--mst]

BATISCAF - BAd conTIg removal SCAFfolder.

optional arguments:
  -h, --help            show this help message and exit
  --fasta FASTA         output fasta file (default
                        output.scaffolds.batiscaf.fasta
  --filter_threshold FILTER_THRESHOLD
                        filter out edges with weight less than this value
                        (default 5)
  --mst                 find MST (minimum spanning tree) before running
                        BATISCAF algorithm

required arguments:
  --graphml SCAFFOLDING_GRAPH
                        scaffolding graph in graphml format
```

### Contact us

In order to get more information on BATISCAF, do not hesitate to contact [Igor Mandric (Georgia State University)](mailto:mandric.igor@gmail.com).
