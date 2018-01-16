# BATISCAF
## BAd conTIg removal SCAFfolding


BATISCAF is a novel repeat aware scaffolding tool. 

The main steps of the algorithms are:

  1. Removal of "bad" contigs, i.e. contigs which are short and which are deemed to be repeated.
  2. Solving the trivial scaffolding problem on the set of reliable contigs
  3. Re-inserting the previously removed contigs into the scaffolds


BATISCAF requires the user to supply a scaffolding graph (in .graphml format) and (optionally) the output filename. The output represents the set of scaffolds written in .fasta format.

In order to run BATISCAF, the user must execute the following command:

```
python batiscaf.py --graph rhodobacter-97-200-stranded.graphml --fasta output.fasta
```
