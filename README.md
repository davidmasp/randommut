# RANDOM-MUT

## Install

How I do it.

I am in a anaconda session. So my libraries are looking there. ( I think )
However, this should work with any python 3 installation.

1. clone the repository
2. install using setup
3. you should be able to run the module script and therefore run the randomization

```bash
git clone URL
cd randommut
python setup.py install
cd ~ ## cd C:\Users\username\whatever
python -m randommut -h
```

In this last chunk if the install was succesful you should see the help message.

## Steps

The **first step** is to serialize your genome and store it as a python object.
This will actually generate a bigger file but that can be used by the software more easily  skipping a likely very long step.
The amount of time saved doing this step will depend in the IO speed.
It is convininent because you can serialize your refseq genome and then use the same to randomize more data sets.

```bash
GENOME="path/to/refseq.fa"
python -m randommut -M serialize -g $GENOME -a hg19
```

It should output the serialized version of the genome in the current directory. Memory wise this processed got a peak of memory when writting the file at 30Gb.

The **second step** is to generate the random positions.

* The `times` argument (`-t`) is equivalent to the number of randomizations.
* The `winlen` argument (`-w`) is equivalent to the length of the windows where the positions will be generated.

```bash
GENOME="path/to/refseq.fa.p"
MUTS="path/to/muts.tsv"
OUTFILE="path/to/outfile.tsv"
python random_genome_classic.py -M randomize -g $GENOME -m $MUTS -a hg19 -o $OUTFILE -t 50 -w 50000
```

## Preprocessing

The input format is defined as...

### Convert ICGC data in the inputed format

```bash
ICGC="path/to/ICGC"
zcat $ICGC | awk 'BEGIN {FS="\t";OFS="\t"} $14~/single/{$9="chr"$9;print($9,$10,$11,$16,$17,$12,$5);}'
```

## Post-processing

Sometimes do you need to perform small actions to adecuate the output.

### Remove NNN positions

```bash
awk '$8!~/N/ {print($0);}' output_file.tsv
```

### Check results are in the winlen range

Open the file with Excel.

Open conditional formating, custom rule

```txt
=AND(H2>($B2-50000);H2<($C2+50000))
```

Change `50000` with your winlen.

## Input format requirement

Should be a `.tsv` with the folloing columns:

* chr
* start
* end
* ref
* alt
* strand
* sample

```text
chr1    10  10  G   A   1   SAMPLE1
chr2    20  20  C   T   1   SAMPLE2
...
```
