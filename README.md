# RANDOM-MUT

## Steps

The **first step** is to serialize your genome and store it as a python object.
This will actually generate a bigger file but that can be used by the software more easily  skipping a likely very long step.
The amount of time saved doing this step will depend in the IO speed.
It is convininent because you can serialize your refseq genome and then use the same to randomize more data sets.

```bash
GENOME="path/to/refseq.fa"
python random_genome_classic.py -M serialize -g $GENOME -a hg19
```

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
