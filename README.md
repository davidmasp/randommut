# README

## Steps

The **first step** is to serialize your genome and store it as a python object.
This will actually generate a bigger file but that can be used by the software more easily  skipping a likely very long step.
The amount of time saving that this step will generate depends in the IO speed.
It is recomended because you can serialize your refseq genome and then use the same to randomize more data sets.

```bash
GENOME=path/to/refseq.fa
python random_genome_classic.py -M serialize -g $GENOME -a hg19
```

The **second step** is to generate the random positions.
The `times` argument (`-t`) is equivalent to the number of randomizations.
The `winlen` argument (`-w`) is equivalent to the length of the windows where the positions will be generated.

```bash
python random_genome_classic.py -M randomize -g chr1.fa.p -m PACA_AU_chr1_small.tsv -a hg19 -o outfile.test -t 50 -w 50000
```

### Check results are in the winlen range

Open the file with Excel.

Open conditional formating, custom rule

```txt
=AND(H2>($B2-50000);H2<($C2+50000))
```

Change `50000` with your winlen.