# randomMut

RandomMut is a python package to randomize mutations within a defined window
size and preserving the trinucleotide context of the mutation.

## Install

Install from PyPi:

```bash
pip install randommut
```

or cloning this repo:

```bash
git clone URL
cd randommut
python setup.py install
```

If the install was succesful you should see the help message when executing:

```
randommut -h
```

## Usage

The **first step** is to serialize your genome and store it as a python object.
This will actually generate a bigger file but it does save time in the long run.
The amount of time saved doing this step will depend in the IO speed.
It is convininent because you can serialize your refseq genome and then use
the same to randomize more data sets.

```bash
GENOME="path/to/refseq.fa"
python -m randommut -M serialize -g $GENOME -a hg19
```

It should output the serialized version of the genome in the current directory.
In terms of memory, this processed has a peak of memory when writting the file
approx. at 30Gb for the hg19 assembly.

The **second step** is to generate the random positions. We will input the
mutations in table format (see below).

* The `times` argument (`-t`) is equivalent to the number of randomizations.
* The `winlen` argument (`-w`) is equivalent to the length of the windows where the positions will be generated.

```bash
GENOME="path/to/refseq.fa.p"
MUTS="path/to/muts.tsv"
OUTFILE="path/to/outfile.tsv"
python random_genome_classic.py -M randomize -g $GENOME -m $MUTS -a hg19 -o $OUTFILE -t 50 -w 50000
```

## Preprocessing

The input format is defined as a table (tsv format) with the folloing columns:

* chr: chromosome name (should match the chromsoma name in the fasta file)
* start: position start (1-based)
* end: position end
* ref: reference allele
* alt: alternative allele
* strand: Not used currently
* sample: Not relevant as each mutation is independent

```text
chr1    10  10  G   A   1   SAMPLE1
chr2    20  20  C   T   1   SAMPLE2
...
```

### Convert to table format

* From ICGC data release files.

```bash
ICGC="path/to/ICGC"
zcat $ICGC | awk 'BEGIN {FS="\t";OFS="\t"} $14~/single/{$9="chr"$9;print($9,$10,$11,$16,$17,$12,$5);}'
```

* From a vcf (unisample) we can use [bcftools](http://samtools.github.io/bcftools/bcftools.html)

```
bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT{0}\t1\tsampleA\n' file.vcf > "file.tsv"
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

Change `50000` with your inputed winlen parameter.

