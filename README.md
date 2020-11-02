# CosntructASCATFiles

This repository contains code for constructing the necessary input files for running ASCAT based on a dbSNP VCF file. This multi-threaded algorithm walks across the genome collecting information on SNPs present in the dbSNP file.

## How do I run it?

ConstructASCATFiles takes a dbSNP VCF file and one or more BAM files (comma-separated list) as input and extracts coverage and b-allele frequency (BAF) information for SNPs exceeding a predefined population threshold.

### Compile code

First clone the GitHub repository with git and run the following line:

```bash
make all
```
The binary 'constructascatfileswgs' is located at:

```bash
./dist/Release/GNU-Linux/constructascatfileswgs
```

The following parameters are available:

- -b/--bam-files*: One or more BAM files for which coverage and BAF file are constructed (comma-separated)
- -s/--sample-names*: Comma-separated list of sample names. Order of names must match BAM files.
- -S/--snp-file*: dbSNP VCF file
- -o/--output-prefix*: Prefix for output files
- -a/--min-alignment-score: Minimal alignment score (Default = 40)
- -B/--min-base-score: Minimal base score (Default = 30)
- -t/--threads: Number of threads to use (Default = 1)
- -g/--gmaf-threshold: Population minor allele frequency (Default = 0.01)
- -c/--count-duplicates: Should duplicates be included (optional)

*Dependencies*
- gcc 4.3+
- g++ 4.3+
