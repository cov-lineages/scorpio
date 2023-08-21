# scorpio <img src="https://anaconda.org/bioconda/scorpio/badges/downloads.svg" alt="" align = "right"/>
(serious constellations of reoccurring phylogenetically-independent origin)

<img src="https://github.com/cov-lineages/scorpio/blob/main/docs/scorpio_logo.png" width="100">

Scorpio provides a set of command line utilities for classifying, haplotyping and defining constellations of mutations for an aligned set of genome sequences. It was developed to enable exploration and classification of variants of concern within the SARS-CoV-2 pandemic and all SARS-CoV-2 specific information can be installed via [constellations](https://github.com/cov-lineages/constellations).

## Wiki
For example commands and FAQ please checkout the [wiki](https://github.com/cov-lineages/scorpio/wiki).

## Installation

You can install scorpio from Bioconda:

`conda install -c bioconda scorpio`

You can also build the contents of this repository locally with:

```
git clone https://github.com/cov-lineages/scorpio.git
cd scorpio
conda env create -f environment.yml
conda activate scorpio
pip install .
```

If you want to check your local installation has been successful, you can install pytest and run the included tests:
```
pip install pytest
pytest .
```
Please note that scorpio installation will always clone the most up-to-date version of the constellations repository, and these tests have been designed to pass with these definitions. Running with older constellations versions is likely to cause the tests to fail.

## Commands

Scorpio currently includes the following commands:
1. `classify` - takes a set of lineage-defining constellations with rules and classifies sequences by them.
2. `haplotype` - takes a set of constellations and writes haplotypes (either as strings or individual columns).
3. `list` - print the `mrca_lineage` and `output_name` of constellations as a single column to stdout.
4. `define` - takes a CSV with a group column and a mutations column and extracts the common mutations within the group, optionally with reference to a specified outgroup

An overview and example commands for each of these can be found in [the wiki](https://github.com/cov-lineages/scorpio/wiki).

## Constellations
The JSON file for an individual constellation (in this case a lineage defining one) would look like this:
```json
{
	"name": "B.1.1.7",
	"description": "B.1.1.7 lineage defining mutations",
	"citation": "https://virological.org/t/563",
	"sites": [
		"nuc:C913T",
		"1ab:T1001I",
		"1ab:A1708D",
		"nuc:C5986T",
		"1ab:I2230T",
		"1ab:SGF3675-",
		"nuc:C14676T",
		"nuc:C15279T",
		"nuc:C16176T",
		"s:HV69-",
		"s:Y144-",
		"s:N501Y",
		"s:A570D",
		"s:P681H",
		"s:T716I",
		"s:S982A",
		"s:D1118H",
		"nuc:T26801C",
		"8:Q27*",
		"8:R52I",
		"8:Y73C",
		"N:D3L",
		"N:S235F"
	],
        "rules": {
                "min_alt": 4,
                "max_ref": 6,
        }
}
```

The general format of a mutation code is:
`gene`:[`ref`]`coordinates`[`alt`]
where `gene` is a gene code (or `nuc` for the genomic nucleotide sequence), `ref` is the nucleotide or amino acids in the reference, `alt` is the specific nucleotide or amino acid for the mutatant. Either of `ref` or `alt` can be missing if no specific state is required.

Rules can either specify [min|max]_[ref|alt|ambig|oth] OR the call required at a mutation e.g. "N:S235F": (not )[ref|alt|ambig|oth]

More information can be found [about constellations](https://github.com/cov-lineages/scorpio/wiki/What-does-a-valid-constellation-look-like%3F) and [mutation definitions](https://github.com/cov-lineages/scorpio/wiki/What-does-a-valid-mutation-site-look-like%3F) on the wiki.
