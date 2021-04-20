# scorpio
serious constellations of reoccurring phylogenetically-independent origin

<img src="https://github.com/cov-lineages/scorpio/blob/main/docs/scorpio_logo.png" width="300">

## Command line options:

### commands
1. `classify` - takes a set of lineage-defining constellations with rules and classifies sequences by them.
2. `haplotype` - takes a set of constellations and writes haplotypes (either as strings or individual columns).
3. `report` - creates a report HTML for a constellation

### general options
* `-i`, `--input` - primary input file (usually the FASTA file)
* `-m`, `--metadata` - the metadata CSV file (required for some commands)
* `-o`, `--output` - the output file or path
* `-p`, `--prefix` - the output prefix (when multiple output files are being produced)
* `-c`, `--constellation` - a file of one or more constellations in JSON format (default to installed file from constellation github?)
* `-n`, `--names` - a list of constellation names to include from the file

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
	]
}
```

The general format of a mutation code is:
`gene`:[`ref`]`coordinates`[`alt`]
where `gene` is a gene code (or `nuc` for the genomic nucleotide sequence), `ref` is the nucleotide or amino acids in the reference, `alt` is the specific nucleotide or amino acid for the mutatant. Either of `ref` or `alt` can be missing if no specific state is required.
