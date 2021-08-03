# scorpio
serious constellations of reoccurring phylogenetically-independent origin

<img src="https://github.com/cov-lineages/scorpio/blob/main/docs/scorpio_logo.png" width="300">

## Command line options:

### commands
1. `classify` - takes a set of lineage-defining constellations with rules and classifies sequences by them.
2. `haplotype` - takes a set of constellations and writes haplotypes (either as strings or individual columns).
3. `report` - creates a report HTML for a constellation
4. `define` - takes a CSV with a group column and a mutations column and extracts the common mutations within the group, optionally with reference to a specified outgroup

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

## Valid Mutation Definitions
The following are valid ways to describe variants of each type. We prefer the definition at the top of each list, but provide alternatives for backwards compatibility. 
* these are case insensitive e.g. S vs s
* genes can be full e.g. orf1ab spike, or shortened e.g. 1ab, s
* protein based definitions may be acceptable if the reference JSON includes them but may not be shortened e.g. NSP2
* all coordinates are 1-based
* for amino acid mutations, reference can be longer than 1 amino acid

SNP:
* nuc:[`ref`]`nucleotide_coordinate`[`alt`]
* snp:[`ref`]`nucleotide_coordinate`[`alt`]

Amino acid mutation: 
* `gene`:[`ref`]`amino_acid_coordinate_relative_to_gene`[`alt`]
* `protein`:[`ref`]`amino_acid_coordinate_relative_to_protein`[`alt`]
* `gene`:[`ref`]`amino_acid_coordinate_relative_to_gene` - this allows any other aa to be called as alt
* aa:`gene`:[`ref`]`amino_acid_coordinate_relative_to_gene`[`alt`]
* aa:`protein`:[`ref`]`amino_acid_coordinate_relative_to_protein`[`alt`]
* aa:`gene`:[`ref`]`amino_acid_coordinate_relative_to_gene` - this allows any other aa to be called as alt

Deletion: 
* del:`nucleotide_coordinate`:`nucleotide_length`
* `gene`:[`ref`]`amino_acid_coordinate`-
* `gene`:[`ref`]`amino_acid_coordinate`del

Insertion (currently parsed but not typed):
* nuc:`nucleotide_coordinate`+`inserted_sequence` 
* snp:`nucleotide_coordinate`+`inserted_sequence` 
* `gene`:`amino_acid_coordinate_relative_to_gene`+`inserted_sequence` 
* aa:`gene`:`amino_acid_coordinate_relative_to_gene`+`inserted_sequence`
