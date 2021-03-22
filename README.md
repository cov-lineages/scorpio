# scorpio
serious constellations of reoccurring phylogenetically-independent origin

## Command line options:

1. `classify` - takes a set of lineage-defining constellations with rules and classifies sequences by them.
2. `haplotype` - takes a set of constellations and writes haplotypes (either as strings or individual columns).
3. `report` - creates a report HTML for a constellation

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
