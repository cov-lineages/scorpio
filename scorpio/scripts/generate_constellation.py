import csv
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Generate constellation barcode from variant calls.')

    parser.add_argument("--variant-ref", action="store", type=str, dest="ref_csv")
    parser.add_argument("--variant-calls", action="store", type=str, dest="calls_csv")
    parser.add_argument("--out-file", action="store", type=str, dest="out_file")

    return parser.parse_args()

def record_to_variant(record):
    if record["type"] == 'deletion':
        return "del:%s:%s" %(record["nuc_location"],record["variants"])
    elif record["type"] == 'replacement':
        return "aa:%s" %record["id"]
    return None

def call_to_barcode(key, call, variant_dict):
    if call == 'del':
        return variant_dict[key]["variants"]
    elif call == 'ref' or call == variant_dict[key]["reference"]:
        return '-'
    else:
        return call

def generate_constellations(infile1, infile2, outfile):
    variant_dict = {}
    with open(infile1, 'r') as in_handle:
        reader = csv.DictReader(in_handle)
        for r in reader:
            variant = record_to_variant(r)
            variant_dict[variant] = r

    constellations = []
    with open(infile2, 'r') as in_handle:
        reader = csv.DictReader(in_handle)
        for read in reader:
            read_name = read['query']
            read_barcode = []
            for variant in variant_dict.keys():
                call = read[variant]
                barcode = call_to_barcode(variant, call, variant_dict)
                if barcode:
                    read_barcode.append(barcode)
            constellations.append({'sequence_name':read_name, 'constellation':''.join(read_barcode)})


    with open(outfile, 'w') as out_handle:
        f = csv.DictWriter(out_handle, fieldnames=['sequence_name','constellation'])
        f.writeheader()
        f.writerows(constellations)

if __name__ == '__main__':
    args = parse_args()
    generate_constellations(args.ref_csv, args.calls_csv, args.out_file)