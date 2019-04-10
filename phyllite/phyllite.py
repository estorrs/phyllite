import argparse
import json
import os
import subprocess

"""
Assumtions:
    - Two files (dna vaf, rna vaf) all
        contain the same sites in the same order
    - The files have the following format:
        <chrom>\t<pos>\t<ref>\t<depth>\t<ref_vaf>\t<minor_vaf>\t<a_vaf>\t<c_vaf>\t<g_vaf>\t\
<t_vaf>\t<n_vaf>
"""

parser = argparse.ArgumentParser()

# file inputs group
file_group = parser.add_argument_group('file_group')
file_group.add_argument('--dna-vaf', type=str,
        default=None, help='dna vaf for sample')
file_group.add_argument('--rna-vaf', type=str,
        default=None, help='rna vaf for sample')

parser.add_argument('--min-depth', type=int,
        default=10, help='Minimum depth threshold. Each of the sites must be >= min depth \
for both samples')
parser.add_argument('--max-dna-minor-count', type=int,
        default=1, help='Maximum number of reads supporting minor allele in dna for a site \
to be considered for rna editing.')
parser.add_argument('--max-dna-minor-vaf', type=float,
        default=.1, help='Maximum minor vaf in both dna normal and dna tumor for site to be \
considered for rna editing')
parser.add_argument('--min-rna-minor-count', type=int,
        default=3, help='Minimum number of reads supporting minor allele in rna for a site \
to be considered for rna editing.')
parser.add_argument('--min-rna-minor-vaf', type=float,
        default=.1, help='Minimum minor vaf in either rna normal or rna tumor for site to be \
considered for rna editing')
parser.add_argument('--table-output', type=str,
        default='output.tsv', help='name of .tsv table output')
parser.add_argument('--json-output', type=str,
        default='output.json', help='name of .json output')

parser.add_argument('--no-vaf-header', action='store_true',
        help='Indicates whether input vaf files have a header line')

args = parser.parse_args()

MAX_DNA_MINOR_VAF = args.max_dna_minor_vaf
MIN_RNA_MINOR_VAF = args.min_rna_minor_vaf
MAX_DNA_MINOR_COUNT = args.max_dna_minor_count
MIN_RNA_MINOR_COUNT = args.min_rna_minor_count
MIN_DEPTH = args.min_depth

def check_arguments():
    if args.dna_vaf is None and args.rna_vaf is None:
        raise ValueError('Must specify --dna-vaf and --rna-vaf')

def get_table_header():
    """Returns header for output table"""
    header = 'CHROM\tPOS\tREF\tALT\t'
    header += 'DNA_DEPTH\tDNA_REF_VAF\tDNA_MINOR_VAF\tDNA_A_VAF\tDNA_C_VAF\tDNA_G_VAF\tDNA_T_VAF\tDNA_N_VAF\t'
    header += 'RNA_DEPTH\tRNA_REF_VAF\tRNA_MINOR_VAF\tRNA_A_VAF\tRNA_C_VAF\tRNA_G_VAF\tRNA_T_VAF\tRNA_N_VAF\n'

    return header

def get_positions(fp, no_input_header=True):
    f = open(fp)
    if not no_input_header:
        f.readline()

    positions = set()
    for line in f:
        pieces = line.split('\t', 2)
        chrom = pieces[0]
        pos = pieces[1]
        positions.add((chrom, pos))

    return positions
        
def process_vaf_line_light(vaf_line):
    """Processes vaf line for key spots

    Returns (chrom, pos, depth, minor_vaf)
    """

    pieces = vaf_line.split('\t', 6)

    chrom = pieces[0]
    pos = pieces[1]
    depth = int(pieces[3])
    minor_vaf = float(pieces[5])
    minor_count = round(minor_vaf * depth)
    
    return chrom, pos, depth, minor_vaf, minor_count

def process_vaf_line(vaf_line):
    """Processes full vaf line into dict

    Returns {
        chrom: ...
        pos: ...
        ref: ...
        depth: ...
        ref_vaf: ...
        minor_vaf: ...
        a_vaf: ...
        c_vaf: ...
        g_vaf: ...
        t_vaf: ...
        n_vaf: ...
    }
    """
    pieces = vaf_line.strip().split('\t')

    d = {
            'chrom': pieces[0],
            'pos': pieces[1],
            'ref': pieces[2],
            'depth': int(pieces[3]),
            'ref_vaf': float(pieces[4]),
            'minor_vaf': float(pieces[5]),
            'a_vaf': float(pieces[6]),
            'c_vaf': float(pieces[7]),
            'g_vaf': float(pieces[8]),
            't_vaf': float(pieces[9]),
            'n_vaf': float(pieces[10])
            }

    return d

def is_editing_site_by_vaf(dna_minor_vaf, rna_minor_vaf):
    """Returns true if vaf thresholds are passed."""
    if dna_minor_vaf > MAX_DNA_MINOR_VAF:
        return False

    if rna_minor_vaf < MIN_RNA_MINOR_VAF:
        return False

    return True

def is_editing_site_by_count(dna_minor_count, rna_minor_count):
    """Returns true if count thresholds are passed."""
    if dna_minor_count > MAX_DNA_MINOR_COUNT:
        return False

    if rna_minor_count < MIN_RNA_MINOR_COUNT:
        return False

    return True

def passes_depth_filter(depths):
    for d in depths:
        if d < MIN_DEPTH:
            return False

    return True

def get_vaf_output_line(line_dict):
    order = ['depth', 'ref_vaf', 'minor_vaf', 'a_vaf', 'c_vaf', 'g_vaf', 't_vaf', 'n_vaf']

    pieces = []
    for k in order:
        pieces.append(str(line_dict[k]))

    return '\t'.join(pieces)

def get_vaf_output_dict(line_dict):
    keep = ['depth', 'ref_vaf', 'minor_vaf', 'a_vaf', 'c_vaf', 'g_vaf', 't_vaf', 'n_vaf']

    return {k:v for k, v in line_dict.items()
            if k in keep}

def get_base_change(ref_base, minor_vaf, a_vaf, c_vaf, g_vaf, t_vaf):
    vafs = [v for v in [a_vaf, c_vaf, g_vaf, t_vaf]
            if v <= minor_vaf]
    max_vaf = max(vafs)

    if 'a' != ref_base.lower() and a_vaf == max_vaf:
        return 'A', a_vaf
    if 'c' != ref_base.lower() and c_vaf == max_vaf:
        return 'C', c_vaf
    if 'g' != ref_base.lower() and g_vaf == max_vaf:
        return 'G', g_vaf
    if 't' != ref_base.lower() and t_vaf == max_vaf:
        return 'T', t_vaf

    return '.', 0.0

def add_to_table(table_output, dna_dict, rna_dict):
    """Add position to table output.
    line format is
        - <chrom>\t<pos>\t<ref>\t
            <dna_depth>\t<dna_minor_vaf>\t<dna_a_vaf>\t<dna_c_vaf>\t<dna_g_vaf>\t<dna_t_vaf>\t<dna_n_vaf>\t
            <rna_depth>\t........
    """
    chrom = dna_dict['chrom']
    pos = dna_dict['pos']
    ref = dna_dict['ref']
    alt, _ = get_base_change(ref, rna_dict['minor_vaf'], rna_dict['a_vaf'],
            rna_dict['c_vaf'], rna_dict['g_vaf'], rna_dict['t_vaf'])

    line = ''
    line += f'{chrom}\t{pos}\t{ref}\t{alt}\t' 
    line += get_vaf_output_line(dna_dict) + '\t'
    line += get_vaf_output_line(rna_dict)

    table_output.append(line)

def add_to_json(json_output, dna_dict, rna_dict):
    """Add to json output"""
    ref = dna_dict['ref']
    alt, _ = get_base_change(ref, rna_dict['minor_vaf'], rna_dict['a_vaf'],
            rna_dict['c_vaf'], rna_dict['g_vaf'], rna_dict['t_vaf'])

    d = {
            'chrom': dna_dict['chrom'],
            'pos': dna_dict['pos'],
            'ref': dna_dict['ref'],
            'alt': alt,
    }

    d['dna'] = get_vaf_output_dict(dna_dict)
    d['rna'] = get_vaf_output_dict(rna_dict)

    json_output.append(d)


def write_editing_sites(dna_vaf_fp, rna_vaf_fp,
        no_input_header=False, table_output_fp='output.tsv', json_output_fp='output.json'):
    """Return editing site tups

    site tup: (chrom, pos)
    """
    table_output = []
    json_output = []

#     fps_dict = {
#             'dna': dna_vaf_fp,
#             'rna': rna_vaf_fp,
#             }

    positions = set()
    for fp in [dna_vaf_fp, rna_vaf_fp]:
        if fp is not None:
            if len(positions) == 0:
                positions = get_positions(fp, no_input_header=no_input_header)
            else:
                positions = positions.intersection(get_positions(fp, no_input_header=no_input_header))

    position_to_lines = {p:{} for p in positions}
#     dna_fps_identifiers = [fp_id for fp_id, fp in fps_dict.items()
#             if fp is not None
#             if 'dna' in fp_id]
#     rna_fps_identifiers = [fp_id for fp_id, fp in fps_dict.items()
#             if fp is not None
#             if 'rna' in fp_id]

    for fp_identifier, fp in zip(['dna', 'rna'], [dna_vaf_fp, rna_vaf_fp]):
        if fp is not None:
            f = open(fp)
            # kill header if needed
            if not no_input_header:
                f.readline()

            for line in f:
                chrom, pos = line.split('\t', 2)[:2]
                if (chrom, pos) in positions:
                    position_to_lines[(chrom, pos)][fp_identifier] = line

                    if fp_identifier == 'dna':
                        position_to_lines[(chrom, pos)]['dna'] = line
                    if fp_identifier == 'rna':
                        position_to_lines[(chrom, pos)]['rna'] = line
            f.close()

    for (chrom, pos), line_dict in position_to_lines.items():
        dna_line = line_dict['dna']
        rna_line = line_dict['rna']

        _, _, dna_depth, dna_minor_vaf, dna_minor_count = process_vaf_line_light(dna_line)
        _, _, rna_depth, rna_minor_vaf, rna_minor_count = process_vaf_line_light(rna_line)
    
        passes_depth = passes_depth_filter([dna_depth, rna_depth])

        is_valid_editing_site_by_vaf = is_editing_site_by_vaf(dna_minor_vaf, rna_minor_vaf)
        is_valid_editing_site_by_count = is_editing_site_by_count(dna_minor_count, rna_minor_count)

        if passes_depth and is_valid_editing_site_by_vaf and is_valid_editing_site_by_count:
            dna_dict = process_vaf_line(dna_line)
            rna_dict = process_vaf_line(rna_line)

            add_to_table(table_output, dna_dict, rna_dict)
            add_to_json(json_output, dna_dict, rna_dict)

    # sort table output
    table_output = sorted(table_output)

    # write output table
    table_out = open(table_output_fp, 'w')
    table_out.write(get_table_header())
    table_out.write('\n'.join(table_output) + '\n')
    table_out.close()

    # write json output
    json.dump(json_output, open(json_output_fp, 'w'))

def main():
    check_arguments()

    print('writing editing sites')
    write_editing_sites(args.dna_vaf, args.rna_vaf, no_input_header=args.no_vaf_header,
            table_output_fp=args.table_output, json_output_fp=args.json_output)


if __name__ == '__main__':
    main()
