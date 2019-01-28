import argparse
import json
import os
import subprocess

"""
Assumtions:
    - All four files (dna normal vaf, rna normal vaf, dna tumor vaf, and rna tumor vaf) all
        contain the same sites in the same order
    - The four files have the following format:
        <chrom>\t<pos>\t<ref>\t<depth>\t<ref_vaf>\t<minor_vaf>\t<a_vaf>\t<c_vaf>\t<g_vaf>\t\
<t_vaf>\t<n_vaf>
"""

parser = argparse.ArgumentParser()

# file inputs group
file_group = parser.add_argument_group('file_group')
file_group.add_argument('--dna-normal-vaf', type=str,
        help='dna normal vaf for sample')
file_group.add_argument('--dna-tumor-vaf', type=str,
        help='dna tumor vaf for sample')
file_group.add_argument('--rna-normal-vaf', type=str,
        help='rna normal vaf for sample')
file_group.add_argument('--rna-tumor-vaf', type=str,
        help='rna tumor vaf for sample')

parser.add_argument('--min-depth', type=int,
        default=10, help='Minimum depth threshold. Each of the 4 sites must be >= min depth.')
parser.add_argument('--max-dna-minor-vaf', type=float,
        default=.01, help='Maximum minor vaf in both dna normal and dna tumor for site to be \
considered for rna editing')
parser.add_argument('--min-rna-minor-vaf', type=float,
        default=.1, help='Minimum minor vaf in either rna normal or rna tumor for site to be \
considered for rna editing')
parser.add_argument('--table-output', type=str,
        default='output.tsv', help='name of .tsv table output')
parser.add_argument('--json-output', type=str,
        default='output.json', help='name of .json output')

parser.add_argument('--no-vaf-header', action='store_true',
        help='Indicates whether input vaf files have a header line')
# parser.add_argument('--input-int', type=int,
#         default=2, help='an input int')

args = parser.parse_args()

MAX_DNA_MINOR_VAF = args.max_dna_minor_vaf
MIN_RNA_MINOR_VAF = args.min_rna_minor_vaf
MIN_DEPTH = args.min_depth

def check_arguments():
    if args.dna_normal_vaf is None:
        raise ValueError('Must specify --dna-normal-vaf')
    if args.dna_tumor_vaf is None:
        raise ValueError('Must specify --dna-tumor-vaf')
    if args.rna_normal_vaf is None:
        raise ValueError('Must specify --rna-normal-vaf')
    if args.rna_tumor_vaf is None:
        raise ValueError('Must specify --rna-tumor-vaf')

def get_table_header():
    """Returns header for output table"""
    header = 'CHROM\tPOS\tREF\t'
    header += 'DNA_NORMAL_DEPTH\tDNA_NORMAL_REF_VAF\tDNA_NORMAL_MINOR_VAF\tDNA_NORMAL_A_VAF\t\
DNA_NORMAL_C_VAF\tDNA_NORMAL_G_VAF\tDNA_NORMAL_T_VAF\tDNA_NORMAL_N_VAF\t'
    header += 'DNA_TUMOR_DEPTH\tDNA_TUMOR_REF_VAF\tDNA_TUMOR_MINOR_VAF\tDNA_TUMOR_A_VAF\t\
DNA_TUMOR_C_VAF\tDNA_TUMOR_G_VAF\tDNA_TUMOR_T_VAF\tDNA_TUMOR_N_VAF\t'
    header += 'RNA_NORMAL_DEPTH\tRNA_NORMAL_REF_VAF\tRNA_NORMAL_MINOR_VAF\tRNA_NORMAL_A_VAF\t\
RNA_NORMAL_C_VAF\tRNA_NORMAL_G_VAF\tRNA_NORMAL_T_VAF\tRNA_NORMAL_N_VAF\t'
    header += 'RNA_TUMOR_DEPTH\tRNA_TUMOR_REF_VAF\tRNA_TUMOR_MINOR_VAF\tRNA_TUMOR_A_VAF\t\
RNA_TUMOR_C_VAF\tRNA_TUMOR_G_VAF\tRNA_TUMOR_T_VAF\tRNA_TUMOR_N_VAF\n'

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
    return pieces[0], pieces[1], int(pieces[3]), float(pieces[5])

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

def is_editing_site(dna_a_minor_vaf, dna_t_minor_vaf, rna_a_minor_vaf, rna_t_minor_vaf):
    """Returns true if thresholds are passed."""
    if dna_a_minor_vaf > MAX_DNA_MINOR_VAF or dna_t_minor_vaf > MAX_DNA_MINOR_VAF:
        return False

    if rna_a_minor_vaf < MIN_RNA_MINOR_VAF and rna_t_minor_vaf < MIN_RNA_MINOR_VAF:
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

def add_to_table(table_output, dna_a_dict, dna_t_dict, rna_a_dict, rna_t_dict):
    """Add position to table output.
    line format is
        - <chrom>\t<pos>\t<ref>\t
            <dna_normal_depth>\t<dna_normal_minor_vaf>\t<dna_normal_a_vaf>\t<dna_normal_c_vaf>\t<dna_normal_g_vaf>\t<dna_normal_g_vaf>\t<dna_normal_t_vaf>\t<dna_normal_n_vaf>\t
            <dna_tumor_depth>\t........\t
            <rna_normal_depth>\t........\t
            <rna_tumor_depth>\t........
    """
    chrom = dna_a_dict['chrom']
    pos = dna_a_dict['pos']
    ref = dna_a_dict['ref']

    line = ''
    line += f'{chrom}\t{pos}\t{ref}\t' 
    line += get_vaf_output_line(dna_a_dict) + '\t'
    line += get_vaf_output_line(dna_t_dict) + '\t'
    line += get_vaf_output_line(rna_a_dict) + '\t'
    line += get_vaf_output_line(rna_t_dict)

    table_output.append(line)

def add_to_json(json_output, dna_a_dict, dna_t_dict, rna_a_dict, rna_t_dict):
    """Add to json output"""
    d = {
            'chrom': dna_a_dict['chrom'],
            'pos': dna_a_dict['pos'],
            'ref': dna_a_dict['ref'],
    }

    d['dna_normal'] = get_vaf_output_dict(dna_a_dict)
    d['dna_tumor'] = get_vaf_output_dict(dna_t_dict)
    d['rna_normal'] = get_vaf_output_dict(rna_a_dict)
    d['rna_tumor'] = get_vaf_output_dict(rna_t_dict)

    json_output.append(d)

def write_editing_sites(dna_a_vaf_fp, dna_t_vaf_fp, rna_a_vaf_fp, rna_t_vaf_fp,
        no_input_header=False, table_output_fp='output.tsv', json_output_fp='output.json'):
    """Return editing site tups

    site tup: (chrom, pos)
    """


    table_output = []
    json_output = []

    positions = get_positions(dna_a_vaf_fp, no_input_header=no_input_header)
    positions = positions.intersection(get_positions(dna_t_vaf_fp, no_input_header=no_input_header))
    positions = positions.intersection(get_positions(rna_a_vaf_fp, no_input_header=no_input_header))
    positions = positions.intersection(get_positions(rna_t_vaf_fp, no_input_header=no_input_header))

    # {{chr1, 10): {'dna_a': line, 'dna_t': line, ..}}
    position_to_lines = {p:{} for p in positions}

    dna_a_f = open(dna_a_vaf_fp)
    dna_t_f = open(dna_t_vaf_fp)
    rna_a_f = open(rna_a_vaf_fp)
    rna_t_f = open(rna_t_vaf_fp)
    # kill header if needed
    if not no_input_header:
        dna_a_f.readline()
        dna_t_f.readline()
        rna_a_f.readline()
        rna_t_f.readline()

    for line in dna_a_f:
        chrom, pos = line.split('\t', 2)[:2]
        if (chrom, pos) in positions:
            position_to_lines[(chrom, pos)]['dna_a'] = line
    dna_a_f.close()
    for line in dna_t_f:
        chrom, pos = line.split('\t', 2)[:2]
        if (chrom, pos) in positions:
            position_to_lines[(chrom, pos)]['dna_t'] = line
    dna_t_f.close()
    for line in rna_a_f:
        chrom, pos = line.split('\t', 2)[:2]
        if (chrom, pos) in positions:
            position_to_lines[(chrom, pos)]['rna_a'] = line
    rna_a_f.close()
    for line in rna_t_f:
        chrom, pos = line.split('\t', 2)[:2]
        if (chrom, pos) in positions:
            position_to_lines[(chrom, pos)]['rna_t'] = line
    rna_t_f.close()

    for (chrom, pos), line_dict in position_to_lines.items():
        dna_a_line = line_dict['dna_a']
        dna_t_line = line_dict['dna_t']
        rna_a_line = line_dict['rna_a']
        rna_t_line = line_dict['rna_t']

        _, _, dna_a_depth, dna_a_minor_vaf = process_vaf_line_light(dna_a_line)
        _, _, dna_t_depth, dna_t_minor_vaf = process_vaf_line_light(dna_t_line)
        _, _, rna_a_depth, rna_a_minor_vaf = process_vaf_line_light(rna_a_line)
        _, _, rna_t_depth, rna_t_minor_vaf = process_vaf_line_light(rna_t_line)

    

        passes_depth = passes_depth_filter([dna_a_depth, dna_t_depth, rna_a_depth, rna_t_depth])
        valid_editing_site = is_editing_site(dna_a_minor_vaf, dna_t_minor_vaf,
                rna_a_minor_vaf, rna_t_minor_vaf)
        if passes_depth and valid_editing_site:
            dna_a_dict = process_vaf_line(dna_a_line)
            dna_t_dict = process_vaf_line(dna_t_line)
            rna_a_dict = process_vaf_line(rna_a_line)
            rna_t_dict = process_vaf_line(rna_t_line)

            add_to_table(table_output, dna_a_dict, dna_t_dict, rna_a_dict, rna_t_dict)
            add_to_json(json_output, dna_a_dict, dna_t_dict, rna_a_dict, rna_t_dict)

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
    write_editing_sites(args.dna_normal_vaf, args.dna_tumor_vaf, args.rna_normal_vaf,
            args.rna_tumor_vaf, no_input_header=args.no_vaf_header,
            table_output_fp=args.table_output, json_output_fp=args.json_output)


if __name__ == '__main__':
    main()
