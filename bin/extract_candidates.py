#!/usr/bin/env python3
"""Extract conservative CDS models from Miniprot GFF3; no orthology claims.

All models enter the TSV. Only models passing identity/coverage and structural
QC enter FASTA/GFF. GFF coordinates are 1-based inclusive; BED-style conversion
is used internally only. Genome access is disk-indexed, not loaded into RAM.
"""
import argparse
import csv
import json
from collections import Counter
from urllib.parse import unquote
from Bio import SeqIO
from Bio.Seq import Seq


def attrs(text):
    return {k: unquote(v) for k, v in
            (item.split('=', 1) for item in text.split(';') if '=' in item)}


def run(args):
    if not 0 < args.min_coverage <= 1 or not 0 <= args.min_identity <= 1:
        raise ValueError('Invalid identity or coverage threshold')
    refs = {}
    for rec in SeqIO.parse(args.references, 'fasta'):
        if rec.id in refs:
            raise ValueError(f'Duplicate reference ID: {rec.id}')
        seq = str(rec.seq).rstrip('*')
        if not seq:
            raise ValueError(f'Empty reference: {rec.id}')
        refs[rec.id] = len(seq)
    if not refs:
        raise ValueError('No reference proteins found')

    models, children = {}, []
    with open(args.gff) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError('Expected nine GFF columns')
            a = attrs(fields[8])
            if fields[2] == 'mRNA':
                mid = a['ID']
                if mid in models:
                    raise ValueError(f'Duplicate model ID: {mid}')
                target = a['Target'].split()
                if len(target) < 3 or target[0] not in refs:
                    raise ValueError(f'Unknown or malformed Target: {a["Target"]}')
                models[mid] = {'fields': fields, 'attrs': a, 'target': target,
                               'cds': [], 'children': []}
            elif fields[2] in ('CDS', 'stop_codon'):
                children.append((fields, a))
    for fields, a in children:
        for parent in a.get('Parent', '').split(','):
            if parent not in models:
                raise ValueError(f'Orphan GFF child: {parent}')
            models[parent]['children'].append(fields)
            if fields[2] == 'CDS':
                models[parent]['cds'].append(fields)

    genome = SeqIO.index(args.genome, 'fasta')
    columns = ['model_id', 'query_id', 'contig', 'start', 'end', 'strand',
               'query_start', 'query_end', 'query_length', 'query_coverage',
               'identity', 'alignment_score', 'cds_count', 'cds_length',
               'protein_length', 'status', 'qc_flags']
    counts = Counter()
    mapped = set()
    with open('candidates.tsv', 'w') as tab, \
         open('candidates.gff3', 'w') as gout, \
         open('candidates.cds.fna', 'w') as cout, \
         open('candidates.proteins.faa', 'w') as pout:
        writer = csv.DictWriter(tab, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        gout.write('##gff-version 3\n')
        for mid, model in models.items():
            f, a, target = model['fields'], model['attrs'], model['target']
            query, qs, qe = target[0], int(target[1]), int(target[2])
            if not 1 <= qs <= qe <= refs[query]:
                raise ValueError(f'Invalid query coordinates: {mid}')
            mapped.add(query)
            cov = (qe - qs + 1) / refs[query]
            identity = float(a.get('Identity', 'nan'))
            flags = []
            if cov < args.min_coverage:
                flags.append('low_query_coverage')
            if not identity >= args.min_identity:
                flags.append('low_or_missing_identity')
            for key in ('Frameshift', 'StopCodon'):
                # Check both transcript and CDS attributes, conservatively.
                values = [int(attrs(x[8]).get(key, 0)) for x in [f] + model['cds']]
                if any(v > 0 for v in values):
                    flags.append(key.lower())
            if qs > 1 or qe < refs[query]:
                flags.append('partial_query_alignment')
            if f[6] not in ('+', '-'):
                raise ValueError(f'Invalid strand: {mid}')
            cds = sorted(model['cds'], key=lambda x: int(x[3]), reverse=f[6] == '-')
            if not cds:
                flags.append('no_cds')
            intervals = sorted((int(x[3]), int(x[4])) for x in cds)
            if any(b[0] <= x[1] for x, b in zip(intervals, intervals[1:])):
                flags.append('overlapping_cds')
            sequence = ''
            for i, feature in enumerate(cds):
                if feature[0] != f[0] or feature[6] != f[6]:
                    raise ValueError(f'Inconsistent CDS location: {mid}')
                if feature[0] not in genome:
                    raise ValueError(f'Contig absent from genome: {feature[0]}')
                start, end = int(feature[3]), int(feature[4])
                if not 1 <= start <= end <= len(genome[feature[0]]):
                    raise ValueError(f'CDS outside contig: {mid}')
                phase = int(feature[7])
                if phase not in (0, 1, 2):
                    raise ValueError(f'Invalid phase: {mid}')
                fragment = genome[feature[0]].seq[start-1:end]
                if f[6] == '-':
                    fragment = fragment.reverse_complement()
                if i == 0:
                    # Only trim the first CDS phase. Later phases describe
                    # codons spanning junctions, not bases to delete.
                    fragment = fragment[phase:]
                    if phase:
                        flags.append('partial_initial_codon')
                elif phase != (3 - len(sequence) % 3) % 3:
                    flags.append('inconsistent_cds_phase')
                sequence += str(fragment).upper()
            if len(sequence) % 3:
                flags.append('incomplete_terminal_codon')
            usable = sequence[:len(sequence) // 3 * 3]
            protein = str(Seq(usable).translate(table=args.translation_table)) if usable else ''
            if '*' in protein:
                flags.append('stop_in_cds')
            if 'X' in protein:
                flags.append('ambiguous_translation')
            if not protein:
                flags.append('empty_protein')
            # Partial but structurally coherent models are retained for review.
            warnings = {'partial_query_alignment', 'partial_initial_codon'}
            flags = sorted(set(flags))
            passed = not (set(flags) - warnings)
            status = 'candidate' if passed else 'rejected'
            counts[status] += 1
            writer.writerow(dict(zip(columns, [mid, query, f[0], f[3], f[4], f[6],
                qs, qe, refs[query], f'{cov:.6f}', f'{identity:.6f}', f[5],
                len(cds), len(sequence), len(protein), status, ','.join(flags) or 'none'])))
            if passed:
                description = f'query={query} status=candidate'
                cout.write(f'>{mid} {description}\n{sequence}\n')
                pout.write(f'>{mid} {description}\n{protein}\n')
                for feature in [f] + model['children']:
                    gout.write('\t'.join(feature) + '\n')
    genome.close()
    with open('qc.json', 'w') as handle:
        json.dump({'reference_count': len(refs), 'model_count': len(models),
                   'candidate_count': counts['candidate'], 'rejected_count': counts['rejected'],
                   'queries_without_reported_alignment': sorted(set(refs) - mapped),
                   'note': 'No alignment is not evidence of biological absence. '
                           'Candidates are protein-guided models, not confirmed orthologues.'},
                  handle, indent=2)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--genome', required=True)
    p.add_argument('--references', required=True)
    p.add_argument('--gff', required=True)
    p.add_argument('--min-coverage', type=float, default=0.5)
    p.add_argument('--min-identity', type=float, default=0.3)
    p.add_argument('--translation-table', type=int, default=1)
    run(p.parse_args())
