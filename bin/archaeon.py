#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

try:
    import bioframe
    import gffutils
    import pandas as pd
    import pysam
    from Bio import AlignIO, SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from pyfaidx import Fasta
except ImportError as exc:
    raise SystemExit(
        f"Missing Python dependency: {exc.name}. Install Archaeon dependencies "
        "with: mamba install -c conda-forge -c bioconda "
        "gffutils biopython bioframe pandas pyfaidx pysam"
    ) from exc


AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Stop": "*"
}


###############################################################################
# Generic utilities
###############################################################################

def open_text(path):
    """Open an uncompressed or gzip-compressed UTF-8 text file."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, "rt", encoding="utf-8-sig", errors="strict", newline="")


def require_columns(reader, required, source):
    """Report missing columns before processing records."""
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(
            f"{source} is missing required columns: "
            + ", ".join(sorted(missing))
        )


def require_executable(name):
    if shutil.which(name) is None:
        raise SystemExit(f"Required executable not found on PATH: {name}")


def run(command, stdout=None):
    try:
        subprocess.run(command, check=True, stdout=stdout)
    except subprocess.CalledProcessError as exc:
        raise SystemExit(
            "Command failed:\n  " + " ".join(map(str, command))
        ) from exc


def read_fasta(path):
    """Read FASTA records with Biopython and reject duplicate identifiers."""
    records = {}
    with open_text(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in records:
                raise SystemExit(f"Duplicate FASTA identifier in {path}: {record.id}")
            records[record.id] = str(record.seq)
    if not records:
        raise SystemExit(f"No FASTA records found in {path}")
    return records

def write_fasta(path, records):
    """Write sequence dictionaries as FASTA with Biopython."""
    seq_records = [
        SeqRecord(Seq(sequence), id=identifier, description="")
        for identifier, sequence in records.items()
    ]
    with open(path, "w", encoding="utf-8") as handle:
        SeqIO.write(seq_records, handle, "fasta")

def parse_gff_attributes(value):
    attributes = {}
    for item in value.rstrip().split(";"):
        if "=" in item:
            key, item_value = item.split("=", 1)
            attributes[key] = item_value
    return attributes


def clean_id(value):
    if value is None:
        return "."
    value = str(value).strip()
    if value in {"", "."}:
        return "."
    return re.sub(
        r"^(?:gene|transcript|mrna|mRNA|rna|cds|CDS|protein)[:-]",
        "",
        value
    )


def id_variants(value):
    value = clean_id(value)
    if value == ".":
        return set()
    variants = {value}
    if re.search(r"\.\d+$", value):
        variants.add(re.sub(r"\.\d+$", "", value))
    return variants


def normalise_aa_change(value):
    """
    Normalise a simple amino-acid substitution to REF-POS-ALT notation.

    Accepted examples:

        G252E       -> G252E
        p.G252E     -> G252E
        Gly252Glu   -> G252E
        252G>252E   -> G252E
        252G>E      -> G252E
        G252>E      -> G252E
        252G>253E   -> None

    bcftools csq commonly reports substitutions as 252G>252E.

    Synonymous changes such as 252G, frameshifts, indels and multi-residue
    changes return None.
    """

    if value is None:
        return None

    value = str(value).strip()

    if value in {"", ".", "NA"}:
        return None

    value = re.sub(r"^p\.", "", value)

    # Archaeon-style one-letter notation: G252E.
    match = re.fullmatch(
        r"([A-Z*])(\d+)([A-Z*])",
        value
    )

    if match:
        ref_aa, position, alt_aa = match.groups()
        return f"{ref_aa}{position}{alt_aa}"

    # Three-letter notation: Gly252Glu.
    match = re.fullmatch(
        r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|Stop)",
        value
    )

    if (
        match
        and match.group(1) in AA3
        and match.group(3) in AA3
    ):
        return (
            f"{AA3[match.group(1)]}"
            f"{match.group(2)}"
            f"{AA3[match.group(3)]}"
        )

    # bcftools csq notation: 252G>252E.
    match = re.fullmatch(
        r"(\d+)([A-Z*])>(\d+)([A-Z*])",
        value
    )

    if match:
        ref_pos, ref_aa, alt_pos, alt_aa = match.groups()

        # A simple substitution must affect the same amino-acid position.
        if ref_pos != alt_pos:
            return None

        return f"{ref_aa}{ref_pos}{alt_aa}"

    # Abbreviated coordinate-first notation: 252G>E.
    match = re.fullmatch(
        r"(\d+)([A-Z*])>([A-Z*])",
        value
    )

    if match:
        position, ref_aa, alt_aa = match.groups()
        return f"{ref_aa}{position}{alt_aa}"

    # Residue-first arrow notation: G252>E.
    match = re.fullmatch(
        r"([A-Z*])(\d+)>([A-Z*])",
        value
    )

    if match:
        ref_aa, position, alt_aa = match.groups()
        return f"{ref_aa}{position}{alt_aa}"

    return None

def vcf_string(value):
    if value is None:
        return "."
    value = str(value).strip()
    if value in {"", "."}:
        return "."
    return (
        value.replace(" ", "_")
        .replace("\t", "_")
        .replace("\n", "_")
        .replace(";", "&")
        .replace(",", "&")
        .replace("=", ":")
    )


def collapse_unique(values):
    unique = dict.fromkeys(
        value for value in map(vcf_string, values) if value != "."
    )
    return ",".join(unique) or "."


def parse_statuses(value):
    if value.lower() == "all":
        return None
    statuses = {item.strip() for item in value.split(",") if item.strip()}
    if not statuses:
        raise SystemExit("--include-status must contain a status or 'all'")
    return statuses


###############################################################################
# Identify target orthologues with DIAMOND
###############################################################################

def run_diamond(archetypes, proteome, output, db_prefix, threads, evalue,
                max_target_seqs):
    run([
        "diamond", "makedb", 
        "--in", str(proteome),
        "--db", str(db_prefix),
        "--quiet"
        ])
    run([
        "diamond", "blastp",
        "--query", str(archetypes),
        "--db", str(db_prefix),
        "--out", str(output),
        "--outfmt", "6",
        "qseqid", "sseqid", "pident", "length", "qlen", "slen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "qcovhsp", "scovhsp",
        "--evalue", str(evalue),
        "--max-target-seqs", str(max_target_seqs),
        "--more-sensitive",
        "--threads", str(threads),
        "--quiet"
    ])


def select_orthologues(hits, output, min_query_cover, min_target_cover):
    """Select the highest-bitscore passing target for each archetype."""
    fieldnames = [
        "archetype_id", "target_protein", "pident", "alignment_length",
        "archetype_length", "target_length", "qstart", "qend", "sstart",
        "send", "evalue", "bitscore", "query_coverage", "target_coverage"
    ]
    best = {}

    with open_text(hits) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for values in reader:
            if not values:
                continue
            row = dict(zip(fieldnames, values))
            if float(row["query_coverage"]) < min_query_cover:
                continue
            if float(row["target_coverage"]) < min_target_cover:
                continue
            query = row["archetype_id"]
            if query not in best or float(row["bitscore"]) > float(best[query]["bitscore"]):
                best[query] = row

    with open(output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames,
                                lineterminator="\n")
        writer.writeheader()
        for query in sorted(best):
            writer.writerow(best[query])

    return {query: row["target_protein"] for query, row in best.items()}


###############################################################################
# Align archetype and target proteins and build residue coordinate maps
###############################################################################

def parse_alignment(path):
    """Read a single MAFFT FASTA alignment with Biopython AlignIO."""
    alignment = AlignIO.read(str(path), "fasta")
    records = {record.id: str(record.seq) for record in alignment}
    expected = {"archetype", "target"}
    if set(records) != expected:
        raise SystemExit(
            f"Unexpected MAFFT records in {path}: expected {sorted(expected)}, "
            f"found {sorted(records)}"
        )
    return records["archetype"], records["target"]

def build_coordinate_map(archetype_id, target_protein, archetype_sequence,
                         target_sequence, workdir):
    pair = workdir / f"{archetype_id}.pair.faa"
    alignment = workdir / f"{archetype_id}.aln.faa"
    write_fasta(pair, {"archetype": archetype_sequence, "target": target_sequence})

    with open(alignment, "w", encoding="utf-8") as output:
        run(["mafft", "--quiet", "--auto", str(pair)], stdout=output)

    archetype_aln, target_aln = parse_alignment(alignment)
    coordinate_map = {}
    archetype_pos = 0
    target_pos = 0

    for archetype_aa, target_aa in zip(archetype_aln, target_aln):
        if archetype_aa != "-":
            archetype_pos += 1
        if target_aa != "-":
            target_pos += 1
        if archetype_aa != "-":
            coordinate_map[archetype_pos] = {
                "target_pos": None if target_aa == "-" else target_pos,
                "archetype_ref": archetype_aa,
                "target_ref": target_aa
            }

    return coordinate_map


###############################################################################
# Project known mutations into target-protein coordinates
###############################################################################
MUTATION_RE = re.compile(r"^([A-Z*])([0-9]+)([A-Z*])$")


def parse_mutation_row(row, line_number):
    """
    Extract a simple amino-acid substitution from a mutation database row.

    The explicit ref_aa, arch_pos and alt_aa columns are preferred.
    source_mutation is used to validate them or repair incomplete fields.
    """

    ref_aa = row["ref_aa"].strip()
    arch_pos = row["arch_pos"].strip()
    alt_aa = row["alt_aa"].strip()
    source_mutation = row.get("source_mutation", "").strip()

    explicit = MUTATION_RE.fullmatch(
        f"{ref_aa}{arch_pos}{alt_aa}"
    )
    source = MUTATION_RE.fullmatch(source_mutation)

    if explicit:
        mutation = explicit.groups()

        if source and source.groups() != mutation:
            print(
                f"[archaeon] warning: line {line_number}: "
                f"explicit mutation {''.join(mutation)} disagrees with "
                f"source_mutation {source_mutation}; using explicit fields",
                file=sys.stderr
            )

        return mutation

    if source:
        print(
            f"[archaeon] warning: line {line_number}: "
            f"using source_mutation={source_mutation} to repair "
            "ref_aa/arch_pos/alt_aa",
            file=sys.stderr
        )

        return source.groups()

    return None

def project_mutations(
    mutations,
    archetypes,
    proteome,
    output,
    workdir,
    threads,
    evalue,
    max_target_seqs,
    min_query_cover,
    min_target_cover
):
    required = {
        "gene_id",
        "archetype_id",
        "ref_aa",
        "arch_pos",
        "alt_aa",
        "source_mutation"
    }

    workdir.mkdir(parents=True, exist_ok=True)

    # Load archetype and target protein sequences.
    archetype_sequences = read_fasta(archetypes)
    target_sequences = read_fasta(proteome)

    hits = workdir / "diamond_hits.tsv"
    orthologue_table = workdir / "provisional_orthologues.tsv"

    # Identify candidate target proteins for each archetype.
    run_diamond(
        archetypes,
        proteome,
        hits,
        workdir / "target_proteome",
        threads,
        evalue,
        max_target_seqs
    )

    # Retain the highest-bitscore hit passing the coverage thresholds.
    orthologues = select_orthologues(
        hits,
        orthologue_table,
        min_query_cover,
        min_target_cover
    )

    # Align each archetype to its selected target protein and build
    # archetype-residue to target-residue coordinate maps.
    coordinate_maps = {}

    for archetype_id, target_protein in orthologues.items():
        if archetype_id not in archetype_sequences:
            raise SystemExit(
                f"Archetype absent from FASTA: {archetype_id}"
            )

        if target_protein not in target_sequences:
            raise SystemExit(
                f"Target protein absent from proteome: {target_protein}"
            )

        coordinate_maps[archetype_id] = build_coordinate_map(
            archetype_id,
            target_protein,
            archetype_sequences[archetype_id],
            target_sequences[target_protein],
            workdir
        )

    fieldnames = [
        "gene_id",
        "archetype_id",
        "ref_aa",
        "arch_pos",
        "alt_aa",
        "source_mutation",
        "target_pos",
        "target_ref",
        "target_mutation",
        "projection_status",
        "target_protein"
    ]

    stats = defaultdict(int)
    seen = set()

    with (
        open_text(mutations) as source,
        open(
            output,
            "w",
            encoding="utf-8",
            newline=""
        ) as destination
    ):
        reader = csv.DictReader(source, delimiter="\t")
        require_columns(reader, required, mutations)

        writer = csv.DictWriter(
            destination,
            delimiter="\t",
            fieldnames=fieldnames,
            lineterminator="\n"
        )
        writer.writeheader()

        for line_number, row in enumerate(reader, start=2):
            stats["mutation_rows"] += 1

            gene_id = row["gene_id"].strip()
            archetype_id = row["archetype_id"].strip()

            parsed = parse_mutation_row(row, line_number)

            if parsed is None:
                stats["invalid_source_mutation"] += 1

                print(
                    f"[archaeon] warning: line {line_number}: "
                    "could not parse a simple amino-acid substitution; "
                    "record skipped",
                    file=sys.stderr
                )

                continue

            ref_aa, arch_pos, alt_aa = parsed
            source_mutation = f"{ref_aa}{arch_pos}{alt_aa}"

            # Project each unique biological mutation only once. Multiple
            # evidence/reference rows for the same mutation do not require
            # separate protein alignments or projected database records.
            mutation_key = (
                gene_id,
                archetype_id,
                ref_aa,
                arch_pos,
                alt_aa
            )

            if mutation_key in seen:
                stats["duplicate_mutation_rows"] += 1
                continue

            seen.add(mutation_key)

            target_protein = orthologues.get(archetype_id, ".")
            target_pos = "."
            target_ref = "."
            target_mutation = "."

            if archetype_id not in orthologues:
                status = "orthologue_not_found"

            else:
                position = int(arch_pos)
                mapped = coordinate_maps[archetype_id].get(position)

                if mapped is None:
                    status = "position_not_in_coord_map"

                elif mapped["target_pos"] is None:
                    status = "target_gap"
                    target_ref = mapped["target_ref"]

                else:
                    target_pos = str(mapped["target_pos"])
                    target_ref = mapped["target_ref"]
                    target_mutation = (
                        f"{target_ref}{target_pos}{alt_aa}"
                    )

                    if ref_aa != mapped["archetype_ref"]:
                        status = "archetype_ref_mismatch"

                    elif ref_aa != target_ref:
                        status = "different_target_reference"

                    else:
                        status = "same_reference_residue"

            writer.writerow({
                "gene_id": gene_id,
                "archetype_id": archetype_id,
                "ref_aa": ref_aa,
                "arch_pos": arch_pos,
                "alt_aa": alt_aa,
                "source_mutation": source_mutation,
                "target_pos": target_pos,
                "target_ref": target_ref,
                "target_mutation": target_mutation,
                "projection_status": status,
                "target_protein": target_protein
            })

            stats["projected_output_rows"] += 1
            stats[status] += 1

    stats["orthologues"] = len(orthologues)

    return dict(stats)


###############################################################################
# Link projected proteins to identifiers used by split-vep
###############################################################################

DEFAULT_PROTEIN_ATTRIBUTES = (
    "Protein_Accession",
    "protein_accession",
    "protein_id",
    "protein",
    "Derives_from"
)


def open_gff_database(gff, database):
    """Create or open a gffutils SQLite database for hierarchical queries."""
    database = Path(database)
    database.parent.mkdir(parents=True, exist_ok=True)
    if not database.exists():
        gffutils.create_db(
            str(gff),
            dbfn=str(database),
            force=False,
            keep_order=True,
            merge_strategy="merge",
            sort_attribute_values=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True
        )
    return gffutils.FeatureDB(str(database), keep_order=True)


def feature_values(feature, keys):
    """Collect values for configured GFF attribute keys."""
    values = []
    for key in keys:
        values.extend(feature.attributes.get(key, []))
    return values


def feature_protein_ids(feature, protein_attributes):
    return {
        clean_id(value)
        for value in feature_values(feature, protein_attributes)
        if clean_id(value) != "."
    }


def feature_related_ids(db, feature):
    """Collect IDs and attributes from a CDS and its transcript/gene parents."""
    related = set(id_variants(feature.id))
    keys = ("ID", "Parent", "gene", "gene_id", "locus_tag", "Name", "accession")
    for value in feature_values(feature, keys):
        related.update(id_variants(value))
    try:
        parents = db.parents(feature, level=None)
    except (ValueError, gffutils.exceptions.FeatureNotFoundError):
        parents = []
    for parent in parents:
        related.update(id_variants(parent.id))
        for value in feature_values(parent, keys):
            related.update(id_variants(value))
    return related


def load_gff_identifiers(gff, database=None, protein_attributes=None):
    """Map protein accessions to CDS/transcript/gene IDs using gffutils."""
    protein_attributes = tuple(protein_attributes or DEFAULT_PROTEIN_ATTRIBUTES)
    database = database or Path(str(gff) + ".archaeon.sqlite")
    db = open_gff_database(gff, database)
    protein_ids = defaultdict(set)

    for cds in db.features_of_type("CDS", order_by=("seqid", "start")):
        proteins = feature_protein_ids(cds, protein_attributes)
        if not proteins:
            continue
        related = feature_related_ids(db, cds)
        for protein in proteins:
            protein_ids[protein].update(id_variants(protein))
            protein_ids[protein].update(related)
    return protein_ids

def load_resistance_database(projected, gff, include_statuses):
    required = {
        "gene_id",
        "archetype_id",
        "source_mutation",
        "target_mutation",
        "projection_status",
        "target_protein"
    }
    gff_ids = load_gff_identifiers(gff)
    database = defaultdict(list)
    stats = defaultdict(int)

    with open_text(projected) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader, required, projected)

        for row in reader:
            stats["projected_input_rows"] += 1
            if include_statuses is not None and row["projection_status"] not in include_statuses:
                stats["skipped_status"] += 1
                continue

            target_mutation = normalise_aa_change(row["target_mutation"])
            if target_mutation is None:
                stats["skipped_invalid_mutation"] += 1
                continue

            protein = clean_id(row["target_protein"])
            identifiers = set()
            identifiers.update(id_variants(row["gene_id"]))
            identifiers.update(id_variants(row["archetype_id"]))
            identifiers.update(id_variants(protein))
            identifiers.update(gff_ids.get(protein, set()))

            database[target_mutation].append({
                "gene_id": row["gene_id"],
                "archetype_id": row["archetype_id"],
                "target_protein": row["target_protein"],
                "source_mutation": row["source_mutation"],
                "target_mutation": target_mutation,
                "projection_status": row["projection_status"],
                "identifiers": identifiers
            })
            stats["database_rows"] += 1

    return database, dict(stats)


###############################################################################
# Match split-vep consequences against projected resistance mutations
###############################################################################

def consequence_identifiers(row):
    identifiers = set()
    for field in ("gene", "transcript"):
        for value in row[field].split(","):
            identifiers.update(id_variants(value))
    return identifiers


def match_consequences(
    consequences,
    database,
    annotations,
    report,
    strict_ids=True,
    annotate_unmatched=False
):
    """
    Match split-vep consequences against projected resistance mutations.

    When annotate_unmatched is enabled, every allele with a parsable simple
    amino-acid substitution is written to the annotation table, including
    alleles that do not match an Archaeon resistance mutation.
    """
    required = {
        "CHROM", "POS", "REF", "ALT", "Consequence", "gene",
        "transcript", "amino_acid_change", "dna_change"
    }
    protein_variants = defaultdict(list)
    matched_variants = defaultdict(list)
    report_rows = []
    stats = defaultdict(int)

    with open_text(consequences) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        # bcftools +split-vep normally writes #CHROM as the first column.
        if reader.fieldnames and reader.fieldnames[0] == "#CHROM":
            reader.fieldnames[0] = "CHROM"

        require_columns(reader, required, consequences)

        for consequence in reader:
            raw_aa_change = consequence["amino_acid_change"].strip()

            if raw_aa_change in {"", "."}:
                stats["missing_amino_acid_change"] += 1
                continue

            observed_mutation = normalise_aa_change(raw_aa_change)

            if observed_mutation is None:
                stats["non_substitution_consequences"] += 1
                continue

            stats["protein_affecting_records"] += 1

            try:
                position = int(consequence["POS"])
            except ValueError as exc:
                raise SystemExit(
                    f"{consequences}: invalid POS value "
                    f"{consequence['POS']!r}"
                ) from exc

            key = (
                consequence["CHROM"], position,
                consequence["REF"], consequence["ALT"]
            )
            protein_variants[key].append({
                "consequence": consequence,
                "observed_mutation": observed_mutation
            })

            consequence_ids = consequence_identifiers(consequence)
            for archaeon in database.get(observed_mutation, []):
                if strict_ids and not (
                    consequence_ids & archaeon["identifiers"]
                ):
                    stats["identifier_mismatches"] += 1
                    continue

                matched_variants[key].append({
                    "consequence": consequence,
                    "observed_mutation": observed_mutation,
                    "archaeon": archaeon
                })
                report_rows.append({
                    "chrom": consequence["CHROM"],
                    "pos": consequence["POS"],
                    "ref": consequence["REF"],
                    "alt": consequence["ALT"],
                    "archaeon_result": "known_resistance_match",
                    "consequence": consequence["Consequence"],
                    "gene": consequence["gene"],
                    "transcript": consequence["transcript"],
                    "observed_mutation": observed_mutation,
                    "dna_change": consequence["dna_change"],
                    "archaeon_gene": archaeon["gene_id"],
                    "archetype_id": archaeon["archetype_id"],
                    "target_protein": archaeon["target_protein"],
                    "source_mutation": archaeon["source_mutation"],
                    "target_mutation": archaeon["target_mutation"],
                    "projection_status": archaeon["projection_status"]
                })

    if annotate_unmatched:
        for key, records in protein_variants.items():
            if key in matched_variants:
                continue
            for record in records:
                consequence = record["consequence"]
                report_rows.append({
                    "chrom": consequence["CHROM"],
                    "pos": consequence["POS"],
                    "ref": consequence["REF"],
                    "alt": consequence["ALT"],
                    "archaeon_result": "protein_affecting_no_match",
                    "consequence": consequence["Consequence"],
                    "gene": consequence["gene"],
                    "transcript": consequence["transcript"],
                    "observed_mutation": record["observed_mutation"],
                    "dna_change": consequence["dna_change"],
                    "archaeon_gene": ".",
                    "archetype_id": ".",
                    "target_protein": ".",
                    "source_mutation": ".",
                    "target_mutation": ".",
                    "projection_status": "."
                })

    output_keys = (
        set(protein_variants) if annotate_unmatched else set(matched_variants)
    )

    with open(annotations, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for chrom, pos, ref, alt in sorted(output_keys):
            key = (chrom, pos, ref, alt)
            protein_records = protein_variants[key]
            matches = matched_variants.get(key, [])
            observed_mutations = collapse_unique(
                record["observed_mutation"] for record in protein_records
            )

            if matches:
                result = "known_resistance_match"
                archaeon_match = collapse_unique(
                    record["archaeon"]["target_mutation"] for record in matches
                )
                archaeon_gene = collapse_unique(
                    record["archaeon"]["gene_id"] for record in matches
                )
                source_mutation = collapse_unique(
                    record["archaeon"]["source_mutation"] for record in matches
                )
                target_protein = collapse_unique(
                    record["archaeon"]["target_protein"] for record in matches
                )
                projection_status = collapse_unique(
                    record["archaeon"]["projection_status"] for record in matches
                )
            else:
                result = "protein_affecting_no_match"
                archaeon_match = "."
                archaeon_gene = "."
                source_mutation = "."
                target_protein = "."
                projection_status = "."

            writer.writerow([
                chrom, pos, ref, alt,
                result,
                observed_mutations,
                archaeon_match,
                archaeon_gene,
                source_mutation,
                target_protein,
                projection_status
            ])

    report_fields = [
        "chrom", "pos", "ref", "alt", "archaeon_result", "consequence",
        "gene", "transcript", "observed_mutation", "dna_change",
        "archaeon_gene", "archetype_id", "target_protein",
        "source_mutation", "target_mutation", "projection_status"
    ]
    with open(report, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=report_fields,
            lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(report_rows)

    stats["matched_sites"] = len(matched_variants)
    stats["matched_records"] = sum(
        len(records) for records in matched_variants.values()
    )
    stats["protein_affecting_sites"] = len(protein_variants)
    stats["unmatched_protein_affecting_sites"] = len(
        set(protein_variants) - set(matched_variants)
    )
    stats["annotation_rows"] = len(output_keys)
    return dict(stats)


###############################################################################
# Write files consumed by the external bcftools annotate step
###############################################################################

def write_vcf_header(path):
    definitions = {
        "ARCHAEON_RESULT": (
            "Archaeon classification of the protein-affecting variant"
        ),
        "ARCHAEON_OBSERVED_MUTATION": (
            "Amino-acid substitution predicted by bcftools csq"
        ),
        "ARCHAEON_MATCH": (
            "Matched Archaeon resistance-associated amino-acid substitution"
        ),
        "ARCHAEON_GENE": (
            "Archaeon resistance-associated gene"
        ),
        "ARCHAEON_SOURCE_MUTATION": (
            "Resistance mutation in archetype protein coordinates"
        ),
        "ARCHAEON_PROTEIN": (
            "Target protein used for mutation projection"
        ),
        "ARCHAEON_STATUS": (
            "Archaeon mutation projection status"
        )
    }

    with open(path, "w", encoding="utf-8") as handle:
        handle.writelines(
            f'##INFO=<ID={tag},Number=.,Type=String,'
            f'Description="{description}">\n'
            for tag, description in definitions.items()
        )


def write_summary(path, *sections):
    with open(path, "w", encoding="utf-8") as handle:
        handle.writelines(
            f"{key}\t{value}\n"
            for section in sections
            for key, value in section.items()
        )


###############################################################################
# Build target-specific Archaeon bundle
###############################################################################

def projected_target_proteins(projected):
    """Return target proteins with a successfully projected mutation."""
    proteins = set()
    with open_text(projected) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(
            reader,
            {"target_protein", "target_mutation", "projection_status"},
            projected
        )
        for row in reader:
            if row["target_protein"] not in {"", "."} and row["target_mutation"] not in {"", "."}:
                proteins.add(clean_id(row["target_protein"]))
    return proteins


def merge_bed_intervals(intervals):
    """Merge overlapping or directly adjacent BED intervals per contig."""
    merged = []
    for chrom, start, end in sorted(intervals):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))
    return merged


def build_calling_bed(
    gff,
    target_proteins,
    output,
    metadata_output,
    padding=50,
    merge_intervals=True,
    database=None,
    protein_attributes=None,
    reference=None
):
    """Create validated, padded CDS BED intervals using gffutils and bioframe."""
    protein_attributes = tuple(protein_attributes or DEFAULT_PROTEIN_ATTRIBUTES)
    database = database or Path(str(gff) + ".archaeon.sqlite")
    db = open_gff_database(gff, database)
    records = []

    available_keys = set()
    matched_proteins = set()
    for cds in db.features_of_type("CDS", order_by=("seqid", "start")):
        available_keys.update(cds.attributes.keys())
        proteins = feature_protein_ids(cds, protein_attributes)
        matched = sorted(proteins & target_proteins)
        if not matched:
            continue
        matched_proteins.update(matched)
        for protein in matched:
            records.append({
                "chrom": cds.seqid,
                "start": max(0, int(cds.start) - 1 - padding),
                "end": int(cds.end) + padding,
                "protein": protein,
                "parent": ",".join(cds.attributes.get("Parent", ["."])),
                "strand": cds.strand,
                "feature": cds.featuretype
            })

    missing = sorted(target_proteins - matched_proteins)
    if missing:
        raise SystemExit(
            "Selected target proteins were not found on CDS features in the GFF: "
            + ", ".join(missing)
            + ". Protein attribute keys searched: "
            + ", ".join(protein_attributes)
            + ". Available CDS attribute keys: "
            + ", ".join(sorted(available_keys))
        )

    metadata = pd.DataFrame.from_records(records)
    if metadata.empty:
        raise SystemExit("No calling intervals were generated from the selected proteins")

    intervals = metadata[["chrom", "start", "end"]].drop_duplicates()

    # Clip padded intervals to reference contig bounds when a FASTA is provided.
    if reference:
        fasta = Fasta(str(reference), as_raw=True, sequence_always_upper=True)
        chromsizes = pd.Series(
            {
                record.name: len(record)
                for record in fasta
            },
            dtype="int64",
            name="length"
        )
        unknown = sorted(set(intervals["chrom"]) - set(chromsizes.index))
        if unknown:
            raise SystemExit(
                "GFF contigs absent from reference FASTA: " + ", ".join(unknown)
            )
        intervals["end"] = intervals.apply(
            lambda row: min(int(row["end"]), int(chromsizes[row["chrom"]])),
            axis=1
        )

    if merge_intervals:
        intervals = bioframe.merge(intervals, min_dist=0)
    intervals = bioframe.sort_bedframe(intervals)

    if not bioframe.is_bedframe(intervals):
        raise SystemExit("Generated calling intervals are not a valid BED frame")

    # NOTE: Bioframe.to_bed() applies a strict chromosome-name validator that rejects
    # accession versions containing dots, such as GWHFHGK00000001.1. These names
    # must be preserved because they match the reference FASTA, GFF and VCF.
    intervals[
        ["chrom", "start", "end"]
    ].to_csv(
        output,
        sep="\t",
        header=False,
        index=False,
        lineterminator="\n"
    )
    metadata.to_csv(metadata_output, sep="\t", index=False)

    return {
        "target_proteins": len(target_proteins),
        "calling_intervals": len(intervals),
        "calling_feature_records": len(metadata)
    }

def write_identifier_table(
    gff,
    target_proteins,
    output,
    database=None,
    protein_attributes=None
):
    """Persist protein-to-gene/transcript identifiers resolved by gffutils."""
    mappings = load_gff_identifiers(
        gff,
        database=database,
        protein_attributes=protein_attributes
    )
    with open(output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["target_protein", "identifier"])
        for protein in sorted(target_proteins):
            for identifier in sorted(mappings.get(protein, set())):
                writer.writerow([protein, identifier])

def load_identifier_table(path):
    mappings = defaultdict(set)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader, {"target_protein", "identifier"}, path)
        for row in reader:
            mappings[clean_id(row["target_protein"])].add(row["identifier"])
    return mappings


def load_built_resistance_database(projected, identifiers, include_statuses):
    """Load a built database without reparsing the reference GFF."""
    required = {
        "gene_id", "archetype_id", "source_mutation", "target_mutation",
        "projection_status", "target_protein"
    }
    identifier_map = load_identifier_table(identifiers)
    database = defaultdict(list)
    stats = defaultdict(int)

    with open_text(projected) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader, required, projected)
        for row in reader:
            stats["projected_input_rows"] += 1
            if include_statuses is not None and row["projection_status"] not in include_statuses:
                stats["skipped_status"] += 1
                continue
            target_mutation = normalise_aa_change(row["target_mutation"])
            if target_mutation is None:
                stats["skipped_invalid_mutation"] += 1
                continue

            protein = clean_id(row["target_protein"])
            ids = set()
            ids.update(id_variants(row["gene_id"]))
            ids.update(id_variants(row["archetype_id"]))
            ids.update(id_variants(protein))
            ids.update(identifier_map.get(protein, set()))

            database[target_mutation].append({
                "gene_id": row["gene_id"],
                "archetype_id": row["archetype_id"],
                "target_protein": row["target_protein"],
                "source_mutation": row["source_mutation"],
                "target_mutation": target_mutation,
                "projection_status": row["projection_status"],
                "identifiers": ids
            })
            stats["database_rows"] += 1

    return database, dict(stats)


###############################################################################
# Streaming VCF/BCF interface
###############################################################################

def materialise_input(bcftools, input_path, output_bcf):
    """Materialise stdin or a named VCF/BCF as an indexed temporary BCF."""

    input_path = str(input_path).strip() or "-"

    run([
        bcftools,
        "view",
        "-Ob",
        "-o", str(output_bcf),
        input_path
    ])

    run([
        bcftools,
        "index",
        "-f",
        str(output_bcf)
    ])


def extract_consequences(bcftools, input_bcf, output_tsv):
    """Validate and expand INFO/BCSQ into one row per consequence."""
    with pysam.VariantFile(str(input_bcf)) as variants:
        if "BCSQ" not in variants.header.info:
            raise SystemExit(
                "Input VCF/BCF does not contain INFO/BCSQ. Run bcftools csq first."
            )
    format_string = (
        "%CHROM\\t%POS\\t%REF\\t%ALT\\t%Consequence\\t%gene\\t"
        "%transcript\\t%amino_acid_change\\t%dna_change\\n"
    )

    with open(output_tsv, "wb") as handle:
        run([
            bcftools,
            "+split-vep",
            str(input_bcf),
            "-a", "BCSQ",
            "-d",
            "-HH",
            "-f", format_string
        ], stdout=handle)


def sort_annotation_table(input_path, output_path):
    """Sort the headerless Archaeon table for tabix indexing."""
    rows = []

    with open_text(input_path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if row:
                rows.append((row[0], int(row[1]), row[2], row[3], row))

    rows.sort(key=lambda item: item[:4])

    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for _, _, _, _, row in rows:
            writer.writerow(row)


def prepare_annotation_table(bgzip, tabix, input_path, output_path):
    """Bgzip-compress and tabix-index an Archaeon annotation table."""
    with open(output_path, "wb") as handle:
        run([bgzip, "-c", str(input_path)], stdout=handle)

    run([
        tabix,
        "-f",
        "-s", "1",
        "-b", "2",
        "-e", "2",
        str(output_path)
    ])


def write_annotated_output(
    bcftools,
    input_bcf,
    annotation_table,
    header,
    output,
    output_type,
    threads
):
    """Add Archaeon INFO fields while retaining upstream BCSQ annotations."""
    columns = (
        "CHROM,POS,REF,ALT,"
        "INFO/ARCHAEON_RESULT,"
        "INFO/ARCHAEON_OBSERVED_MUTATION,"
        "INFO/ARCHAEON_MATCH,"
        "INFO/ARCHAEON_GENE,"
        "INFO/ARCHAEON_SOURCE_MUTATION,"
        "INFO/ARCHAEON_PROTEIN,"
        "INFO/ARCHAEON_STATUS"
    )

    run([
        bcftools,
        "annotate",
        "-a", str(annotation_table),
        "-h", str(header),
        "-c", columns,
        "--threads", str(threads),
        f"-O{output_type}",
        "-o", str(output),
        str(input_bcf)
    ])


def write_header_only_output(
    bcftools,
    input_bcf,
    header,
    output,
    output_type,
    threads
):
    """Preserve BCSQ and add Archaeon header definitions when no rows match."""
    run([
        bcftools,
        "annotate",
        "-h", str(header),
        "--threads", str(threads),
        f"-O{output_type}",
        "-o", str(output),
        str(input_bcf)
    ])


###############################################################################
# Main workflow
###############################################################################

def add_output_options(parser):
    parser.add_argument("-i", "--input", default="-",
                        help="Input BCSQ-annotated VCF/BCF. Default: stdin.")
    parser.add_argument("-o", "--output", default="-",
                        help="Output Archaeon-annotated VCF/BCF. Default: stdout.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-Ob", dest="output_type", action="store_const", const="b")
    group.add_argument("-Ou", dest="output_type", action="store_const", const="u")
    group.add_argument("-Oz", dest="output_type", action="store_const", const="z")
    group.add_argument("-Ov", dest="output_type", action="store_const", const="v")
    parser.set_defaults(output_type="u")


def run_build(args):
    for executable in ("diamond", "mafft"):
        require_executable(executable)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    workdir = output_dir / "work"
    projected = output_dir / "projected_mutations.tsv"
    orthologues = output_dir / "orthologues.tsv"
    calling_bed = output_dir / "calling_regions.bed"
    calling_metadata = output_dir / "calling_regions.metadata.tsv"
    identifiers = output_dir / "target_identifiers.tsv"
    summary = output_dir / "build_summary.tsv"
    metadata = output_dir / "metadata.json"
    gff_database = output_dir / "features.sqlite"

    projection_stats = project_mutations(
        mutations=args.mutations,
        archetypes=args.archetypes,
        proteome=args.proteome,
        output=projected,
        workdir=workdir,
        threads=args.threads,
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
        min_query_cover=args.min_query_cover,
        min_target_cover=args.min_target_cover
    )

    source_orthologues = workdir / "provisional_orthologues.tsv"
    shutil.copyfile(source_orthologues, orthologues)

    targets = projected_target_proteins(projected)
    bed_stats = build_calling_bed(
        gff=args.gff,
        target_proteins=targets,
        output=calling_bed,
        metadata_output=calling_metadata,
        padding=args.padding,
        merge_intervals=not args.no_merge_intervals,
        database=gff_database,
        protein_attributes=args.protein_attribute,
        reference=args.reference
    )
    write_identifier_table(
        args.gff,
        targets,
        identifiers,
        database=gff_database,
        protein_attributes=args.protein_attribute
    )
    write_summary(summary, projection_stats, bed_stats)

    metadata.write_text(json.dumps({
        "format_version": 1,
        "reference": str(Path(args.reference).resolve()) if args.reference else None,
        "gff": str(Path(args.gff).resolve()),
        "proteome": str(Path(args.proteome).resolve()),
        "archetypes": str(Path(args.archetypes).resolve()),
        "mutations": str(Path(args.mutations).resolve()),
        "padding": args.padding,
        "merged_intervals": not args.no_merge_intervals
    }, indent=2) + "\n", encoding="utf-8")

    print(f"[archaeon build] database: {output_dir}", file=sys.stderr)
    print(f"[archaeon build] projected mutations: {projected}", file=sys.stderr)
    print(f"[archaeon build] calling BED: {calling_bed}", file=sys.stderr)
    print(f"[archaeon build] target proteins: {len(targets)}", file=sys.stderr)


def run_annotate(args):
    for executable in (args.bcftools, args.bgzip, args.tabix):
        require_executable(executable)

    database_dir = Path(args.database)
    projected = database_dir / "projected_mutations.tsv"
    identifiers = database_dir / "target_identifiers.tsv"
    for required in (projected, identifiers):
        if not required.is_file():
            raise SystemExit(f"Missing Archaeon database file: {required}")

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    annotations = Path(f"{prefix}.matches.tsv")
    report = Path(f"{prefix}.match_report.tsv")
    header = Path(f"{prefix}.hdr")
    summary = Path(f"{prefix}.summary.tsv")

    owner = None
    if args.work_dir:
        workdir = Path(args.work_dir)
        workdir.mkdir(parents=True, exist_ok=True)
    else:
        owner = tempfile.TemporaryDirectory(prefix="archaeon.")
        workdir = Path(owner.name)

    input_bcf = workdir / "input.csq.bcf"
    consequences = workdir / "consequences.tsv"
    sorted_annotations = workdir / "archaeon.matches.sorted.tsv"
    compressed_annotations = workdir / "archaeon.matches.sorted.tsv.gz"

    try:
        args.input = str(args.input).strip() or "-"
        args.output = str(args.output).strip() or "-"
        materialise_input(args.bcftools, args.input, input_bcf)
        extract_consequences(args.bcftools, input_bcf, consequences)

        resistance_db, database_stats = load_built_resistance_database(
            projected=projected,
            identifiers=identifiers,
            include_statuses=parse_statuses(args.include_status)
        )
        match_stats = match_consequences(
            consequences=consequences,
            database=resistance_db,
            annotations=annotations,
            report=report,
            strict_ids=not args.relaxed_id_match,
            annotate_unmatched=args.annotate_unmatched
        )
        write_vcf_header(header)
        write_summary(summary, database_stats, match_stats)

        if annotations.stat().st_size == 0:
            write_header_only_output(
                args.bcftools, input_bcf, header, args.output,
                args.output_type, args.threads
            )
        else:
            sort_annotation_table(annotations, sorted_annotations)
            prepare_annotation_table(
                args.bgzip, args.tabix, sorted_annotations,
                compressed_annotations
            )
            write_annotated_output(
                args.bcftools, input_bcf, compressed_annotations, header,
                args.output, args.output_type, args.threads
            )

        if args.output != "-" and args.output_type in {"b", "z"}:
            run([args.bcftools, "index", "-f", str(args.output)])

        print(f"[archaeon annotate] matched VCF alleles: {match_stats.get('matched_sites', 0)}",
              file=sys.stderr)
        print(f"[archaeon annotate] annotation table: {annotations}", file=sys.stderr)
        print(f"[archaeon annotate] report: {report}", file=sys.stderr)
    finally:
        if owner is not None:
            if args.keep_work_dir:
                retained = Path(f"{prefix}.work")
                if retained.exists():
                    shutil.rmtree(retained)
                shutil.copytree(workdir, retained)
                print(f"[archaeon annotate] retained work directory: {retained}",
                      file=sys.stderr)
            owner.cleanup()


def main():
    parser = argparse.ArgumentParser(
        prog="archaeon",
        description="Build target-specific resistance databases and annotate BCSQ VCF/BCF streams."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser(
        "build",
        help="Project mutations and build a reusable calling BED/database."
    )
    build.add_argument("--reference", default=None,
                       help="Reference FASTA recorded in database metadata.")
    build.add_argument("--gff", required=True)
    build.add_argument("--proteome", required=True)
    build.add_argument("--archetypes", required=True)
    build.add_argument("--mutations", required=True)
    build.add_argument(
        "--protein-attribute",
        action="append",
        default=None,
        help=(
            "CDS attribute containing target protein accessions. Repeat for "
            "multiple keys. Defaults: Protein_Accession, protein_accession, "
            "protein_id, protein, Derives_from."
        )
    )
    build.add_argument("--output-dir", required=True)
    build.add_argument("--threads", type=int, default=1)
    build.add_argument("--evalue", default="1e-10")
    build.add_argument("--max-target-seqs", type=int, default=20)
    build.add_argument("--min-query-cover", type=float, default=70.0)
    build.add_argument("--min-target-cover", type=float, default=70.0)
    build.add_argument("--padding", type=int, default=50,
                       help="Bases added around each CDS interval. Default: 50.")
    build.add_argument("--no-merge-intervals", action="store_true",
                       help="Do not merge overlapping or adjacent padded CDS intervals.")
    build.set_defaults(func=run_build)

    annotate = subparsers.add_parser(
        "annotate",
        help="Annotate a BCSQ VCF/BCF using a built Archaeon database."
    )
    add_output_options(annotate)
    annotate.add_argument("--database", required=True,
                          help="Directory created by `archaeon build`.")
    annotate.add_argument("--out-prefix", default="archaeon")
    annotate.add_argument("--work-dir", default=None)
    annotate.add_argument("--keep-work-dir", action="store_true")
    annotate.add_argument("--threads", type=int, default=1)
    annotate.add_argument("--include-status", default="same_reference_residue")
    annotate.add_argument("--relaxed-id-match", action="store_true")
    annotate.add_argument("--annotate-unmatched",
                          action="store_true")
    annotate.add_argument("--bcftools", default="bcftools")
    annotate.add_argument("--bgzip", default="bgzip")
    annotate.add_argument("--tabix", default="tabix")
    annotate.set_defaults(func=run_annotate)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
