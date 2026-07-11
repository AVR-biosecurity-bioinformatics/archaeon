# archaeon
`archaeon` projects curated reference mutations from archetype proteins onto annotated genomes and generates genome-specific annotation tracks for downstream variant analysis.

This pipeline is currently **IN EARLY DEVELOPMENT** and is not yet publicly usable.

# Scope
`archaeon` is a bioinformatics framework for projecting curated reference mutations (such as insecticide mutations) from archetype proteins onto annotated genomes.

`archaeon` translates biologically characterised mutations described in the literature into genome-specific coordinates within a target assembly. The resulting annotations can be used with existing variant analysis tools such as bcftools, bcftools csq, SnpEff, VEP, IGV, and JBrowse.

To achieve this, `archaeon` uses curated archetype proteins as a common coordinate framework and projects literature-derived mutations through sequence alignments onto these archetypes. These archetypes are then compared for orthology to a query reference assembly to locate the genomic coordinates of the known mutations, and a GFF file is output.

The framework of `archaeon` is intentially generic, and while it is currently developed for **insecticide resistance mutations**, this could be expanded in future to other type of resistance mutations or experimentally validated functional mutations.

# Background:
A major challenge in insecticide resistance surveillance is that resistance mutations are reported using protein coordinates from many different species. For example:
* Bactrocera dorsalis `α6 G275E`
* Drosophila melanogaster `Rdl A301S`
* Musca domestica `para L1014F`

These coordinates are often not directly transferable because:
* proteins differ in length
* orthologous proteins contain insertions/deletions
* paralogues may exist
* numbering schemes vary among studies
  
Consequently, genomic analyses frequently require manual curation and alignment before resistance loci can be screened. Archaeon aims to automates this process.

# Core Concept:
`archaeon` adopts an archetype-based coordinate system as proposed by [Mair et. al. 2016](https://scijournals.onlinelibrary.wiley.com/doi/full/10.1002/ps.4301). Rather than storing mutations relative to individual species, mutations are stored relative to curated archetype proteins representing a gene family. For example:
```
Literature:
B. dorsalis α6 G275E
D. melanogaster α6 G277E
↓
Archetype:
nAChR_alpha6 G275E
```

`archaeon` maintains a reference database of known insecticide mutations from a range of insect species mapped to archetypes.

When applied to a new genome, Archaeon:
1. Identifies the orthologous gene.
2. Aligns the target protein to the archetype.
3. Projects archetype mutation coordinates onto the target protein.
4. Converts projected residues to transcript and genomic codon coordinates.
5. Generates genome-specific annotation tracks.
This creates a stable and reproducible framework for transferring mutation knowledge across taxa.

# Supported resistance mutations and archetypes
The initial `archaeon` release will focus on a small set of well-characterised target-site resistance genes with established mutation nomenclature and clear orthology relationships. For these, we have aimed to select archetype sequences around which the field has already converged:

| Gene | Archetype Species | Archetype Protein | Rationale |
|--------|---------------------|-------------------|-----------|
| RyR | *Plutella xylostella* | Ryanodine receptor (RyR) | Canonical reference for diamide resistance mutations, including well-established coordinates such as G4946E and I4790M. |
| VGSC (*para*) | *Musca domestica* | Voltage-gated sodium channel (para) | Widely used kdr coordinate system for pyrethroid and DDT resistance (e.g. L1014F, M918T). |
| Ace1 | *Anopheles gambiae* | Acetylcholinesterase-1 (Ace-1) | Established coordinate system for organophosphate and carbamate resistance mutations such as G119S. |
| nAChR α6 | *Frankliniella occidentalis* | Nicotinic acetylcholine receptor α6 subunit | G275E has become the predominant coordinate system used across the spinosyn resistance literature and surveillance programs. |
| Rdl | *Drosophila melanogaster* | Resistance to dieldrin (Rdl) GABA receptor | Canonical coordinate system for cyclodiene and fipronil resistance mutations such as A301S and A301G. |


# Relationship to Existing Tools
`archaeon` is conceptually similar to [FRAST](https://www.frast.com.au/) in that it uses curated archetype proteins as a common coordinate framework and projects literature-derived mutations through sequence alignments. FRAST aligns protein sequences to curated archetypes and reports mutations relative to standardised archetype coordinates, helping resolve inconsistencies in species-specific residue numbering. 

`archaeon` instead focusses on generating genome annotation (GFF) files that an be used with existing variant analysis tools such as bcftools, bcftools csq, SnpEff, VEP, IGV, and JBrowse.

## References
* Mair, Wesley, et al. "Proposal for a unified nomenclature for target‐site mutations associated with resistance to fungicides." Pest management science 72.8 (2016): 1449-1459.
