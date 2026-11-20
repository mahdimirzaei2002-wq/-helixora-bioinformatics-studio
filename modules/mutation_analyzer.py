from Bio.Seq import Seq

AA_THREE_LETTER = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter"
}


def analyze_mutation(sequence, cds_start, cds_end, position,
                     reference, alternate, coordinate_system, accession):
    """Analyze one single-nucleotide substitution on a transcript sequence."""
    sequence = str(sequence).upper().replace("U", "T")
    reference = reference.upper()
    alternate = alternate.upper()

    if coordinate_system not in ("mRNA", "CDS"):
        raise ValueError("Coordinate system must be mRNA or CDS.")
    if reference not in "ACGT" or alternate not in "ACGT":
        raise ValueError("Reference and alternate must each be A, C, G, or T.")
    if reference == alternate:
        raise ValueError("Reference and alternate nucleotides must differ.")
    if not (1 <= cds_start <= cds_end <= len(sequence)):
        raise ValueError("CDS boundaries are invalid for this sequence.")

    cds_length = cds_end - cds_start + 1
    if cds_length % 3 != 0:
        raise ValueError(
            "The annotated CDS length is not divisible by 3. "
            "Check the selected transcript/CDS annotation."
        )

    if coordinate_system == "mRNA":
        mrna_position = int(position)
        cds_position = mrna_position - cds_start + 1
    else:
        cds_position = int(position)
        mrna_position = cds_start + cds_position - 1

    if not (1 <= cds_position <= cds_length):
        raise ValueError("The selected position is outside the CDS.")

    observed = sequence[mrna_position - 1]
    if observed != reference:
        raise ValueError(
            f"Reference mismatch at mRNA position {mrna_position}: "
            f"the reference sequence contains {observed}, not {reference}."
        )

    cds_sequence = sequence[cds_start - 1:cds_end]
    codon_number = (cds_position - 1) // 3 + 1
    codon_index = (codon_number - 1) * 3
    original_codon = cds_sequence[codon_index:codon_index + 3]
    offset = (cds_position - 1) % 3

    if len(original_codon) != 3:
        raise ValueError("Could not extract a complete codon.")

    mutated_codon = (
        original_codon[:offset] + alternate + original_codon[offset + 1:]
    )
    original_aa = str(Seq(original_codon).translate())
    mutated_aa = str(Seq(mutated_codon).translate())

    if original_aa == mutated_aa:
        mutation_type = "Synonymous"
        protein_change = f"p.{AA_THREE_LETTER[original_aa]}{codon_number}="
    elif original_aa == "*" and mutated_aa != "*":
        mutation_type = "Stop lost"
        protein_change = f"p.Ter{codon_number}{AA_THREE_LETTER[mutated_aa]}"
    elif mutated_aa == "*":
        mutation_type = "Stop gained"
        protein_change = f"p.{AA_THREE_LETTER[original_aa]}{codon_number}Ter"
    else:
        mutation_type = "Missense"
        protein_change = (
            f"p.{AA_THREE_LETTER[original_aa]}{codon_number}"
            f"{AA_THREE_LETTER[mutated_aa]}"
        )

    hgvs = f"{accession}:c.{cds_position}{reference}>{alternate}"

    return {
        "mRNA position": mrna_position,
        "CDS position": cds_position,
        "HGVS c.": f"c.{cds_position}{reference}>{alternate}",
        "HGVS (transcript)": hgvs,
        "Codon number": codon_number,
        "Reference": reference,
        "Alternate": alternate,
        "Original codon": original_codon,
        "Mutated codon": mutated_codon,
        "Original amino acid": original_aa,
        "Mutated amino acid": mutated_aa,
        "Protein change (predicted)": protein_change,
        "Consequence": mutation_type,
        "dbSNP / ClinVar query": hgvs
    }
