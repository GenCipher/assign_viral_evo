import math

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    valid_nucleotides = {"A", "C", "G", "T"}
    mismatches = 0
    valid_length = 0


    ref = reference_sequence.upper()
    dist = distant_sequence.upper()


    for ref_nucleotide, dist_nucleotide in zip(ref, dist):
        # Skip gaps and ambiguous bases
        if ref_nucleotide not in valid_nucleotides or dist_nucleotide not in valid_nucleotides:
            continue


        valid_length += 1
        if ref_nucleotide != dist_nucleotide:
            mismatches += 1


    if valid_length == 0:
        return 0.0


    # Proportion of mismatches
    mm_prop = mismatches / valid_length


    # The JC formula is undefined once p reaches 0.75 (sequences saturated)
    if mm_prop >= 0.75:
        return float("inf")


    # Jukes-Cantor formula
    dist_jc = -0.75 * math.log(1.0 - (4.0 / 3.0) * mm_prop)


    # Scale back to mutation count
    return dist_jc * valid_length



def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    standard_nucleotides = {"A", "C", "G", "T"}
    purines = {"A", "G"}
    pyrimidines = {"C", "T"}


    transitions = 0
    transversions = 0
    valid_length = 0


    ref = reference_sequence.upper()
    dist = distant_sequence.upper()


    for r, d in zip(ref, dist):
        # Skip gaps and ambiguous bases
        if r not in standard_nucleotides or d not in standard_nucleotides:
            continue


        valid_length += 1
        if r != d:
            # Transition: purine <-> purine or pyrimidine <-> pyrimidine
            if (r in purines and d in purines) or (
                r in pyrimidines and d in pyrimidines
            ):
                transitions += 1
            else:
                transversions += 1


    if valid_length == 0:
        return 0.0


    # Frequencies of transitions (TS) and transversions (TV)
    TS = transitions / valid_length
    TV = transversions / valid_length


    term1 = 1.0 - (2.0 * TS) - TV
    term2 = 1.0 - (2.0 * TV)


    # The K2P formula is undefined if either log argument is <= 0 (saturation)
    if term1 <= 0 or term2 <= 0:
        return float("inf")


    # Kimura 2-Parameter formula
    dist_kimura2p = -0.5 * math.log(term1) - 0.25 * math.log(term2)


    # Scale back to mutation count
    return dist_kimura2p * valid_length

