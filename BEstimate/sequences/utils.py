

def find_pam_protospacer(
    sequence: str,
    pam_sequence: str,
    searched_nucleotide: str,
    activity_window: tuple,
    pam_window: tuple,
    protospacer_length: int,
) -> list:
    """
    Find all possible PAM and protospacer regions on the gene sequence.

    Searches for PAM (Protospacer Adjacent Motif) sequences and their corresponding
    protospacers within the given gene sequence. Identifies potential base editing
    sites within the activity window.

    :param sequence: The DNA sequence of the gene of interest
    :param pam_sequence: PAM sequence pattern (e.g., 'NGG', 'NG')
    :param searched_nucleotide: Target nucleotide to be edited by base editor
    :param activity_window: Tuple defining the activity window positions on protospacer (1-indexed)
    :param pam_window: Tuple defining PAM position relative to protospacer start (1-indexed)
    :param protospacer_length: Length of the protospacer sequence

    :return: List of dictionaries containing gRNA sequences, locations, and editing information

    .. note::
        The function uses regular expressions to find PAM patterns and considers
        both forward and reverse orientations. Positions are converted to 0-indexed
        for internal processing but returned as 1-indexed for biological relevance.
    """
    # Since python index starts from 0, decrease the start position index given from the user
    activity_window = [activity_window[0] - 1, activity_window[1]]
    pam_window = [pam_window[0] - 1, pam_window[1]]

    # Using Regular Expressions, specify PAM pattern
    pam_pattern = ""
    for nuc in list(pam_sequence):
        if nuc != "N":
            pam_pattern += nuc + "{1}"
        else:
            pam_pattern += "[ATCG]{1}"

    # Search protospacer length of nucleatide sequence and add PAM pattern after that
    pattern = r"[ATCG]{%s}%s" % (protospacer_length, pam_pattern)
    print("Pattern is created!")

    crisprs = []
    print("Pattern is searching through the sequence...")
    for nuc_index in range(0, len(sequence)):

        # One by one in the given sequence
        if nuc_index + pam_window[1] <= len(sequence):

            # Add started nucleotide index total length of targeted base editing site (PAM index)
            sub_sequence = sequence[nuc_index : nuc_index + pam_window[1]]

            # Search regex pattern inside the sub sequence
            for match_sequence in re.finditer(pattern, sub_sequence):

                # If there is match, then check if searched nucleotide inside the activity window
                if searched_nucleotide in list(
                    match_sequence.group()[activity_window[0] : activity_window[1]]
                ):
                    activity_sequence = match_sequence.group()[
                        activity_window[0] : activity_window[1]
                    ]

                    # If searched nucleotide is also there, add the sequence inside crisprs!
                    crisprs.append(
                        {
                            "index": [nuc_index, nuc_index + pam_window[1]],
                            "crispr": match_sequence.group(),
                            "activity_seq": activity_sequence,
                        }
                    )
                else:
                    crisprs.append(
                        {
                            "index": [nuc_index, nuc_index + pam_window[1]],
                            "crispr": match_sequence.group(),
                            "activity_seq": "No window",
                        }
                    )

    if crisprs is not []:
        print("CRISPRs are found!")
    return crisprs


def add_genomic_location(
    sequence_range: list, crispr_dict: dict, crispr_direction: str, strand: int
) -> tuple:
    """
    Adding genomic location info on crisprs found by extract_activity_window() function.

    :param sequence_range: The range of the sequence on the genome (from Ensembl)
    :param crispr_dict: The CRISPR dictionary created by extract_activity_window()
    :param crispr_direction: The direction of the created CRISPR ('left' or 'right')
    :param strand: The strand of the given gene (-1 or 1)
    :return crispr_seq: The sequence of the CRISPR

    :return genomic_location: The genomic coordinate of the above CRISPR on genome
    """
    # Prepare the sequence range according to the strand of the gene

    new_range = []
    if strand == 1:
        new_range = sequence_range
    elif strand == -1:
        new_range = [sequence_range[1], sequence_range[0]]

    # Look for both direction
    crispr_seq, genomic_location = crispr_dict["crispr"], ""

    if crispr_direction == "right":
        genomic_start = new_range[0] + crispr_dict["index"][0]
        genomic_end = (genomic_start + len(crispr_seq)) - 1
        genomic_location = str(genomic_start) + "-" + str(genomic_end)

    elif crispr_direction == "left":
        genomic_start = new_range[1] - crispr_dict["index"][0]
        genomic_end = (genomic_start - len(crispr_seq)) + 1
        genomic_location = str(genomic_end) + "-" + str(genomic_start)

    return crispr_seq, genomic_location


def extract_grna_sites(
    hugo_symbol: str,
    pam_sequence: str,
    searched_nucleotide: str,
    activity_window: list,
    pam_window: list,
    protospacer_length: int,
    flan: bool,
    flan_3: str,
    flan_5: str,
    ensembl_object: Ensembl,
) -> pandas.DataFrame:
    """
    Extract gRNA target sites containing editable nucleotides for base editing.

    Identifies all potential gRNA target sites within the gene sequence that contain
    the specified nucleotide within the activity window. Creates a comprehensive
    dataframe with gRNA sequences, locations, and relevant genomic annotations.

    :param hugo_symbol: HGNC gene symbol for the target gene
    :param pam_sequence: PAM sequence pattern (e.g., 'NGG', 'NG')
    :param searched_nucleotide: Nucleotide to be edited by base editor
    :param activity_window: Activity window positions on protospacer (1-indexed)
    :param pam_window: PAM position relative to protospacer start (1-indexed)
    :param protospacer_length: Length of the protospacer sequence
    :param flan: Whether to include flanking sequences for gRNAs
    :param flan_3: Number of nucleotides in 3' flanking region
    :param flan_5: Number of nucleotides in 5' flanking region
    :param ensembl_object: Ensembl object containing gene sequence and annotation data

    :return: DataFrame containing gRNA sites with annotations

    .. note::
        This function processes both forward and reverse strand orientations
        and includes transcript, exon, and CDS annotations for each gRNA site.
    """

    print("Sequence is preparing...")

    left_sequence = ensembl_object.flan_left_sequence_analysis
    right_sequence = ensembl_object.flan_right_sequence_analysis
    strand, seq_range = ensembl_object.strand, ensembl_object.flan_gene_range
    chromosome = ensembl_object.chromosome

    # Right CRISPRs 5'-->3' : reversed and base changed sequence

    # Editted should be used
    print("Protospacer and PAM regions are searching for right direction...")
    right_crisprs = find_pam_protospacer(
        sequence=right_sequence,
        pam_sequence=pam_sequence,
        searched_nucleotide=searched_nucleotide,
        activity_window=activity_window,
        pam_window=pam_window,
        protospacer_length=protospacer_length,
    )

    # Left CRISPRs: raw sequence will be used
    print("\nProtospacer and PAM regions are searching for left direction...")
    left_crisprs = find_pam_protospacer(
        sequence=left_sequence,
        pam_sequence=pam_sequence,
        searched_nucleotide=searched_nucleotide,
        activity_window=activity_window,
        pam_window=pam_window,
        protospacer_length=protospacer_length,
    )

    crisprs_dict = {"left": left_crisprs, "right": right_crisprs}

    print("\nCRISPR df is filling...")

    crisprs_df = pandas.DataFrame(
        columns=[
            "Hugo_Symbol",
            "CRISPR_PAM_Sequence",
            "gRNA_Target_Sequence",
            "Location",
            "Direction",
            "Gene_ID",
            "Transcript_ID",
            "Exon_ID",
            "guide_in_CDS",
        ]
    )

    for direction, crispr in crisprs_dict.items():

        for cr in crispr:
            crispr_seq, genomic_location = add_genomic_location(
                sequence_range=seq_range,
                strand=strand,
                crispr_dict=cr,
                crispr_direction=direction,
            )
            # Transcript & Exon Info
            transcript_exon = ensembl_object.check_range_info(
                int(genomic_location.split("-")[0]), int(genomic_location.split("-")[1])
            )
            if transcript_exon is not None:
                for transcript, exon_list in transcript_exon.items():
                    if exon_list:
                        for exon in exon_list:
                            df = pandas.DataFrame(
                                [
                                    [
                                        crispr_seq,
                                        crispr_seq[: -len(pam_sequence)],
                                        chromosome + ":" + genomic_location,
                                        direction,
                                        ensembl_object.gene_id,
                                        transcript,
                                        exon,
                                    ]
                                ],
                                columns=[
                                    "CRISPR_PAM_Sequence",
                                    "gRNA_Target_Sequence",
                                    "Location",
                                    "Direction",
                                    "Gene_ID",
                                    "Transcript_ID",
                                    "Exon_ID",
                                ],
                            )
                            crisprs_df = pandas.concat([crisprs_df, df])
                    else:
                        df = pandas.DataFrame(
                            [
                                [
                                    crispr_seq,
                                    crispr_seq[: -len(pam_sequence)],
                                    chromosome + ":" + genomic_location,
                                    direction,
                                    ensembl_object.gene_id,
                                    transcript,
                                    None,
                                ]
                            ],
                            columns=[
                                "CRISPR_PAM_Sequence",
                                "gRNA_Target_Sequence",
                                "Location",
                                "Direction",
                                "Gene_ID",
                                "Transcript_ID",
                                "Exon_ID",
                            ],
                        )
                        crisprs_df = pandas.concat([crisprs_df, df])
            else:
                df = pandas.DataFrame(
                    [
                        [
                            crispr_seq,
                            crispr_seq[: -len(pam_sequence)],
                            chromosome + ":" + genomic_location,
                            direction,
                            ensembl_object.gene_id,
                            None,
                            None,
                        ]
                    ],
                    columns=[
                        "CRISPR_PAM_Sequence",
                        "gRNA_Target_Sequence",
                        "Location",
                        "Direction",
                        "Gene_ID",
                        "Transcript_ID",
                        "Exon_ID",
                    ],
                )
                crisprs_df = pandas.concat([crisprs_df, df])

    crisprs_df["Hugo_Symbol"] = hugo_symbol
    crisprs_df["guide_in_CDS"] = crisprs_df.apply(
        lambda x: (
            ensembl_object.check_cds(
                x["Transcript_ID"],
                int(x["Location"].split(":")[1].split("-")[0]),
                int(x["Location"].split(":")[1].split("-")[1]) + 1,
            )
            if int(x["Location"].split(":")[1].split("-")[0])
            < int(x["Location"].split(":")[1].split("-")[1])
            else ensembl_object.check_cds(
                x["Transcript_ID"],
                int(x["Location"].split(":")[1].split("-")[1]),
                int(x["Location"].split(":")[1].split("-")[0]) + 1,
            )
        ),
        axis=1,
    )

    if flan:
        crisprs_df["gRNA_flanking_sequences"] = crisprs_df.apply(
            lambda x: ensembl_object.extract_gRNA_flan_sequence(
                location=x.Location, direction=x.Direction, fivep=flan_5, threep=flan_3
            ),
            axis=1,
        )
    else:
        crisprs_df["gRNA_flanking_sequences"] = None

    return crisprs_df


def collect_mutation_location(mutations: list | None) -> list | None:
    """
    TODO documentation
    """
    if mutations:
        locations = list()
        for mutation in mutations:
            alteration = mutation.split(":")[1].split(".")[1]
            mutation_location = re.match(
                "([0-9]+)([a-z]+)", alteration.split(">")[0], re.I
            ).groups()[0]
            if int(mutation_location) not in locations:
                locations.append(int(mutation_location))
        if locations:
            return locations
        else:
            return None
    else:
        return None


def check_genome_for_mutation(genomic_range, direction, mutations, window_type, window):
    """
    TODO documentation
    """
    yes_mutation = False

    if direction == "left":
        end = int(genomic_range.split("-")[0])
        start = int(genomic_range.split("-")[1])

    elif direction == "right":
        end = int(genomic_range.split("-")[1])
        start = int(genomic_range.split("-")[0])

    if window_type == "gRNA":
        if mutations:
            for loc in mutations:
                if (
                    int(genomic_range.split("-")[0])
                    <= loc
                    <= int(genomic_range.split("-")[1])
                ):
                    yes_mutation = True

    elif window_type == "activity":
        if direction == "right":
            act_start = start + window[0]
            act_end = start + window[1]
            activity_sites = list(range(act_start, act_end))
        elif direction == "left":
            act_start = start - window[0]
            act_end = start - window[1]
            activity_sites = list(range(act_end, act_start))
        if mutations:
            for loc in mutations:
                if loc in activity_sites:
                    yes_mutation = True

    elif window_type == "PAM":
        if direction == "right":
            pam_start = start + window[0]
            pam_end = start + window[1]
            pam_sites = list(range(pam_start, pam_end))
        elif direction == "left":
            pam_start = start - window[0]
            pam_end = start - window[1]
            pam_sites = list(range(pam_end, pam_start))
        if mutations:
            for loc in mutations:
                if loc in pam_sites:
                    yes_mutation = True

    return yes_mutation


def find_editable_nucleotide(
    crispr_df: pandas.DataFrame,
    searched_nucleotide: str,
    activity_window: list,
    pam_window: list,
    ensembl_object: Ensembl,
    mutations: list | None,
) -> pandas.DataFrame:
    """
    Find editable nucleotides within gRNA activity windows and map to genomic coordinates.

    Identifies specific nucleotides within the activity window of each gRNA that can
    be targeted for base editing. Maps these positions to genomic coordinates and
    incorporates mutation information if provided.

    :param crispr_df: DataFrame containing gRNA sequence, location and direction information
    :param searched_nucleotide: Target nucleotide to be edited by base editor
    :param activity_window: Activity window positions on protospacer (1-indexed)
    :param pam_window: PAM position relative to protospacer start (1-indexed)
    :param ensembl_object: Ensembl object containing gene sequence and annotation data
    :param mutations: List of user-provided genomic mutations to consider

    :return: DataFrame with editable nucleotide positions and genomic coordinates

    .. note::
        This function processes each gRNA from the input DataFrame and identifies
        all editable positions within the activity window, handling both strand
        orientations and mutation contexts.
    """

    actual_seq_range = ensembl_object.gene_range
    if actual_seq_range[0] > actual_seq_range[1]:
        actual_seq_range = [actual_seq_range[1], actual_seq_range[0]]

    actual_locations = list(range(actual_seq_range[0], actual_seq_range[1]))

    activity_window = [activity_window[0] - 1, activity_window[1]]
    pam_window = [pam_window[0] - 1, pam_window[1]]

    print("Edit Data Frame is filling...")
    edit_df = pandas.DataFrame(
        columns=[
            "Hugo_Symbol",
            "CRISPR_PAM_Sequence",
            "gRNA_Target_Sequence",
            "Location",
            "Edit_Location",
            "Direction",
            "Strand",
            "Gene_ID",
            "Transcript_ID",
            "Exon_ID",
            "guide_in_CDS",
            "gRNA_flanking_sequences",
            "Edit_in_Exon",
            "Edit_in_CDS",
            "GC%",
            "# Edits/guide",
            "Poly_T",
            "mutation_on_guide",
            "guide_change_mutation",
            "mutation_on_window",
            "mutation_on_PAM",
        ]
    )

    for ind, row in crispr_df.iterrows():
        # Check only with the sequence having PAM since it only has the searched nucleotide!
        try:
            searched_ind = [
                nuc_ind
                for nuc_ind in range(0, len(row["gRNA_Target_Sequence"]))
                if nuc_ind in list(range(activity_window[0], activity_window[1]))
                and row["gRNA_Target_Sequence"][nuc_ind] == searched_nucleotide
            ]
        except TypeError:
            print("TypeError on ", row)
        if searched_ind is not []:
            # If there is an editable nucleotide in the activity sites
            actual_inds = []
            if row["Direction"] == "left":
                for nuc_ind in searched_ind:
                    if (
                        int(row["Location"].split(":")[1].split("-")[1]) - nuc_ind
                        in actual_locations
                    ):
                        actual_inds.append(
                            int(row["Location"].split(":")[1].split("-")[1]) - nuc_ind
                        )

            elif row["Direction"] == "right":
                for nuc_ind in searched_ind:
                    if (
                        int(row["Location"].split(":")[1].split("-")[0]) + nuc_ind
                        in actual_locations
                    ):
                        actual_inds.append(
                            int(row["Location"].split(":")[1].split("-")[0]) + nuc_ind
                        )

            for actual_ind in actual_inds:

                transcript_exon = ensembl_object.check_range_info(
                    actual_ind, actual_ind + 1
                )

                if transcript_exon is not None:
                    for transcript, exon_list in transcript_exon.items():
                        if row["Exon_ID"] is not None and row["Exon_ID"] in exon_list:
                            edit_in_exon = True
                        else:
                            edit_in_exon = False
                else:
                    edit_in_exon = False

                edit_in_cds = ensembl_object.check_cds(
                    row["Transcript_ID"], actual_ind, actual_ind + 1
                )

                df = pandas.DataFrame(
                    [
                        [
                            row["Hugo_Symbol"],
                            row["CRISPR_PAM_Sequence"],
                            row["gRNA_Target_Sequence"],
                            row["Location"],
                            actual_ind,
                            row["Direction"],
                            ensembl_object.strand,
                            ensembl_object.gene_id,
                            row["Transcript_ID"],
                            row["Exon_ID"],
                            row["guide_in_CDS"],
                            row["gRNA_flanking_sequences"],
                            edit_in_exon,
                            edit_in_cds,
                        ]
                    ],
                    columns=[
                        "Hugo_Symbol",
                        "CRISPR_PAM_Sequence",
                        "gRNA_Target_Sequence",
                        "Location",
                        "Edit_Location",
                        "Direction",
                        "Strand",
                        "Gene_ID",
                        "Transcript_ID",
                        "Exon_ID",
                        "guide_in_CDS",
                        "gRNA_flanking_sequences",
                        "Edit_in_Exon",
                        "Edit_in_CDS",
                    ],
                )
                edit_df = pandas.concat([edit_df, df])

        else:
            # If not --> no edit
            df = pandas.DataFrame(
                [
                    [
                        row["Hugo_Symbol"],
                        row["CRISPR_PAM_Sequence"],
                        row["gRNA_Target_Sequence"],
                        row["Location"],
                        "No edit",
                        row["Direction"],
                        ensembl_object.strand,
                        ensembl_object.gene_id,
                        row["Transcript_ID"],
                        row["Exon_ID"],
                        row["guide_in_CDS"],
                        row["gRNA_flanking_sequences"],
                        False,
                        False,
                    ]
                ],
                columns=[
                    "Hugo_Symbol",
                    "CRISPR_PAM_Sequence",
                    "gRNA_Target_Sequence",
                    "Location",
                    "Edit_Location",
                    "Direction",
                    "Strand",
                    "Gene_ID",
                    "Transcript_ID",
                    "Exon_ID",
                    "guide_in_CDS",
                    "gRNA_flanking_sequences",
                    "Edit_in_Exon",
                    "Edit_in_CDS",
                ],
            )
            edit_df = pandas.concat([edit_df, df])

    edit_df["# Edits/guide"] = 0
    for guide, g_df in edit_df.groupby("gRNA_Target_Sequence"):
        unique_edits_per_guide = len(set(list(g_df["Edit_Location"])))
        edit_df.loc[edit_df.gRNA_Target_Sequence == guide, "# Edits/guide"] = (
            unique_edits_per_guide
        )

    edit_df["Poly_T"] = edit_df.apply(
        lambda x: (
            True if re.search("T{4,}", x.gRNA_Target_Sequence) is not None else False
        ),
        axis=1,
    )

    edit_df["GC%"] = edit_df.apply(
        lambda x: (
            x.gRNA_Target_Sequence.count("C") + x.gRNA_Target_Sequence.count("G")
        )
        * 100.0
        / len(x.gRNA_Target_Sequence),
        axis=1,
    )

    mutation_locations = collect_mutation_location(mutations=mutations)

    edit_df["mutation_on_guide"] = edit_df.apply(
        lambda x: check_genome_for_mutation(
            genomic_range=x.Location.split(":")[1],
            direction=x.Direction,
            mutations=mutation_locations,
            window_type="gRNA",
            window=None,
        ),
        axis=1,
    )
    edit_df["guide_change_mutation"] = edit_df.apply(
        lambda x: (
            True
            if mutation_locations is not None
            and int(x.Edit_Location) in mutation_locations
            else False
        ),
        axis=1,
    )

    edit_df["mutation_on_window"] = edit_df.apply(
        lambda x: check_genome_for_mutation(
            genomic_range=x.Location.split(":")[1],
            direction=x.Direction,
            mutations=mutation_locations,
            window_type="activity",
            window=activity_window,
        ),
        axis=1,
    )

    edit_df["mutation_on_PAM"] = edit_df.apply(
        lambda x: check_genome_for_mutation(
            genomic_range=x.Location.split(":")[1],
            direction=x.Direction,
            mutations=mutation_locations,
            window_type="PAM",
            window=pam_window,
        ),
        axis=1,
    )

    return edit_df


def extract_hgvs_df(
    edit_df: pandas.DataFrame,
    ensembl_object: Ensembl,
    transcript_id: str,
    edited_nucleotide: str,
    new_nucleotide: str,
    activity_window: list,
    mutations: list | None,
) -> pandas.DataFrame:
    """
    Generate HGVS nomenclature for base editing variants.

    Creates standardized HGVS (Human Genome Variation Society) nomenclature
    for all potential base editing sites identified in the analysis. Handles
    both genomic and transcript-level coordinate systems.

    :param edit_df: DataFrame containing editable nucleotide positions
    :param ensembl_object: Ensembl object with gene and transcript information
    :param transcript_id: Specific Ensembl transcript ID for analysis
    :param edited_nucleotide: Original nucleotide to be edited
    :param new_nucleotide: Target nucleotide after base editing
    :param activity_window: Activity window positions on protospacer (1-indexed)
    :param mutations: List of user-provided genomic mutations

    :return: DataFrame containing HGVS nomenclature for all variants

    .. note::
        This function handles strand orientation conversions and generates
        both genomic (g.) and coding sequence (c.) HGVS nomenclature as
        appropriate for VEP analysis.
    """
    # For (-) direction crisprs, base reversion should be done.
    nucleotide_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    # Collect chromosome
    chromosome, strand = ensembl_object.chromosome, ensembl_object.strand
    activity_window = [activity_window[0] - 1, activity_window[1]]

    # Collect mutations
    mutation_locations = collect_mutation_location(mutations)

    # Transcript filtration
    loc_edit_df = None
    if transcript_id is not None:
        loc_edit_df = edit_df[edit_df.Transcript_ID == transcript_id]
    else:
        for transcript, transcript_dict in ensembl_object.info_dict.items():
            loc_edit_df = edit_df[edit_df.Transcript_ID == transcript]

    # Each gRNA at a time
    row_dicts = list()
    for direction, direction_df in loc_edit_df.groupby("Direction"):
        if direction == "left" or "left" in list(direction):
            # Base reversion of the (-) direction crisprs
            rev_edited_nucleotide, rev_new_nucleotide = (
                nucleotide_dict[edited_nucleotide],
                nucleotide_dict[new_nucleotide],
            )

            for grna, grna_df in direction_df.groupby("gRNA_Target_Sequence"):

                total_edit = len(set(list(grna_df["Edit_Location"].values)))

                # For individual edits
                for edit_loc, grna_edit_df in grna_df.groupby("Edit_Location"):

                    if True not in grna_edit_df.mutation_on_window.unique():

                        hgvs = "%s:g.%s%s>%s" % (
                            str(chromosome),
                            str(edit_loc),
                            rev_edited_nucleotide,
                            rev_new_nucleotide,
                        )

                    elif (
                        len(list(grna_edit_df.mutation_on_window.unique())) == 1
                        and list(grna_edit_df.mutation_on_window.unique())[0] is None
                    ):

                        hgvs = "%s:g.%s%s>%s" % (
                            str(chromosome),
                            str(edit_loc),
                            rev_edited_nucleotide,
                            rev_new_nucleotide,
                        )

                    else:
                        start = (
                            int(
                                list(grna_df["Location"].values)[0]
                                .split(":")[1]
                                .split("-")[1]
                            )
                            - activity_window[1]
                            + 1
                        )
                        end = (
                            int(
                                list(grna_df["Location"].values)[0]
                                .split(":")[1]
                                .split("-")[1]
                            )
                            - activity_window[0]
                        )

                        mutations_on_window, mutation_edited = list(), False
                        for mutation in mutation_locations:
                            if start <= mutation <= end:
                                mutations_on_window.append(mutation)
                            if mutation == edit_loc:
                                mutation_edited = True

                        if mutation_edited:
                            for mut in mutations:
                                alteration = mut.split(".")[1]
                                if int(
                                    re.match(
                                        "([0-9]+)([a-z]+)",
                                        alteration.split(">")[0],
                                        re.I,
                                    ).groups()[0]
                                ) == int(edit_loc):
                                    ref_nuc = re.match(
                                        "([0-9]+)([a-z]+)",
                                        alteration.split(">")[1],
                                        re.I,
                                    ).groups()[0]
                                    hgvs = "%s:g.%s%s>%s" % (
                                        str(chromosome),
                                        str(edit_loc),
                                        ref_nuc,
                                        rev_new_nucleotide,
                                    )
                        else:
                            hgvs = "%s:g.%s%s>%s" % (
                                str(chromosome),
                                str(edit_loc),
                                rev_edited_nucleotide,
                                rev_new_nucleotide,
                            )

                    d = {
                        "Hugo_Symbol": list(grna_edit_df["Hugo_Symbol"].values)[0],
                        "Edit_Type": "individual",
                        "CRISPR_PAM_Sequence": grna_edit_df[
                            "CRISPR_PAM_Sequence"
                        ].values[0],
                        "CRISPR_PAM_Location": grna_edit_df["Location"].values[0],
                        "gRNA_Target_Sequence": grna,
                        "gRNA_Target_Location": grna_edit_df["Location"]
                        .values[0]
                        .split(":")[0]
                        + ":"
                        + str(
                            int(
                                grna_edit_df["Location"]
                                .values[0]
                                .split(":")[1]
                                .split("-")[0]
                            )
                            - 3
                        )
                        + "-"
                        + grna_edit_df["Location"]
                        .values[0]
                        .split(":")[1]
                        .split("-")[1],
                        "Total_Edit": total_edit,
                        "Edit_Location": edit_loc,
                        "Direction": direction,
                        "Transcript_ID": grna_edit_df["Transcript_ID"].values[0],
                        "Exon_ID": grna_edit_df["Exon_ID"].values[0],
                        "guide_in_CDS": grna_edit_df["guide_in_CDS"].values[0],
                        "gRNA_flanking_sequences": grna_edit_df[
                            "gRNA_flanking_sequences"
                        ].values[0],
                        "Edit_in_Exon": grna_edit_df["Edit_in_Exon"].values[0],
                        "Edit_in_CDS": grna_edit_df["Edit_in_CDS"].values[0],
                        "mutation_on_guide": grna_edit_df["mutation_on_guide"].values[
                            0
                        ],
                        "guide_change_mutation": grna_edit_df[
                            "guide_change_mutation"
                        ].values[0],
                        "mutation_on_window": grna_edit_df["mutation_on_window"].values[
                            0
                        ],
                        "mutation_on_PAM": grna_edit_df["mutation_on_PAM"].values[0],
                        "# Edits/guide": grna_edit_df["# Edits/guide"].values[0],
                        "Poly_T": grna_edit_df["Poly_T"].values[0],
                        "GC%": grna_edit_df["GC%"].values[0],
                        "HGVS": hgvs,
                    }
                    row_dicts.append(d)

                if total_edit > 1:
                    # For multiple edits
                    start = (
                        int(
                            list(grna_df["Location"].values)[0]
                            .split(":")[1]
                            .split("-")[1]
                        )
                        - activity_window[1]
                        + 1
                    )
                    end = (
                        int(
                            list(grna_df["Location"].values)[0]
                            .split(":")[1]
                            .split("-")[1]
                        )
                        - activity_window[0]
                    )
                    position = str(start) + "_" + str(end)

                    activity_sites = grna[activity_window[0] : activity_window[1]]
                    activity_sites = "".join(
                        [nucleotide_dict[n] for n in activity_sites[::-1]]
                    )
                    edited_activity_sites = activity_sites.replace(
                        rev_edited_nucleotide, rev_new_nucleotide
                    )
                    hgvs = "%s:g.%sdelins%s" % (
                        str(chromosome),
                        position,
                        edited_activity_sites,
                    )

                    d = {
                        "Hugo_Symbol": list(grna_edit_df["Hugo_Symbol"].values)[0],
                        "Edit_Type": "multiple",
                        "CRISPR_PAM_Sequence": grna_edit_df[
                            "CRISPR_PAM_Sequence"
                        ].values[0],
                        "CRISPR_PAM_Location": grna_edit_df["Location"].values[0],
                        "gRNA_Target_Sequence": grna,
                        "gRNA_Target_Location": grna_edit_df["Location"]
                        .values[0]
                        .split(":")[0]
                        + ":"
                        + str(
                            int(
                                grna_edit_df["Location"]
                                .values[0]
                                .split(":")[1]
                                .split("-")[0]
                            )
                            - 3
                        )
                        + "-"
                        + grna_edit_df["Location"]
                        .values[0]
                        .split(":")[1]
                        .split("-")[1],
                        "Total_Edit": total_edit,
                        "Edit_Location": position.split("_")[0]
                        + "-"
                        + position.split("_")[1],
                        "Direction": direction,
                        "Transcript_ID": grna_edit_df["Transcript_ID"].values[0],
                        "Exon_ID": grna_edit_df["Exon_ID"].values[0],
                        "guide_in_CDS": grna_edit_df["guide_in_CDS"].values[0],
                        "gRNA_flanking_sequences": grna_edit_df[
                            "gRNA_flanking_sequences"
                        ].values[0],
                        "Edit_in_Exon": grna_edit_df["Edit_in_Exon"].values[0],
                        "Edit_in_CDS": grna_edit_df["Edit_in_CDS"].values[0],
                        "mutation_on_guide": grna_edit_df["mutation_on_guide"].values[
                            0
                        ],
                        "guide_change_mutation": grna_edit_df[
                            "guide_change_mutation"
                        ].values[0],
                        "mutation_on_window": grna_edit_df["mutation_on_window"].values[
                            0
                        ],
                        "mutation_on_PAM": grna_edit_df["mutation_on_PAM"].values[0],
                        "# Edits/guide": grna_edit_df["# Edits/guide"].values[0],
                        "Poly_T": grna_edit_df["Poly_T"].values[0],
                        "GC%": grna_edit_df["GC%"].values[0],
                        "HGVS": hgvs,
                    }
                    row_dicts.append(d)

        elif direction == "right" or "right" in list(direction):

            for grna, grna_df in direction_df.groupby("gRNA_Target_Sequence"):

                total_edit = len(set(list(grna_df["Edit_Location"].values)))

                # For individual edits

                for edit_loc, grna_edit_df in grna_df.groupby("Edit_Location"):

                    if True not in grna_edit_df.mutation_on_window.unique():

                        hgvs = "%s:g.%s%s>%s" % (
                            str(chromosome),
                            str(edit_loc),
                            edited_nucleotide,
                            new_nucleotide,
                        )

                    elif (
                        len(list(grna_edit_df.mutation_on_window.unique())) == 1
                        and list(grna_edit_df.mutation_on_window.unique())[0] is False
                    ):

                        hgvs = "%s:g.%s%s>%s" % (
                            str(chromosome),
                            str(edit_loc),
                            edited_nucleotide,
                            new_nucleotide,
                        )

                    else:
                        end = (
                            int(
                                list(grna_df["Location"].values)[0]
                                .split(":")[1]
                                .split("-")[0]
                            )
                            + activity_window[1]
                            - 1
                        )
                        start = (
                            int(
                                list(grna_df["Location"].values)[0]
                                .split(":")[1]
                                .split("-")[0]
                            )
                            + activity_window[0]
                        )

                        mutations_on_window, mutation_edited = list(), False
                        for mutation in mutation_locations:
                            if start <= mutation <= end:
                                mutations_on_window.append(mutation)
                            if mutation == edit_loc:
                                mutation_edited = True

                        if mutation_edited:
                            for mut in mutations:
                                alteration = mut.split(".")[1]
                                if int(
                                    re.match(
                                        "([0-9]+)([a-z]+)",
                                        alteration.split(">")[0],
                                        re.I,
                                    ).groups()[0]
                                ) == int(edit_loc):
                                    ref_nuc = re.match(
                                        "([0-9]+)([a-z]+)",
                                        alteration.split(">")[1],
                                        re.I,
                                    ).groups()[0]
                                    hgvs = "%s:g.%s%s>%s" % (
                                        str(chromosome),
                                        str(edit_loc),
                                        ref_nuc,
                                        new_nucleotide,
                                    )
                        else:
                            hgvs = "%s:g.%s%s>%s" % (
                                str(chromosome),
                                str(edit_loc),
                                edited_nucleotide,
                                new_nucleotide,
                            )

                    d = {
                        "Hugo_Symbol": list(grna_edit_df["Hugo_Symbol"].values)[0],
                        "Edit_Type": "individual",
                        "CRISPR_PAM_Sequence": grna_edit_df[
                            "CRISPR_PAM_Sequence"
                        ].values[0],
                        "CRISPR_PAM_Location": grna_edit_df["Location"].values[0],
                        "gRNA_Target_Sequence": grna,
                        "gRNA_Target_Location": grna_edit_df["Location"]
                        .values[0]
                        .split(":")[0]
                        + ":"
                        + grna_edit_df["Location"].values[0].split(":")[1].split("-")[0]
                        + "-"
                        + str(
                            int(
                                grna_edit_df["Location"]
                                .values[0]
                                .split(":")[1]
                                .split("-")[1]
                            )
                            - 3
                        ),
                        "Total_Edit": total_edit,
                        "Edit_Location": edit_loc,
                        "Direction": direction,
                        "Transcript_ID": grna_edit_df["Transcript_ID"].values[0],
                        "Exon_ID": grna_edit_df["Exon_ID"].values[0],
                        "guide_in_CDS": grna_edit_df["guide_in_CDS"].values[0],
                        "gRNA_flanking_sequences": grna_edit_df[
                            "gRNA_flanking_sequences"
                        ].values[0],
                        "Edit_in_Exon": grna_edit_df["Edit_in_Exon"].values[0],
                        "Edit_in_CDS": grna_edit_df["Edit_in_CDS"].values[0],
                        "mutation_on_guide": grna_edit_df["mutation_on_guide"].values[
                            0
                        ],
                        "guide_change_mutation": grna_edit_df[
                            "guide_change_mutation"
                        ].values[0],
                        "mutation_on_window": grna_edit_df["mutation_on_window"].values[
                            0
                        ],
                        "mutation_on_PAM": grna_edit_df["mutation_on_PAM"].values[0],
                        "# Edits/guide": grna_edit_df["# Edits/guide"].values[0],
                        "Poly_T": grna_edit_df["Poly_T"].values[0],
                        "GC%": grna_edit_df["GC%"].values[0],
                        "HGVS": hgvs,
                    }

                    row_dicts.append(d)

                if total_edit > 1:
                    # For multiple edits
                    end = (
                        int(
                            list(grna_df["Location"].values)[0]
                            .split(":")[1]
                            .split("-")[0]
                        )
                        + activity_window[1]
                        - 1
                    )
                    start = (
                        int(
                            list(grna_df["Location"].values)[0]
                            .split(":")[1]
                            .split("-")[0]
                        )
                        + activity_window[0]
                    )
                    position = str(start) + "_" + str(end)

                    activity_sites = grna[activity_window[0] : activity_window[1]]
                    edited_activity_sites = activity_sites.replace(
                        edited_nucleotide, new_nucleotide
                    )
                    hgvs = "%s:g.%sdelins%s" % (
                        str(chromosome),
                        position,
                        edited_activity_sites,
                    )

                    d = {
                        "Hugo_Symbol": grna_edit_df["Hugo_Symbol"].values[0],
                        "Edit_Type": "multiple",
                        "CRISPR_PAM_Sequence": grna_edit_df[
                            "CRISPR_PAM_Sequence"
                        ].values[0],
                        "CRISPR_PAM_Location": grna_edit_df["Location"].values[0],
                        "gRNA_Target_Sequence": grna,
                        "gRNA_Target_Location": grna_edit_df["Location"]
                        .values[0]
                        .split(":")[0]
                        + ":"
                        + grna_edit_df["Location"].values[0].split(":")[1].split("-")[0]
                        + "-"
                        + str(
                            int(
                                grna_edit_df["Location"]
                                .values[0]
                                .split(":")[1]
                                .split("-")[1]
                            )
                            - 3
                        ),
                        "Total_Edit": total_edit,
                        "Edit_Location": position.split("_")[0]
                        + "-"
                        + position.split("_")[1],
                        "Direction": direction,
                        "Transcript_ID": grna_edit_df["Transcript_ID"].values[0],
                        "Exon_ID": grna_edit_df["Exon_ID"].values[0],
                        "guide_in_CDS": grna_edit_df["guide_in_CDS"].values[0],
                        "gRNA_flanking_sequences": grna_edit_df[
                            "gRNA_flanking_sequences"
                        ].values[0],
                        "Edit_in_Exon": grna_edit_df["Edit_in_Exon"].values[0],
                        "Edit_in_CDS": grna_edit_df["Edit_in_CDS"].values[0],
                        "mutation_on_guide": grna_edit_df["mutation_on_guide"].values[
                            0
                        ],
                        "guide_change_mutation": grna_edit_df[
                            "guide_change_mutation"
                        ].values[0],
                        "mutation_on_window": grna_edit_df["mutation_on_window"].values[
                            0
                        ],
                        "mutation_on_PAM": grna_edit_df["mutation_on_PAM"].values[0],
                        "# Edits/guide": grna_edit_df["# Edits/guide"].values[0],
                        "Poly_T": grna_edit_df["Poly_T"].values[0],
                        "GC%": grna_edit_df["GC%"].values[0],
                        "HGVS": hgvs,
                    }
                    row_dicts.append(d)

    hgvs_df = pandas.DataFrame(row_dicts)
    return hgvs_df



def retrieve_vep_info(
    hgvs_df: pandas.DataFrame,
    ensembl_object: Ensembl,
    uniprot: str | None,
    transcript_id: str | None = None,
) -> pandas.DataFrame:
    """
    Retrieve variant effect predictions using Ensembl VEP API.

    Collects comprehensive variant annotation data from the Ensembl Variant Effect
    Predictor (VEP) for all base edits identified in the analysis. Includes
    consequence predictions, protein effects, and clinical significance.

    :param hgvs_df: DataFrame containing HGVS nomenclature for all variants
    :param ensembl_object: Ensembl object containing gene and transcript information
    :param uniprot: User-specified UniProt accession ID
    :param transcript_id: Specific Ensembl transcript ID to analyze, default None

    :return: DataFrame enriched with VEP annotation data

    .. note::
        This function makes batch requests to the VEP API for efficient processing.
        It handles rate limiting and includes comprehensive variant consequence
        predictions, pathogenicity scores, and clinical annotations.
    """

    chromosome, strand = ensembl_object.chromosome, ensembl_object.strand

    if transcript_id is None:
        transcript_id = list(ensembl_object.info_dict.keys())[0]

    vep_columns = [
        "Protein_ID",
        "VEP_input",
        "allele",
        "variant_classification",
        "most_severe_consequence",
        "consequence_terms",
        "variant_biotype",
        "Regulatory_ID",
        "Motif_ID",
        "TFs_on_motif",
        "cDNA_Change",
        "Edited_Codon",
        "New_Codon",
        "CDS_Position",
        "Protein_Position_ensembl",
        "Protein_Change",
        "Edited_AA",
        "Edited_AA_Prop",
        "New_AA",
        "New_AA_Prop",
        "is_Synonymous",
        "is_Stop",
        "proline_addition",
        "swissprot_vep",
        "uniprot_provided",
        "polyphen_score",
        "polyphen_prediction",
        "sift_score",
        "sift_prediction",
        "cadd_phred",
        "cadd_raw",
        "lof",
        "impact",
        "blosum62",
        "is_clinical",
        "clinical_id",
        "clinical_significance",
        "cosmic_id",
        "clinvar_id",
        "ancestral_populations",
    ]

    for c in vep_columns:
        hgvs_df[c] = None

    print("VEP Data Frame is filling with VEP API.")
    # Decide the server
    server = (
        "http://grch37.rest.ensembl.org"
        if ensembl_object.assembly == "hg19"
        else "https://rest.ensembl.org"
    )
    ext = "/vep/human/hgvs"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    params = {
        "AncestralAllele": 1,
        "Blosum62": 1,
        "Conservation": 1,
        "LoF": 1,
        "CADD": 1,
        "protein": 1,
        "variant_class": 1,
        "hgvs": 1,
        "uniprot": 1,
        "transcript_id": transcript_id,
    }

    hgvs_index = pandas.DataFrame(
        columns=["HGVS"], index=list(range(len(list(hgvs_df["HGVS"].unique()))))
    )
    count = 0
    for hgvs in list(hgvs_df["HGVS"].unique()):
        hgvs_index.loc[count, "HGVS"] = hgvs
        count += 1

    t = count // 200
    r = count - (200 * t)
    hgvs_obj = dict()
    for i in range(t):
        x = 200 * i
        hgvs_list = list(hgvs_index.loc[x : x + 199]["HGVS"].values)
        hgvs_json = json.dumps(hgvs_list)

        check_point = 0
        max_retry = 3
        while check_point < max_retry:
            try:
                vep_request = requests.post(
                    server + ext,
                    headers=headers,
                    params=params,
                    data='{ "hgvs_notations" : %s }' % hgvs_json,
                )

                if vep_request.status_code == requests.codes.ok:
                    break
                else:
                    print("No response in %s" % x)
            except requests.exceptions.RequestException as e:
                # time.sleep(2 ** check_point)
                check_point += 1
                if check_point >= max_retry:
                    raise e

        if vep_request.status_code != requests.codes.ok:
            print("No response from VEP %d - %d" % (x, x + 199))
        else:
            try:
                whole_vep = json.loads(vep_request.text)
                for hgvs in hgvs_list:
                    obj = Variant(
                        hgvs=hgvs,
                        gene=ensembl_object.hugo_symbol,
                        transcript=transcript_id,
                        strand=strand,
                    )
                    obj.extract_vep_obj(vep_json=whole_vep)
                    obj.extract_consequences()
                    hgvs_obj[hgvs] = obj
            except json.decoder.JSONDecodeError:
                print("No retrieval for %s" % hgvs)

    hgvs_list = list(
        hgvs_index.loc[(200 * (t - 1)) + 200 : (200 * (t - 1)) + 200 + r]["HGVS"].values
    )
    hgvs_json = json.dumps(hgvs_list)

    check_point = 0
    max_retry = 3
    while check_point < max_retry:
        try:
            vep_request = requests.post(
                server + ext,
                headers=headers,
                params=params,
                data='{ "hgvs_notations" : %s }' % hgvs_json,
            )
            if vep_request.status_code == requests.codes.ok:
                break

        except requests.exceptions.RequestException as e:
            # time.sleep(2 ** check_point)
            check_point += 1
            if check_point >= max_retry:
                raise e

    if vep_request.status_code != requests.codes.ok:
        print("No response from VEP %d - %d" % (x + 200, x + 200 + r))
    else:
        whole_vep = vep_request.json()
        for hgvs in hgvs_list:
            obj = Variant(
                hgvs=hgvs,
                gene=ensembl_object.hugo_symbol,
                transcript=transcript_id,
                strand=strand,
            )
            obj.extract_vep_obj(vep_json=whole_vep)
            obj.extract_consequences()
            hgvs_obj[hgvs] = obj

    for hgvs, obj in hgvs_obj.items():
        ind = list(hgvs_df[hgvs_df.HGVS == hgvs].index)
        hgvs_df.loc[ind, "Protein_ID"] = obj.protein
        hgvs_df.loc[ind, "VEP_input"] = obj.hgvs
        hgvs_df.loc[ind, "allele"] = obj.allele
        hgvs_df.loc[ind, "variant_classification"] = obj.variant_class
        hgvs_df.loc[ind, "most_severe_consequence"] = obj.most_severe_consequence
        hgvs_df.loc[ind, "consequence_terms"] = obj.consequence_terms
        hgvs_df.loc[ind, "variant_biotype"] = obj.biotype
        hgvs_df.loc[ind, "Regulatory_ID"] = obj.regulatory
        hgvs_df.loc[ind, "Motif_ID"] = obj.motif
        hgvs_df.loc[ind, "TFs_on_motif"] = obj.motif_TFs
        hgvs_df.loc[ind, "cDNA_Change"] = obj.cdna_change
        hgvs_df.loc[ind, "Edited_Codon"] = obj.old_codon
        hgvs_df.loc[ind, "New_Codon"] = obj.new_codon
        hgvs_df.loc[ind, "CDS_Position"] = obj.cds_position
        hgvs_df.loc[ind, "Protein_Position_ensembl"] = obj.protein_position
        hgvs_df.loc[ind, "Protein_Change"] = obj.protein_change
        hgvs_df.loc[ind, "Edited_AA"] = obj.old_aa
        hgvs_df.loc[ind, "Edited_AA_Prop"] = obj.old_aa_chem
        hgvs_df.loc[ind, "New_AA"] = obj.new_aa
        hgvs_df.loc[ind, "New_AA_Prop"] = obj.new_aa_chem
        hgvs_df.loc[ind, "is_Synonymous"] = obj.synonymous
        hgvs_df.loc[ind, "is_Stop"] = obj.stop
        hgvs_df.loc[ind, "proline_addition"] = obj.proline
        hgvs_df.loc[ind, "swissprot_vep"] = obj.swissprot
        hgvs_df.loc[ind, "uniprot_provided"] = uniprot
        hgvs_df.loc[ind, "polyphen_score"] = obj.polyphen_score
        hgvs_df.loc[ind, "polyphen_prediction"] = obj.polyphen_prediction
        hgvs_df.loc[ind, "sift_score"] = obj.sift_score
        hgvs_df.loc[ind, "sift_prediction"] = obj.sift_prediction
        hgvs_df.loc[ind, "cadd_phred"] = obj.cadd_phred
        hgvs_df.loc[ind, "cadd_raw"] = obj.cadd_raw
        hgvs_df.loc[ind, "lof"] = obj.lof
        hgvs_df.loc[ind, "impact"] = obj.impact
        hgvs_df.loc[ind, "blosum62"] = obj.blosum62
        hgvs_df.loc[ind, "is_clinical"] = obj.clinical
        hgvs_df.loc[ind, "clinical_id"] = obj.clinical_id
        hgvs_df.loc[ind, "clinical_significance"] = obj.clinical_sig
        hgvs_df.loc[ind, "cosmic_id"] = obj.cosmic_id
        hgvs_df.loc[ind, "clinvar_id"] = obj.clinvar_id
        hgvs_df.loc[ind, "ancestral_populations"] = obj.ancestral_populations

    vep_df = hgvs_df.copy()
    vep_df = vep_df.drop_duplicates()

    return vep_df


def annotate_edits(
    ensembl_object: Ensembl, vep_df: pandas.DataFrame, uniprot_id: str | None
) -> pandas.DataFrame:
    """
    Annotate base edits with UniProt protein domain and PTM information.

    Enriches the VEP (Variant Effect Predictor) dataframe with additional protein
    annotations from UniProt including protein domains, post-translational
    modifications, and sequence mapping information.

    :param ensembl_object: Ensembl object containing gene and protein sequence data
    :param vep_df: DataFrame containing VEP annotation results
    :param uniprot_id: User-specified UniProt accession ID (if provided)

    :return: DataFrame enriched with UniProt protein domain and PTM annotations

    .. note::
        This function performs sequence alignment between Ensembl and UniProt
        protein sequences when direct mapping is not available. It adds columns
        for protein domains, curated domains, and post-translational modifications.
    """

    uniprot_df = vep_df.copy()
    uniprot_df["Domain"] = None
    uniprot_df["curated_Domain"] = None
    uniprot_df["PTM"] = None
    uniprot = None
    if uniprot_id is not None:
        uniprot = uniprot_id
    else:
        uniprot_list = [
            x for x in list(vep_df["swissprot_vep"].unique()) if not pandas.isna(x)
        ]
        if len(uniprot_list) == 1:
            uniprot = uniprot_list[0]

    ensembl_p = [
        x for x in list(vep_df["Protein_ID"].unique()) if pandas.isna(x) is not True
    ][0]
    seq_mapping = ensembl_object.extract_uniprot_info(
        ensembl_pid=ensembl_p, uniprot=uniprot
    )

    if seq_mapping:
        if uniprot:
            if len(uniprot.split(".")) > 1:
                uniprot = uniprot.split(".")[0].split("-")[0]
            else:
                print("Uniprot ID: %s" % uniprot)

            smap = seq_mapping[uniprot]
            obj = Uniprot(uniprotid=uniprot)
            obj.extract_uniprot()

            uniprot_df["Protein_Position"] = None
            for ind, row in uniprot_df.iterrows():
                ptm, domain, c_domain = None, None, None
                if (
                    row["Protein_Position_ensembl"] is not None
                    and pandas.isna(row["Protein_Position_ensembl"]) is False
                ):

                    # First check if ensembl and uniprot sequences have same indices
                    if len(row["Protein_Position_ensembl"].split(";")) == 1:
                        position = int(row["Protein_Position_ensembl"])
                        if position in smap.keys():
                            uniprot_df.loc[ind, "Protein_Position"] = str(
                                smap[position]
                            )

                    elif len(row["Protein_Position_ensembl"].split(";")) > 1:
                        pos_text_list = list()
                        for position in row["Protein_Position_ensembl"].split(";"):
                            position = int(position)
                            if position in smap.keys():
                                pos_text_list.append(str(smap[position]))

                        uniprot_df.loc[ind, "Protein_Position"] = ";".join(
                            pos_text_list
                        )

                    ptms, domains = list(), list()
                    if uniprot_df.loc[ind, "Protein_Position"] is not None:
                        for position in str(
                            uniprot_df.loc[ind, "Protein_Position"]
                        ).split(";"):
                            if (
                                position is not None
                                and pandas.isna(position) is False
                                and position != "None"
                                and position != ""
                                and type(position) != float
                            ):
                                dom = obj.find_domain(int(position), row["Edited_AA"])
                                phos = obj.find_ptm_site(
                                    "phosphorylation", int(position), row["Edited_AA"]
                                )
                                meth = obj.find_ptm_site(
                                    "methylation", int(position), row["Edited_AA"]
                                )
                                ubi = obj.find_ptm_site(
                                    "ubiquitination", int(position), row["Edited_AA"]
                                )
                                acet = obj.find_ptm_site(
                                    "acetylation", int(position), row["Edited_AA"]
                                )
                                if dom is not None:
                                    domains.append(dom + "-" + position)
                                if phos is not None:
                                    ptms.append(phos + "-" + position)
                                if meth is not None:
                                    ptms.append(meth + "-" + position)
                                if ubi is not None:
                                    ptms.append(ubi + "-" + position)
                                if acet is not None:
                                    ptms.append(acet + "-" + position)
                    if ptms:
                        ptm = ";".join([i for i in ptms])
                    if domains:
                        domain = ";".join([i for i in domains])
                        c_domain = ";".join(
                            ["-".join(i.split("-")[:-1]) for i in domains]
                        )
                uniprot_df.loc[ind, "Domain"] = domain
                uniprot_df.loc[ind, "curated_Domain"] = c_domain
                uniprot_df.loc[ind, "PTM"] = ptm
        else:
            print("No Uniprot ID can be found.")
    return uniprot_df


def extract_pis(pis: str) -> list[int] | None:
    """
    Extract protein interaction site positions from YULab data format.

    Parses protein interaction site strings from the YULab database format
    to extract individual amino acid positions involved in protein-protein
    interactions. Handles various formatting patterns including ranges.

    :param pis: Protein interaction site string from YULab data

    :return: List of amino acid positions involved in protein interactions

    .. note::
        This function handles multiple string formats including individual
        positions, ranges (e.g., "15-20"), and bracketed notations from
        the YULab protein interaction database.
    """
    sites: list[int] = list()
    if pis != "[]":
        for site in pis.split(","):
            if site[0] == "[" and site[-1] != "]":
                s_first = site[1:]
                if len(s_first.split("-")) > 1:
                    for s in list(
                        range(
                            int(s_first.split("-")[0]), int(s_first.split("-")[1]) + 1
                        )
                    ):
                        sites.append(int(s))
                else:
                    sites.append(int(s_first))
            elif site[-1] == "]" and site[0] != "[":
                s_last = site[:-1]
                if len(s_last.split("-")) > 1:
                    for s in list(
                        range(int(s_last.split("-")[0]), int(s_last.split("-")[1]) + 1)
                    ):
                        sites.append(int(s))
                else:
                    sites.append(int(s_last))
            elif site[-1] == "]" and site[0] == "[":
                s_only = site[1:-1]
                if len(s_only.split("-")) > 1:
                    for s in list(
                        range(int(s_only.split("-")[0]), int(s_only.split("-")[1]) + 1)
                    ):
                        sites.append(int(s))
                else:
                    sites.append(int(s_only))
            else:
                if len(site.split("-")) > 1:
                    for s in list(
                        range(int(site.split("-")[0]), int(site.split("-")[1]) + 1)
                    ):
                        sites.append(int(s))
                else:
                    sites.append(int(site))
        sites.sort()
        return sites
    else:
        return None


def collect_pis(
    uniprot: str, yulab_df: pandas.DataFrame
) -> dict[int, list[dict[str, str]]]:
    """
    Collecting protein interaction position for a given uniprot id

    :param uniprot: Uniprot ID
    :param yulab_df: A DataFrame corresponding to the YULab protein interaction data

    :return: positional dictionary specify the position and their source and partner (PDB/I3D/ECLAIR)
    """
    pis_dict = dict()
    df1 = yulab_df[yulab_df.P1 == uniprot]
    df2 = yulab_df[yulab_df.P2 == uniprot]

    if len(df1.index) != 0:
        for partner, partner_df in df1.groupby("P2"):
            for p_ind, p_row in partner_df.iterrows():
                source = p_row["Source"]
                interface_indices = extract_pis(p_row["P1_IRES"])
                if interface_indices is not None:
                    for ind in interface_indices:
                        if ind not in pis_dict.keys():
                            pis_dict[ind] = [{"partner": partner, "source": source}]
                        else:
                            t = pis_dict[ind]
                            if {"partner": partner, "source": source} not in t:
                                t.append({"partner": partner, "source": source})
                            pis_dict[ind] = t
    if len(df2.index) != 0:
        for partner, partner_df in df1.groupby("P1"):
            for p_ind, p_row in partner_df.iterrows():
                source = p_row["Source"]
                interface_indices = extract_pis(p_row["P2_IRES"])
                if interface_indices is not None:
                    for ind in interface_indices:
                        if ind not in pis_dict.keys():
                            pis_dict[ind] = [{"partner": partner, "source": source}]
                        else:
                            t = pis_dict[ind]
                            if {"partner": partner, "source": source} not in t:
                                t.append({"partner": partner, "source": source})
                            pis_dict[ind] = t
    return pis_dict


def disrupt_interface(
    uniprot: str, pos: int, yulab_df: pandas.DataFrame
) -> tuple[str | None, str | None, str | None]:
    """
    Checking if the given position disrupts the interfaces in the given uniprot

    :param uniprot: Uniprot ID
    :param pos: Uniprot index
    :param yulab_df: A DataFrame corresponding to the YULab protein interaction data

    :return: A 3-tuple of effected PDB partners, effected I3D partners, effected
        ECLAIR partners, if any
    """
    d = collect_pis(uniprot, yulab_df)
    if pos in d.keys():
        pdb_partner_list = list()
        i3d_partner_list = list()
        eclair_partner_list = list()
        for k in d[pos]:
            if k["source"] == "PDB":
                if k["partner"] not in pdb_partner_list:
                    pdb_partner_list.append(k["partner"])
            elif k["source"] == "I3D":
                if k["partner"] not in i3d_partner_list:
                    i3d_partner_list.append(k["partner"])
            elif k["source"] == "ECLAIR":
                if k["partner"] not in eclair_partner_list:
                    eclair_partner_list.append(k["partner"])
        if len(pdb_partner_list) == 0:
            pdb = None
        else:
            pdb = ",".join(pdb_partner_list)
        if len(i3d_partner_list) == 0:
            i3d = None
        else:
            i3d = ",".join(i3d_partner_list)
        if len(eclair_partner_list) == 0:
            eclair = None
        else:
            eclair = ",".join(eclair_partner_list)
        return pdb, i3d, eclair
    else:
        return None, None, None


def annotate_interface(
    annotated_edit_df: pandas.DataFrame,
    uniprot_id: t.Optional[str],
    yulab_df: pandas.DataFrame,
) -> pandas.DataFrame:
    """
    Add Interactome Insider protein interface information for edgetic perturbation.

    :param annotated_edit_df: DataFrame created with annotate_edits
    :param uniprot_id: UniProt Accession ID, if specified
    :param yulab_df: A DataFrame corresponding to the YULab protein interaction data

    :return: Added interface annotation on edit table
    """
    server_url = "https://www.ebi.ac.uk/proteins/api/proteins?"
    df = annotated_edit_df.copy()
    df["is_disruptive_interface_EXP"] = None
    df["is_disruptive_interface_MOD"] = None
    df["is_disruptive_interface_PRED"] = None
    df["disrupted_PDB_int_partners"] = None
    df["disrupted_I3D_int_partners"] = None
    df["disrupted_Eclair_int_partners"] = None
    df["disrupted_PDB_int_genes"] = None
    df["disrupted_I3D_int_genes"] = None
    df["disrupted_Eclair_int_genes"] = None

    if uniprot_id is not None:
        group_cols = ["uniprot_provided", "Protein_Position"]
    else:
        group_cols = ["swissprot_vep", "Protein_Position"]
    for group, group_df in df.groupby(group_cols):
        if (
            group[1] is not None
            and group[1] != "None"
            and pandas.isna(group[1]) == False
            and group[1] != ""
        ):
            if group[0] in list(yulab_df.P1) or group[0] in list(yulab_df.P2):
                all_pdb_partners, all_i3d_partners, all_eclair_partners = (
                    list(),
                    list(),
                    list(),
                )
                if len(group[1].split(";")) == 1:
                    pdb_partners, i3d_partners, eclair_partners = disrupt_interface(
                        uniprot=group[0], pos=int(group[1]), yulab_df=yulab_df
                    )
                    if pdb_partners is not None:
                        all_pdb_partners.append(pdb_partners)
                    if i3d_partners is not None:
                        all_i3d_partners.append(i3d_partners)
                    if eclair_partners is not None:
                        all_eclair_partners.append(eclair_partners)
                else:
                    for pos in group[1].split(";"):
                        pdb_partners, i3d_partners, eclair_partners = disrupt_interface(
                            uniprot=group[0], pos=int(pos), yulab_df=yulab_df
                        )
                        if pdb_partners is not None:
                            all_pdb_partners.append(pdb_partners)
                        if i3d_partners is not None:
                            all_i3d_partners.append(i3d_partners)
                        if eclair_partners is not None:
                            all_eclair_partners.append(eclair_partners)

                if all_pdb_partners:
                    df.loc[list(group_df.index), "is_disruptive_interface_EXP"] = True
                    df.loc[list(group_df.index), "disrupted_PDB_int_partners"] = (
                        ";".join(all_pdb_partners)
                    )
                    for uniprot in all_pdb_partners:
                        api_url = (
                            "offset=0&size=-1&accession=%s&reviewed=true&isoform=0"
                            % uniprot
                        )

                        api_request = requests.get(
                            server_url + api_url, headers={"Accept": "application/json"}
                        )

                        # Check the response of the server for the request
                        genes = list()
                        if api_request.status_code == 200:
                            for i in api_request.json():
                                for k in i["gene"]:
                                    gene_name = k["name"]["value"]
                                if gene_name not in genes:
                                    genes.append(gene_name)
                            genes = ";".join(genes)
                        else:
                            genes = None
                    df.loc[list(group_df.index), "disrupted_PDB_int_genes"] = genes
                else:
                    df.loc[list(group_df.index), "is_disruptive_interface_EXP"] = False
                    df.loc[list(group_df.index), "disrupted_PDB_int_partners"] = None
                    df.loc[list(group_df.index), "disrupted_PDB_int_genes"] = None
                if all_i3d_partners:
                    df.loc[list(group_df.index), "is_disruptive_interface_MOD"] = True
                    df.loc[list(group_df.index), "disrupted_I3D_int_partners"] = (
                        ";".join(all_i3d_partners)
                    )
                    for uniprot in all_i3d_partners:
                        api_url = (
                            "offset=0&size=-1&accession=%s&reviewed=true&isoform=0"
                            % uniprot
                        )

                        api_request = requests.get(
                            server_url + api_url, headers={"Accept": "application/json"}
                        )

                        # Check the response of the server for the request
                        genes = list()
                        if api_request.status_code == 200:
                            for i in api_request.json():
                                for k in i["gene"]:
                                    gene_name = k["name"]["value"]
                                if gene_name not in genes:
                                    genes.append(gene_name)
                            genes = ";".join(genes)
                        else:
                            genes = None
                    df.loc[list(group_df.index), "disrupted_I3D_int_genes"] = genes
                else:
                    df.loc[list(group_df.index), "is_disruptive_interface_MOD"] = False
                    df.loc[list(group_df.index), "disrupted_I3D_int_partners"] = None
                    df.loc[list(group_df.index), "disrupted_I3D_int_genes"] = None
                if all_eclair_partners:
                    df.loc[list(group_df.index), "is_disruptive_interface_PRED"] = True
                    df.loc[list(group_df.index), "disrupted_Eclair_int_partners"] = (
                        ";".join(all_eclair_partners)
                    )
                    for uniprot in all_eclair_partners:
                        api_url = (
                            "offset=0&size=-1&accession=%s&reviewed=true&isoform=0"
                            % uniprot
                        )

                        api_request = requests.get(
                            server_url + api_url, headers={"Accept": "application/json"}
                        )

                        # Check the response of the server for the request
                        genes = list()
                        if api_request.status_code == 200:
                            for i in api_request.json():
                                for k in i["gene"]:
                                    gene_name = k["name"]["value"]
                                if gene_name not in genes:
                                    genes.append(gene_name)
                            genes = ";".join(genes)
                        else:
                            genes = None
                    df.loc[list(group_df.index), "disrupted_Eclair_int_genes"] = genes
                else:
                    df.loc[list(group_df.index), "is_disruptive_interface_PRED"] = False
                    df.loc[list(group_df.index), "disrupted_Eclair_int_partners"] = None
                    df.loc[list(group_df.index), "disrupted_Eclair_int_genes"] = None
    return df


def rename_mutational_consequences(mutation_consequence):
    """
    TODO documentation
    """
    consequence = list()
    for consq in mutation_consequence.split(";"):
        if consq == "missense_variant":
            consequence.append("missense")
        if consq == "missense_mutation":
            consequence.append("missense")
        elif consq == "missense_variant_splice_region_variant":
            consequence.append("splice variant")
        elif consq == "splice_region_variant":
            consequence.append("splice variant")
        elif consq == "stop_retained_variant":
            consequence.append("synonymous")
        elif consq == "synonymous_variant":
            consequence.append("synonymous")
        elif consq == "splice_region_variant_synonymous_variant":
            consequence.append("splice variant")
        elif consq == "splice_acceptor_variant":
            consequence.append("splice variant")
        elif consq == "splice_donor_variant":
            consequence.append("splice variant")
        elif consq == "splice_region_variant_intron_variant":
            consequence.append("splice variant")
        elif consq == "splice_region_variant,intron_variant":
            consequence.append("splice variant")
        elif consq == "splice_donor_region_variant_intron_variant":
            consequence.append("splice variant")
        elif consq == "splice_polypyrimidine_tract_variant_intron_variant":
            consequence.append("splice variant")
        elif (
            consq
            == "splice_polypyrimidine_tract_variant_splice_region_variant_intron_variant"
        ):
            consequence.append("splice variant")
        elif consq == "splice_donor_5th_base_variant_intron_variant":
            consequence.append("splice variant")
        elif consq == "downstream_gene_variant":
            consequence.append("UTR")
        elif consq == "stop_gained_splice_region_variant":
            consequence.append("stop codon")
        elif consq == "stop_gained,splice_region_variant":
            consequence.append("stop codon")
        elif consq == "start_lost":
            consequence.append("start lost")
        elif consq == "stop_gained_start_lost":
            consequence.append("stop codon")
        elif consq == "upstream_gene_variant":
            consequence.append("promoter")
        elif consq == "intron_variant":
            consequence.append("intron")
        elif consq == "5_prime_UTR_variant":
            consequence.append("5'UTR")
        elif consq == "stop_gained":
            consequence.append("stop codon")
    consequence = ";".join(consequence)
    return consequence


def select_severe_effects(mutation_consequence):
    """
    TODO documentation
    """
    if (
        mutation_consequence is None
        or pandas.isna(mutation_consequence)
        or mutation_consequence == ""
    ):
        return ""
    elif "stop codon" in mutation_consequence.split(";"):
        return "stop codon"
    elif "start lost" in mutation_consequence.split(";"):
        return "start lost"
    elif "splice variant" in mutation_consequence.split(";"):
        return "splice variant"
    elif "missense" in mutation_consequence.split(";"):
        return "missense"
    elif "UTR" in mutation_consequence.split(";"):
        return "UTR"
    elif "intron" in mutation_consequence.split(";"):
        return "intron"
    elif "synonymous" in mutation_consequence.split(";"):
        return "synonymous"


def summarise_3di(list_of_partners):
    """
    TODO documentation
    """
    all_partners = list()
    if list_of_partners is not None:
        for partner_list in list_of_partners:
            if (
                partner_list is not None
                and pandas.isna(partner_list) == False
                and partner_list != []
            ):
                for partner in partner_list.split(";"):
                    if partner not in all_partners:
                        all_partners.append(partner)
    if len(all_partners) > 0:
        return ";".join(all_partners)
    else:
        return None


def summarise_guides(last_df: pandas.DataFrame) -> pandas.DataFrame:
    """
    Create summary report of gRNA guides with aggregated annotations.

    Consolidates all base editing analysis results by gRNA guide, aggregating
    multiple edit sites per guide and their associated annotations including
    protein effects, clinical significance, and functional predictions.

    :param last_df: Complete DataFrame with all edit annotations

    :return: Summary DataFrame with one row per unique gRNA guide

    .. note::
        This function groups edits by CRISPR_PAM_Sequence and aggregates
        all associated annotations, consequence predictions, and functional
        effects into a comprehensive summary for each guide.
    """
    summary_df = pandas.DataFrame(
        index=list(range(0, len(last_df.groupby(["CRISPR_PAM_Sequence"])))),
        columns=[
            "Hugo_Symbol",
            "CRISPR_PAM_Sequence",
            "CRISPR_PAM_Location",
            "gRNA_Target_Sequence",
            "gRNA_Target_Location",
            "gRNA_flanking_sequences",
            "Edit_Location",
            "Direction",
            "Transcript_ID",
            "Exon_ID",
            "Protein_ID",
            "guide_in_CDS",
            "Edit_in_Exon",
            "Edit_in_CDS",
            "mutation_on_guide",
            "guide_change_mutation",
            "mutation_on_window",
            "mutation_on_PAM",
            "# Edits/guide",
            "Poly_T",
            "GC%",
            "allele",
            "cDNA_Change",
            "CDS_Position",
            "Protein_Position_ensembl",
            "Protein_Position",
            "Protein_Change",
            "Edited_AA",
            "Edited_AA_Prop",
            "New_AA",
            "New_AA_Prop",
            "is_stop",
            "is_synonymous",
            "proline_addition",
            "variant_classification",
            "consequence_terms",
            "most_severe_consequence",
            "variant_biotype",
            "Regulatory_ID",
            "Motif_ID",
            "TFs_on_motif",
            "polyphen_prediction",
            "sift_prediction",
            "impact",
            "is_clinical",
            "clinical_id",
            "clinical_significance",
            "cosmic_id",
            "clinvar_id",
            "ancestral_populations",
            "swissprot_vep",
            "uniprot_provided",
            "Domain",
            "curated_Domain",
            "PTM",
            "is_disruptive_interface_EXP",
            "disrupted_PDB_int_partners",
            "disrupted_PDB_int_genes",
            "is_disruptive_interface_MOD",
            "disrupted_I3D_int_partners",
            "disrupted_I3D_int_genes",
            "is_disruptive_interface_PRED",
            "disrupted_Eclair_int_partners",
            "disrupted_Eclair_int_genes",
        ],
    )
    # cosmic_freq

    i = 0
    for guide, guide_df in last_df.groupby("CRISPR_PAM_Sequence"):

        summary_df.loc[i, "Hugo_Symbol"] = ";".join(
            [
                str(x)
                for x in list(guide_df.Hugo_Symbol.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )

        summary_df.loc[i, "CRISPR_PAM_Sequence"] = ";".join(
            [
                str(x)
                for x in list(guide_df.CRISPR_PAM_Sequence.unique())
                if x is not None and pandas.isna(x) is False
            ]
            if guide_df.CRISPR_PAM_Sequence.unique() is not None
            else ""
        )

        summary_df.loc[i, "CRISPR_PAM_Location"] = ";".join(
            [
                str(x)
                for x in list(guide_df.CRISPR_PAM_Location.unique())
                if x is not None and pandas.isna(x) is False
            ]
            if guide_df.CRISPR_PAM_Location.unique() is not None
            else ""
        )

        summary_df.loc[i, "gRNA_Target_Sequence"] = ";".join(
            [
                str(x)
                for x in list(guide_df.gRNA_Target_Sequence.unique())
                if x is not None and pandas.isna(x) is False
            ]
            if guide_df.gRNA_Target_Sequence.unique() is not None
            else ""
        )

        summary_df.loc[i, "gRNA_flanking_sequences"] = ";".join(
            [
                str(x)
                for x in list(guide_df.gRNA_flanking_sequences.unique())
                if x is not None and pandas.isna(x) is False
            ]
            if pandas.isna(guide_df.gRNA_flanking_sequences.unique()) is False
            else ""
        )

        summary_df.loc[i, "gRNA_Target_Location"] = ";".join(
            [
                str(x)
                for x in list(guide_df.gRNA_Target_Location.unique())
                if x is not None and pandas.isna(x) is False
            ]
            if guide_df.gRNA_Target_Location.unique() is not None
            else ""
        )

        summary_df.loc[i, "Edit_Location"] = ";".join(
            [
                str(x)
                for x in list(guide_df.Edit_Location.unique())
                if x is not None and pandas.isna(x) is False
            ]
            if guide_df.Edit_Location.unique() is not None
            else ""
        )

        summary_df.loc[i, "Direction"] = ";".join(
            [
                x
                for x in list(guide_df.Direction.unique())
                if x is not None and pandas.isna(x) is False
            ]
            if guide_df.Direction.unique() is not None
            else ""
        )

        summary_df.loc[i, "Transcript_ID"] = ";".join(
            [
                x
                for x in list(guide_df.Transcript_ID.unique())
                if x is not None and pandas.isna(x) is False
            ]
            if guide_df.Transcript_ID.unique() is not None
            else ""
        )

        summary_df.loc[i, "Exon_ID"] = ";".join(
            [
                x
                for x in list(guide_df.Exon_ID.unique())
                if x is not None and pandas.isna(x) is False
            ]
            if guide_df.Exon_ID.unique() is not None
            else ""
        )

        summary_df.loc[i, "Protein_ID"] = ";".join(
            [
                x
                for x in list(guide_df.Protein_ID.unique())
                if x is not None and pandas.isna(x) is False
            ]
            if guide_df.Protein_ID.unique() is not None
            else ""
        )

        if (
            guide_df[~pandas.isna(guide_df.Regulatory_ID)].Regulatory_ID.unique()
            is not None
            and type(guide_df.Regulatory_ID) != float
            and list(
                guide_df[~pandas.isna(guide_df.Regulatory_ID)].Regulatory_ID.unique()
            )
        ):
            summary_df.loc[i, "Regulatory_ID"] = ";".join(
                [
                    x
                    for x in list(guide_df.Regulatory_ID.unique())
                    if x is not None and pandas.isna(x) is False
                ]
            )
        else:
            summary_df.loc[i, "Regulatory_ID"] = None

        if (
            guide_df[~pandas.isna(guide_df.Motif_ID)].Motif_ID.unique() is not None
            and type(guide_df.Motif_ID) != float
            and list(guide_df[~pandas.isna(guide_df.Motif_ID)].Motif_ID.unique())
        ):
            summary_df.loc[i, "Motif_ID"] = ";".join(
                [
                    x
                    for x in list(guide_df.Motif_ID.unique())
                    if x is not None and pandas.isna(x) is False
                ]
            )
        else:
            summary_df.loc[i, "Motif_ID"] = None

        if (
            guide_df[~pandas.isna(guide_df.TFs_on_motif)].TFs_on_motif.unique()
            is not None
            and type(guide_df.TFs_on_motif) != float
            and list(
                guide_df[~pandas.isna(guide_df.TFs_on_motif)].TFs_on_motif.unique()
            )
        ):
            summary_df.loc[i, "TFs_on_motif"] = ";".join(
                [
                    x
                    for x in list(guide_df.TFs_on_motif.unique())
                    if x is not None and pandas.isna(x) is False
                ]
            )
        else:
            summary_df.loc[i, "TFs_on_motif"] = None

        summary_df.loc[i, "guide_in_CDS"] = (
            True if True in guide_df.guide_in_CDS.unique() else False
        )
        summary_df.loc[i, "Edit_in_Exon"] = (
            True if True in guide_df.Edit_in_Exon.unique() else False
        )
        summary_df.loc[i, "Edit_in_CDS"] = (
            True if True in guide_df.Edit_in_CDS.unique() else False
        )
        summary_df.loc[i, "mutation_on_guide"] = (
            True if True in guide_df.mutation_on_guide.unique() else False
        )
        summary_df.loc[i, "guide_change_mutation"] = (
            True if True in guide_df.guide_change_mutation.unique() else False
        )
        summary_df.loc[i, "mutation_on_window"] = (
            True if True in guide_df.mutation_on_window.unique() else False
        )
        summary_df.loc[i, "mutation_on_PAM"] = (
            True if True in guide_df.mutation_on_PAM.unique() else False
        )

        summary_df.loc[i, "# Edits/guide"] = guide_df["# Edits/guide"].unique()[0]

        summary_df.loc[i, "Poly_T"] = (
            True if True in guide_df.Poly_T.unique() else False
        )
        summary_df.loc[i, "GC%"] = "".join([str(x) for x in guide_df["GC%"].unique()])

        summary_df.loc[i, "allele"] = ";".join(
            [
                x
                for x in list(guide_df.allele.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )

        summary_df.loc[i, "cDNA_Change"] = ";".join(
            [
                x
                for x in list(guide_df.cDNA_Change.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )

        summary_df.loc[i, "CDS_Position"] = ";".join(
            [
                x
                for x in list(guide_df.CDS_Position.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )

        summary_df.loc[i, "Protein_Position"] = ";".join(
            [
                str(x)
                for x in list(guide_df.Protein_Position.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )

        summary_df.loc[i, "Protein_Position_ensembl"] = ";".join(
            [
                str(x)
                for x in list(guide_df.Protein_Position_ensembl.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )
        summary_df.loc[i, "Protein_Change"] = ";".join(
            [
                x
                for x in list(guide_df.Protein_Change.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )

        summary_df.loc[i, "Edited_AA"] = ";".join(
            [
                x
                for x in list(guide_df.Edited_AA.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )

        summary_df.loc[i, "New_AA"] = ";".join(
            [
                x
                for x in list(guide_df.New_AA.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )

        summary_df.loc[i, "Edited_AA_Prop"] = ";".join(
            [
                x
                for x in list(guide_df.Edited_AA_Prop.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )

        summary_df.loc[i, "New_AA_Prop"] = ";".join(
            [
                x
                for x in list(guide_df.New_AA_Prop.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )

        summary_df.loc[i, "swissprot_vep"] = ";".join(
            [
                x
                for x in list(guide_df.swissprot_vep.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )

        if (
            guide_df.uniprot_provided.unique() is not None
            and pandas.isna(guide_df.uniprot_provided) is False
            and list(guide_df.uniprot_provided.unique())
        ):
            summary_df.loc[i, "uniprot_provided"] = ";".join(
                [
                    x
                    for x in list(guide_df.uniprot_provided.unique())
                    if x is not None and pandas.isna(x) is False
                ]
            )
        else:
            summary_df.loc[i, "uniprot_provided"] = None

        summary_df.loc[i, "variant_classification"] = ";".join(
            [
                x
                for x in list(guide_df.variant_classification.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )

        summary_df.loc[i, "variant_biotype"] = ";".join(
            [
                x
                for x in list(guide_df.variant_biotype.unique())
                if x is not None and pandas.isna(x) is False
            ]
        )

        summary_df.loc[i, "consequence_terms"] = ";".join(
            [
                select_severe_effects(x)
                for x in [
                    select_severe_effects(i)
                    for i in [
                        rename_mutational_consequences(c)
                        for c in [
                            x
                            for x in list(guide_df.consequence_terms.unique())
                            if x is not None and pandas.isna(x) == False
                        ]
                    ]
                ]
            ]
        )

        summary_df.loc[i, "most_severe_consequence"] = select_severe_effects(
            ";".join(
                [
                    rename_mutational_consequences(x)
                    for x in guide_df.most_severe_consequence.unique()
                    if x is not None and pandas.isna(x) == False
                ]
            )
        )

        summary_df.loc[i, "is_stop"] = (
            True if "stop codon" in guide_df.most_severe_consequence.unique() else False
        )
        summary_df.loc[i, "is_synonymous"] = (
            True if "synonymous" in guide_df.most_severe_consequence.unique() else False
        )
        summary_df.loc[i, "proline_addition"] = (
            True if True in guide_df.proline_addition.unique() else False
        )

        if (
            guide_df[
                ~pandas.isna(guide_df.polyphen_prediction)
            ].polyphen_prediction.unique()
            is not None
            and type(guide_df.polyphen_prediction) != float
            and list(
                guide_df[
                    ~pandas.isna(guide_df.polyphen_prediction)
                ].polyphen_prediction.unique()
            )
        ):
            summary_df.loc[i, "polyphen_prediction"] = ";".join(
                [
                    str(x)
                    for x in list(guide_df.polyphen_prediction.unique())
                    if x is not None and pandas.isna(x) is False
                ]
            )
        else:
            summary_df.loc[i, "polyphen_prediction"] = None

        summary_df.loc[i, "sift_prediction"] = ";".join(
            [
                x
                for x in list(guide_df.sift_prediction.unique())
                if x is not None and pandas.isna(x) == False and type(x) != float
            ]
        )

        summary_df.loc[i, "impact"] = ";".join(
            [
                x
                for x in list(guide_df.impact.unique())
                if x is not None and pandas.isna(x) == False and type(x) != float
            ]
        )

        summary_df.loc[i, "is_clinical"] = (
            True if True in guide_df.is_clinical.unique() else False
        )

        if (
            guide_df[~pandas.isna(guide_df.clinical_id)].clinical_id.unique()
            is not None
            and type(guide_df.clinical_id) != float
            and list(guide_df[~pandas.isna(guide_df.clinical_id)].clinical_id.unique())
        ):
            summary_df.loc[i, "clinical_id"] = ";".join(
                [
                    x
                    for x in list(guide_df.clinical_id.unique())
                    if x is not None and pandas.isna(x) == False and type(x) != float
                ]
            )
        else:
            summary_df.loc[i, "clinical_id"] = None

        if (
            guide_df[
                ~pandas.isna(guide_df.clinical_significance)
            ].clinical_significance.unique()
            is not None
            and type(guide_df.clinical_significance) != float
            and list(
                guide_df[
                    ~pandas.isna(guide_df.clinical_significance)
                ].clinical_significance.unique()
            )
        ):
            summary_df.loc[i, "clinical_significance"] = ";".join(
                [
                    x
                    for x in list(guide_df.clinical_significance.unique())
                    if x is not None and type(x) != float
                ]
            )
        else:
            summary_df.loc[i, "clinical_significance"] = None

        if (
            guide_df[~pandas.isna(guide_df.cosmic_id)].cosmic_id.unique() is not None
            and type(guide_df.cosmic_id) != float
            and list(guide_df[~pandas.isna(guide_df.cosmic_id)].cosmic_id.unique())
        ):
            summary_df.loc[i, "cosmic_id"] = ";".join(
                [
                    str(x)
                    for x in list(guide_df.cosmic_id.unique())
                    if x is not None and type(x) != float
                ]
            )
        else:
            summary_df.loc[i, "cosmic_id"] = None

        if (
            guide_df[~pandas.isna(guide_df.clinvar_id)].clinvar_id.unique() is not None
            and type(guide_df.clinvar_id) != float
            and list(guide_df[~pandas.isna(guide_df.clinvar_id)].clinvar_id.unique())
        ):
            summary_df.loc[i, "clinvar_id"] = ";".join(
                [
                    str(x)
                    for x in list(guide_df.clinvar_id.unique())
                    if x is not None and type(x) != float
                ]
            )
        else:
            summary_df.loc[i, "clinvar_id"] = None
        if (
            guide_df[
                ~pandas.isna(guide_df.ancestral_populations)
            ].ancestral_populations.unique()
            is not None
            and type(guide_df.ancestral_populations) != float
            and list(
                guide_df[
                    ~pandas.isna(guide_df.ancestral_populations)
                ].ancestral_populations.unique()
            )
        ):
            summary_df.loc[i, "ancestral_populations"] = ";".join(
                [
                    x
                    for x in list(guide_df.ancestral_populations.unique())
                    if x is not None and type(x) != float
                ]
            )
        else:
            summary_df.loc[i, "ancestral_populations"] = None
        if (
            guide_df[~pandas.isna(guide_df.Domain)].Domain.unique() is not None
            and type(guide_df.Domain) != float
            and list(guide_df[~pandas.isna(guide_df.Domain)].Domain.unique())
        ):
            summary_df.loc[i, "Domain"] = ";".join(
                [
                    x
                    for x in list(guide_df.Domain.unique())
                    if x is not None and type(x) != float
                ]
            )
        else:
            summary_df.loc[i, "Domain"] = None
        if (
            guide_df[~pandas.isna(guide_df.curated_Domain)].curated_Domain.unique()
            is not None
            and type(guide_df.curated_Domain) != float
            and list(
                guide_df[~pandas.isna(guide_df.curated_Domain)].curated_Domain.unique()
            )
        ):
            summary_df.loc[i, "curated_Domain"] = ";".join(
                [
                    x
                    for x in list(guide_df.curated_Domain.unique())
                    if x is not None and type(x) != float
                ]
            )
        else:
            summary_df.loc[i, "curated_Domain"] = None
        if (
            guide_df[~pandas.isna(guide_df.PTM)].PTM.unique() is not None
            and type(guide_df.PTM) != float
            and list(guide_df[~pandas.isna(guide_df.PTM)].PTM.unique())
        ):
            summary_df.loc[i, "PTM"] = ";".join(
                [
                    x
                    for x in list(guide_df.PTM.unique())
                    if x is not None and type(x) != float
                ]
            )
        else:
            summary_df.loc[i, "PTM"] = None

        summary_df.loc[i, "is_disruptive_interface_EXP"] = (
            True if True in guide_df.is_disruptive_interface_EXP.unique() else False
        )
        summary_df.loc[i, "disrupted_PDB_int_partners"] = summarise_3di(
            list(guide_df.disrupted_PDB_int_partners.unique())
        )
        summary_df.loc[i, "disrupted_PDB_int_genes"] = summarise_3di(
            list(guide_df.disrupted_PDB_int_genes.unique())
        )
        summary_df.loc[i, "is_disruptive_interface_MOD"] = (
            True if True in guide_df.is_disruptive_interface_MOD.unique() else False
        )
        summary_df.loc[i, "disrupted_I3D_int_partners"] = summarise_3di(
            list(guide_df.disrupted_I3D_int_partners.unique())
        )
        summary_df.loc[i, "disrupted_I3D_int_genes"] = summarise_3di(
            list(guide_df.disrupted_I3D_int_genes.unique())
        )
        summary_df.loc[i, "is_disruptive_interface_PRED"] = (
            True if True in guide_df.is_disruptive_interface_PRED.unique() else False
        )
        summary_df.loc[i, "disrupted_Eclair_int_partners"] = summarise_3di(
            list(guide_df.disrupted_Eclair_int_partners.unique())
        )
        summary_df.loc[i, "disrupted_Eclair_int_genes"] = summarise_3di(
            list(guide_df.disrupted_Eclair_int_genes.unique())
        )
        i += 1

    return summary_df
