# Copyright (C) 2025 Genome Research Ltd.

import pandas, re, requests
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

from BEstimate import constants


class Ensembl:
    """
    A class for interacting with the Ensembl REST API to retrieve genomic information.

    This class handles gene and transcript data retrieval from Ensembl including
    genomic sequences, transcript information, exon boundaries, and coding sequences.
    It provides methods to extract and analyze genomic features for base editor analysis.

    :param hugo_symbol: HGNC gene symbol (e.g., 'TP53')
    :type hugo_symbol: str
    :param assembly: Genome assembly version ('GRCh38' or 'hg19')
    :type assembly: str

    :ivar hugo_symbol: HGNC gene symbol
    :ivar assembly: Genome assembly version
    :ivar server: Base URL for Ensembl REST API
    :ivar gene_id: Ensembl gene identifier
    :ivar info_dict: Dictionary containing transcript and exon information
    :ivar sequence: Gene genomic sequence
    :ivar flan_sequence: Gene sequence with flanking regions
    :ivar chromosome: Chromosome where the gene is located
    :ivar strand: Gene strand orientation (1 or -1)
    :ivar gene_range: List containing gene start and end positions
    :ivar flan_gene_range: List containing gene positions with flanking regions
    """

    def __init__(self, hugo_symbol: str, assembly: str) -> None:
        """
        Initialize Ensembl object with gene symbol and assembly.

        :param hugo_symbol: HGNC gene symbol
        :param assembly: Genome assembly version
        """
        self.hugo_symbol = hugo_symbol
        self.assembly = assembly
        self.server = (
            "http://grch37.rest.ensembl.org"
            if self.assembly == "hg19"
            else "https://rest.ensembl.org"
        )
        self.gene_id = ""
        self.info_dict = dict()
        self.sequence, self.flan_sequence = None, None
        self.right_sequence_analysis, self.flan_right_sequence_analysis = None, None
        self.left_sequence_analysis, self.flan_left_sequence_analysis = None, None
        self.chromosome, self.strand = None, None
        self.gene_range, self.flan_gene_range = list(), list()
        self.p_sequence = None

    def extract_gene_id(self) -> int:
        """
        Extract Ensembl gene identifier from gene symbol.

        Retrieves the Ensembl gene ID corresponding to the HGNC gene symbol
        using the Ensembl REST API. Validates that the gene is on a standard
        chromosome (1-22, X, Y).

        :return: 1 if gene ID found successfully, 0 if not found

        .. note::
            This method prints status messages to stdout and sets the gene_id attribute.
        """

        hugo_ensembl = "/xrefs/symbol/homo_sapiens/%s?" % self.hugo_symbol

        print("Request to Ensembl REST API for Ensembl Gene ID:")
        gene_request = requests.get(
            self.server + hugo_ensembl, headers={"Content-Type": "application/json"}
        )

        if gene_request.status_code != 200:
            print("No response from ensembl!\n")

        for x in gene_request.json():
            if x["id"][:4] == "ENSG":

                info_ensembl = self.server + "/lookup/id/%s?expand=1" % x["id"]
                info_request = requests.get(
                    info_ensembl, headers={"Content-Type": "application/json"}
                )

                if info_request.json()["display_name"] == self.hugo_symbol:
                    seq_ensembl = self.server + "/sequence/id/%s?" % x["id"]
                    seq_request = requests.get(
                        seq_ensembl, headers={"Content-Type": "text/x-fasta"}
                    )
                    print(seq_request.text.split("\n")[0])
                    chr = seq_request.text.split("\n")[0].split(":")[2].strip()
                    try:
                        if chr != "X" and chr != "Y":
                            int(chr)
                            self.gene_id = x["id"]
                        elif chr == "X" or chr == "Y":
                            self.gene_id = x["id"]
                    except ValueError:
                        print(" ")

        if self.gene_id != "":
            print("Ensembl Gene ID: %s\n" % self.gene_id)
            return 1
        else:
            return 0

    def extract_sequence(self, gene_id: str, mutations: list | None) -> None:
        """
        Extract genomic sequence for the gene from Ensembl.

        Retrieves the gene sequence and processes it for analysis, including
        handling strand orientation and incorporating user-provided mutations.

        :param gene_id: Ensembl gene identifier
        :param mutations: List of genomic mutations to incorporate into the sequence

        .. note::
            This method sets multiple sequence-related attributes and handles
            reverse complement for negative strand genes. Mutations are applied
            if provided and validated against the reference sequence.
        """

        seq_ensembl = self.server + "/sequence/id/%s?" % gene_id
        seq_flan_ensembl = (
            self.server + "/sequence/id/%s?expand_3prime=23;expand_5prime=23" % gene_id
        )

        print("Request to Ensembl REST API for sequence information")
        seq_request = requests.get(
            seq_ensembl, headers={"Content-Type": "text/x-fasta"}
        )
        seq_flan_request = requests.get(
            seq_flan_ensembl, headers={"Content-Type": "text/x-fasta"}
        )

        if seq_request.status_code != 200 and seq_flan_request.status_code != 200:
            print("No response from ensembl sequence!\n")

        # Sequence
        label_line = seq_request.text.split("\n")[0]
        if label_line[0] != "{":
            print(
                "The location of the interested gene: %s\n" % label_line.split(" ")[1]
            )
            flan_label_line = seq_flan_request.text.split("\n")[0]
            self.sequence = "".join(seq_request.text.split("\n")[1:])
            self.flan_sequence = "".join(seq_flan_request.text.split("\n")[1:])
            self.gene_range = [
                int(label_line.split(":")[-3]),
                int(label_line.split(":")[-2]),
            ]
            self.flan_gene_range = [
                int(flan_label_line.split(":")[-3]),
                int(flan_label_line.split(":")[-2]),
            ]
            self.strand = int(label_line.split(":")[-1].strip())
            self.chromosome = label_line.split(":")[2].strip()

            # If strand is -1, the sequence has been reversed to be in 5'->3' direction
            # The genomic location should be reverse to match with the sequence too.
            if self.strand == -1:
                self.gene_range = [self.gene_range[1], self.gene_range[0]]
            if self.strand == -1:
                self.flan_gene_range = [
                    self.flan_gene_range[1],
                    self.flan_gene_range[0],
                ]

            # Preparation of the Ensembl sequence for analysis

            nucleotide_dict = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

            if mutations is None:

                if self.strand == 1:
                    self.right_sequence_analysis = self.sequence
                    self.flan_right_sequence_analysis = self.flan_sequence
                    self.left_sequence_analysis = "".join(
                        [nucleotide_dict[n] for n in self.sequence[::-1]]
                    )
                    self.flan_left_sequence_analysis = "".join(
                        [nucleotide_dict[n] for n in self.flan_sequence[::-1]]
                    )

                elif self.strand == -1:
                    self.left_sequence_analysis = self.sequence
                    self.flan_left_sequence_analysis = self.flan_sequence
                    self.right_sequence_analysis = "".join(
                        [nucleotide_dict[n] for n in self.left_sequence_analysis[::-1]]
                    )
                    self.flan_right_sequence_analysis = "".join(
                        [
                            nucleotide_dict[n]
                            for n in self.flan_left_sequence_analysis[::-1]
                        ]
                    )

            else:
                for mutation in mutations:
                    print(mutation)
                    if int(mutation.split(":")[0]) == int(self.chromosome):
                        if mutation.split(":")[1].split(".")[0] == "g":
                            alteration = mutation.split(":")[1].split(".")[1]
                            mutation_location = int(
                                re.match(
                                    "([0-9]+)([a-z]+)", alteration.split(">")[0], re.I
                                ).groups()[0]
                            )
                            altered_nuc = re.match(
                                "([0-9]+)([a-z]+)", alteration.split(">")[0], re.I
                            ).groups()[1]
                            new_nuc = alteration.split(">")[1]

                            # Check if altered nucleotide is in the given location
                            if self.strand == 1:
                                genomic_start = int(self.gene_range[0])
                                genomic_flan_start = int(self.flan_gene_range[0])

                                if (
                                    self.sequence[mutation_location - genomic_start]
                                    == altered_nuc
                                ):
                                    # Altered nucleotide in the given mutation fits with the sequence
                                    s = list(self.sequence)
                                    s[mutation_location - genomic_start] = new_nuc
                                    self.sequence = "".join(s)
                                else:
                                    print(
                                        "\nGiven mutation location does not fit with the sequence."
                                        "Nucleotides are different.\n"
                                    )

                                if (
                                    self.flan_sequence[
                                        mutation_location - genomic_flan_start
                                    ]
                                    == altered_nuc
                                ):
                                    # Altered nucleotide in the given mutation fits with the sequence
                                    s = list(self.flan_sequence)
                                    s[mutation_location - genomic_flan_start] = new_nuc
                                    self.flan_sequence = "".join(s)
                                else:
                                    print(
                                        "\nGiven mutation location does not fit with the sequence."
                                        "Nucleotides are different.\n"
                                    )

                                self.right_sequence_analysis = self.sequence
                                self.flan_right_sequence_analysis = self.flan_sequence
                                self.left_sequence_analysis = "".join(
                                    [nucleotide_dict[n] for n in self.sequence[::-1]]
                                )
                                self.flan_left_sequence_analysis = "".join(
                                    [
                                        nucleotide_dict[n]
                                        for n in self.flan_sequence[::-1]
                                    ]
                                )

                            elif self.strand == -1:
                                genomic_end = int(self.gene_range[1])
                                genomic_flan_end = int(self.flan_gene_range[1])

                                if (
                                    self.sequence[genomic_end - (mutation_location + 1)]
                                    == altered_nuc
                                ):
                                    # Altered nucleotide in the given mutation fits with the sequence
                                    s = list(self.sequence)
                                    s[genomic_end - (mutation_location + 1)] = new_nuc
                                    self.sequence = "".join(s)
                                else:
                                    print(
                                        "\nGiven mutation location does not fit with the sequence."
                                        "\nNucleotides are different.\n"
                                    )
                                if (
                                    self.flan_sequence[
                                        genomic_flan_end - (mutation_location + 1)
                                    ]
                                    == altered_nuc
                                ):
                                    # Altered nucleotide in the given mutation fits with the sequence
                                    s = list(self.flan_sequence)
                                    s[genomic_flan_end - (mutation_location + 1)] = (
                                        new_nuc
                                    )
                                    self.flan_sequence = "".join(s)
                                else:
                                    print(
                                        "\nGiven mutation location does not fit with the sequence."
                                        "\nNucleotides are different.\n"
                                    )

                                self.left_sequence_analysis = self.sequence
                                self.flan_left_sequence_analysis = self.flan_sequence
                                self.right_sequence_analysis = "".join(
                                    [
                                        nucleotide_dict[n]
                                        for n in self.left_sequence_analysis[::-1]
                                    ]
                                )
                                self.flan_right_sequence_analysis = "".join(
                                    [
                                        nucleotide_dict[n]
                                        for n in self.flan_left_sequence_analysis[::-1]
                                    ]
                                )
        if self.sequence is not None:
            print("Sequence was retrieved.")

    def extract_gRNA_flan_sequence(
        self, location: int, direction: str, fivep: int, threep: int
    ) -> str:
        """
        Annotating flanking regions of the gRNAs

        :param location: Location of the gRNA target sequence
        :param direction: direction of the strand 'left' or 'right'
        :param fivep: Number of the nucleotides in the 5' flanking region
        :param threep: Number of the nucleotides in the 3' flanking region

        :return: flanking sequence
        """

        if direction == "left":
            grna_strand = "-1"
        else:
            grna_strand = "1"

        grna_flan_ensembl = (
            self.server
            + "/sequence/region/human/%s:%s?expand_3prime=%s;expand_5prime=%s;content-type=text/plain"
            % (location, grna_strand, threep, fivep)
        )
        grna_flan_request = requests.get(
            grna_flan_ensembl, headers={"Content-Type": "text/plain"}
        )

        if grna_flan_request.status_code != 200:
            print("No response from ensembl sequence!\n")

        return grna_flan_request.text

    def extract_info(
        self, chromosome: str, loc_start: int, loc_end: int, transcript: str | None
    ) -> int | None:
        """
        TODO documentation
        """
        if loc_start < loc_end:
            ensembl = (
                "/overlap/region/human/%s:%s-%s?feature=transcript;feature=exon;feature=mane;"
                "feature=cds" % (chromosome, str(loc_start), str(loc_end))
            )
        else:
            ensembl = (
                "/overlap/region/human/%s:%s-%s?feature=transcript;feature=exon;feature=mane;"
                "feature=cds" % (chromosome, str(loc_end), str(loc_start))
            )

        request = requests.get(
            self.server + ensembl, headers={"Content-Type": "application/json"}
        )

        info_dict = dict()

        if request.status_code != 200:
            print("No response from ensembl!")
        else:
            canonicals = list()
            refseq = ""
            for output in request.json():
                if transcript is None:
                    if (
                        output["feature_type"] == "mane"
                        and output["Parent"] == self.gene_id
                    ):
                        if "refseq_match" in output.keys():
                            if output["type"] != "MANE_Plus_Clinical":
                                if output["id"].split(".")[0] not in canonicals:
                                    canonicals.append(output["id"].split(".")[0])
                                    refseq = output["id"].split(".")[0]

                                if output["id"].split(".")[0] in canonicals:
                                    if output["id"] not in info_dict.keys():
                                        info_dict[output["id"]] = [
                                            {
                                                "start": output["start"],
                                                "end": output["end"],
                                            }
                                        ]

                                    else:
                                        old_val = info_dict[output["id"]]
                                        if {
                                            "start": output["start"],
                                            "end": output["end"],
                                        } not in old_val:
                                            old_val.append(
                                                {
                                                    "start": output["start"],
                                                    "end": output["end"],
                                                }
                                            )
                                            info_dict[output["id"]] = old_val

                    elif (
                        output["feature_type"] == "transcript"
                        and output["Parent"] == self.gene_id
                    ):
                        if "is_canonical" in output.keys():
                            if output["is_canonical"] == 1:
                                if output["id"].split(".")[0] not in canonicals:
                                    canonicals.append(output["id"].split(".")[0])
                        if "source" in output.keys():
                            if output["source"] == "ensembl_havana":
                                if output["id"].split(".")[0] not in canonicals:
                                    canonicals.append(output["id"].split(".")[0])
                        if output["id"].split(".")[0] in canonicals:
                            if output["id"] not in info_dict.keys():
                                info_dict[output["id"]] = [
                                    {"start": output["start"], "end": output["end"]}
                                ]

                            else:
                                old_val = info_dict[output["id"]]
                                if {
                                    "start": output["start"],
                                    "end": output["end"],
                                } not in old_val:
                                    old_val.append(
                                        {"start": output["start"], "end": output["end"]}
                                    )
                                    info_dict[output["id"]] = old_val
                else:
                    # Selected transcript
                    if (
                        output["Parent"] == self.gene_id
                        and output["id"].split(".")[0] == transcript
                    ):
                        transcript_info = "/lookup/id/%s?expand=1;mane=1" % output["id"]
                        transcript_request = requests.get(
                            self.server + transcript_info,
                            headers={"Content-Type": "application/json"},
                        )
                        if transcript_request.status_code != 200:
                            print("No response from ensembl for transcript id!")
                        else:
                            transcript_output = transcript_request.json()
                            if transcript not in info_dict.keys():
                                info_dict[transcript] = [
                                    {
                                        "start": transcript_output["start"],
                                        "end": transcript_output["end"],
                                    }
                                ]

                            else:
                                old_val = info_dict[transcript]
                                if {
                                    "start": transcript_output["start"],
                                    "end": transcript_output["end"],
                                } not in old_val:
                                    old_val.append(
                                        {
                                            "start": transcript_output["start"],
                                            "end": transcript_output["end"],
                                        }
                                    )
                                    info_dict[transcript] = old_val

            for output in request.json():
                if (
                    output["feature_type"] == "cds"
                    and info_dict != {}
                    and output["Parent"] in info_dict.keys()
                ):
                    for k in range(len(info_dict[output["Parent"]])):
                        d = info_dict[output["Parent"]][k]
                        coding_pos = list(range(output["start"], output["end"] + 1))
                        if "cds" not in d.keys():
                            d["cds"] = {output["protein_id"]: coding_pos}
                        else:
                            if output["protein_id"] not in d["cds"].keys():
                                d["cds"][output["protein_id"]] = coding_pos
                            else:
                                t = d["cds"][output["protein_id"]]
                                for i in coding_pos:
                                    if i not in t:
                                        t.append(i)
                                d["cds"][output["protein_id"]] = t
                        if d not in info_dict[output["Parent"]]:
                            info_dict[output["Parent"]][k] = d
            if refseq:
                selected_transcript = [refseq]
            else:
                protein_ids = list()
                swiss_protein_ids = list()
                selected_transcript = list()
                for ids in info_dict.keys():
                    for d in info_dict[ids]:
                        if "cds" in d.keys():
                            protein_ids.extend(d["cds"].keys())

                if protein_ids:
                    for p in protein_ids:
                        protein_ensembl = (
                            "/xrefs/id/{0}?external_db=Uniprot/SWISSPROT%".format(p)
                        )
                        protein_request = requests.get(
                            self.server + protein_ensembl,
                            headers={"Content-Type": "application/json"},
                        )
                        for i in protein_request.json():
                            if i["dbname"] == "Uniprot/SWISSPROT":
                                swiss_protein_ids.append(p)
                    if swiss_protein_ids:
                        for p in swiss_protein_ids:
                            for ids in info_dict.keys():
                                for d in info_dict[ids]:
                                    if "cds" in d.keys():
                                        if p in d["cds"].keys():
                                            if ids not in selected_transcript:
                                                selected_transcript.append(ids)

                    else:
                        selected_transcript = list(info_dict.keys())
                else:
                    selected_transcript = list(info_dict.keys())

            if selected_transcript:
                info_dict2 = {
                    key: info_dict[key]
                    for key in info_dict.keys()
                    if key in selected_transcript
                }

                for output in request.json():
                    if (
                        output["feature_type"] == "exon"
                        and info_dict2 != {}
                        and output["Parent"] in info_dict2.keys()
                    ):
                        for d in info_dict2[output["Parent"]]:
                            if "exon" not in d.keys():
                                d["exon"] = {
                                    output["exon_id"]: {
                                        "start": output["start"],
                                        "end": output["end"],
                                    }
                                }
                            else:
                                if output["exon_id"] not in d["exon"].keys():
                                    d["exon"][output["exon_id"]] = {
                                        "start": output["start"],
                                        "end": output["end"],
                                    }
                if info_dict2 != {}:
                    self.info_dict = info_dict2
                    return 1
            else:
                self.info_dict = None
                return 0

    def check_range_info(self, start: int, end: int) -> dict | None:
        """
        TODO documentation
        """
        range_locations = list(range(start, end))
        transcripts_exons = dict()
        if self.info_dict != {} and self.info_dict is not None:
            for transcript, transcript_dict_list in self.info_dict.items():
                for transcript_dict in transcript_dict_list:
                    in_transcript = False
                    for loc in range_locations:
                        if transcript_dict["start"] <= loc <= transcript_dict["end"]:
                            in_transcript = True
                    if in_transcript:
                        if "exon" in transcript_dict.keys():
                            for exon, exon_dict in transcript_dict["exon"].items():
                                in_exon = False
                                for loc in range_locations:
                                    if exon_dict["start"] <= loc <= exon_dict["end"]:
                                        in_exon = True
                                if in_exon:
                                    if transcript not in transcripts_exons.keys():
                                        transcripts_exons[transcript] = [exon]
                                    else:
                                        if exon not in transcripts_exons[transcript]:
                                            transcripts_exons[transcript].append(exon)
                        else:
                            if transcript not in transcripts_exons.keys():
                                transcripts_exons[transcript] = list()
        if transcripts_exons != {}:
            return transcripts_exons
        else:
            return None

    def check_cds(self, transcript_id: str, start: int, end: int) -> bool:
        """
        TODO documentation
        """
        range_locations = list(range(start, end))
        in_cds = False
        if self.info_dict != {} and self.info_dict is not None:
            if transcript_id in self.info_dict.keys():
                for transcript_dict in self.info_dict[transcript_id]:
                    if "cds" in transcript_dict.keys():
                        for protein, protein_pos in transcript_dict["cds"].items():
                            for loc in range_locations:
                                if loc in protein_pos:
                                    in_cds = True

        return in_cds

    def extract_uniprot_info(self, ensembl_pid, uniprot):
        """
        TODO documentation
        """
        uniprot = uniprot.split(".")[0]
        protein_ensembl = "/xrefs/id/{0}?external_db=Uniprot/SWISSPROT%".format(
            ensembl_pid
        )
        protein_request = requests.get(
            self.server + protein_ensembl, headers={"Content-Type": "application/json"}
        )
        if protein_request.status_code != 200:
            print("No response from ensembl!")
            return 0
        else:
            seq_mapping = dict()
            for i in protein_request.json():
                if i["primary_id"] == uniprot:
                    if i["dbname"] == "Uniprot/SWISSPROT":
                        if (
                            "ensembl_end" in i.keys()
                            and "ensembl_start" in i.keys()
                            and "xref_end" in i.keys()
                            and "xref_start" in i.keys()
                        ):
                            if int(i["ensembl_end"]) - int(i["ensembl_start"]) == int(
                                i["xref_end"]
                            ) - int(i["xref_start"]):
                                # Otherwise, there is an inconsistency --> Not take it
                                seq_mapping[uniprot] = {
                                    i["ensembl_start"] + k: i["xref_start"] + k
                                    for k in range(
                                        int(i["ensembl_end"])
                                        - int(i["ensembl_start"])
                                        + 1
                                    )
                                }
                            break

                    elif i["dbname"] == "Uniprot/SPTREMBL":
                        if (
                            "ensembl_end" in i.keys()
                            and "ensembl_start" in i.keys()
                            and "xref_end" in i.keys()
                            and "xref_start" in i.keys()
                        ):
                            if int(i["ensembl_end"]) - int(i["ensembl_start"]) == int(
                                i["xref_end"]
                            ) - int(i["xref_start"]):
                                # Otherwise, there is an inconsistency --> Not take it
                                seq_mapping[uniprot] = {
                                    i["ensembl_start"] + k: i["xref_start"] + k
                                    for k in range(
                                        int(i["ensembl_end"])
                                        - int(i["ensembl_start"])
                                        + 1
                                    )
                                }

        if seq_mapping == dict() or seq_mapping == 0:
            # Alignment
            uniprot_obj = Uniprot(uniprotid=uniprot)
            uniprot_obj.extract_uniprot()
            uniprot_seq = uniprot_obj.sequence

            if self.p_sequence is not None:
                ensembl_seq = self.p_sequence
            else:
                seq_ensembl = self.server + "/sequence/id/%s?" % ensembl_pid
                seq_request = requests.get(
                    seq_ensembl, headers={"Content-Type": "text/x-fasta"}
                )

                if (
                    seq_request.status_code != 200
                    and seq_flan_request.status_code != 200
                ):
                    print("No response from ensembl protein sequence!\n")

                else:
                    label_line = seq_request.text.split("\n")[0]
                if label_line[0] != "{":
                    self.p_sequence = "".join(seq_request.text.split("\n")[1:])
                    ensembl_seq = self.p_sequence

            # TODO: This can't be correct?
            # OUTPUT_PATH + args[constants.ARGS_KEY_OUTPUT_PATH] === OUTPUT_PATH + OUTPUT_PATH
            # The string concatenation seems wrong
            alignment_f = open(
                OUTPUT_PATH
                + CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH]
                + "_%s_alignment.txt" % uniprot,
                "w",
            )
            alignments = list()
            for a in pairwise2.align.globalms(
                ensembl_seq, uniprot_seq, 2, -1, -1, -0.1
            ):
                alignment_f.write(format_alignment(*a, full_sequences=True))
                alignments.append(format_alignment(*a, full_sequences=True))
            alignment_f.close()

            ens = alignments[0].split("\n")[:3]
            alignment_df = pandas.DataFrame(
                0, index=list(range(len(ens[0]))), columns=["e_aa", "e_index"]
            )
            alignment_df["e_aa"] = [aa for aa in ens[0]]

            alignment_dict = {i: alignments[i] for i in range(len(alignments))}
            for a_num, align in alignment_dict.items():
                alignment_df["u_aa_%d" % a_num] = None
                alignment_df["u_index_%d" % a_num] = None
                alignment_df["mismatch_%d" % a_num] = None
                align_list = align.split("\n")[:3]
                uni_index = 1
                ens_index = 1
                for i in range(len(align_list[0])):
                    ens = align_list[0][i]
                    a_style = align_list[1][i]
                    uni = align_list[2][i]

                    if a_style == "|":
                        # Match
                        alignment_df.loc[i, "e_index"] = ens_index
                        alignment_df.loc[i, "u_aa_%d" % a_num] = uni
                        alignment_df.loc[i, "u_index_%d" % a_num] = uni_index
                        alignment_df.loc[i, "mismatch_%d" % a_num] = False
                        uni_index += 1
                        ens_index += 1
                    elif a_style == ".":
                        # Mismatch
                        alignment_df.loc[i, "e_index"] = ens_index
                        alignment_df.loc[i, "u_aa_%d" % a_num] = uni
                        alignment_df.loc[i, "u_index_%d" % a_num] = uni_index
                        alignment_df.loc[i, "mismatch_%d" % a_num] = True
                        uni_index += 1
                        ens_index += 1
                    elif a_style == " ":
                        # Gap
                        if ens == "-":
                            alignment_df.loc[i, "e_index"] = ens_index
                            alignment_df.loc[i, "u_aa_%d" % a_num] = uni
                            alignment_df.loc[i, "u_index_%d" % a_num] = uni_index
                            alignment_df.loc[i, "mismatch_%d" % a_num] = False
                            uni_index += 1
                        elif uni == "-":
                            alignment_df.loc[i, "e_index"] = ens_index
                            alignment_df.loc[i, "u_aa_%d" % a_num] = "-"
                            alignment_df.loc[i, "u_index_%d" % a_num] = uni_index
                            alignment_df.loc[i, "mismatch_%d" % a_num] = False
                            ens_index += 1

            alignment_df["inconsistent"] = False
            mismatch_indices = list()
            for i in range(len(alignments)):
                mismatch_indices.extend(
                    list(alignment_df[alignment_df["mismatch_%d" % a_num]].index)
                )
            mismatch_indices = list(set(mismatch_indices))
            different_indices = list()
            missing_indices = list()
            uni_ind_cols = [col for col in alignment_df if col.startswith("u_index")]
            uni_aa_cols = [col for col in alignment_df if col.startswith("u_aa")]
            for ind, row in alignment_df.iterrows():
                indices = list()
                aas = list()
                for col in uni_ind_cols:
                    indices.append(row[col])
                for col in uni_aa_cols:
                    aas.append(row[col])
                    if row[col] == "-":
                        if ind not in missing_indices:
                            missing_indices.append(ind)
                if len(list(set(indices))) > 1:
                    different_indices.append(ind)
                if len(list(set(aas))) > 1:
                    if ind not in different_indices:
                        different_indices.append(ind)

            flagged_ind = list(
                set(mismatch_indices).union(
                    set(different_indices).union(set(missing_indices))
                )
            )
            alignment_df.loc[flagged_ind, "inconsistent"] = True

            # TODO: This can't be correct?
            # OUTPUT_PATH + args[constants.ARGS_KEY_OUTPUT_PATH] === OUTPUT_PATH + OUTPUT_PATH
            # The string concatenation seems wrong
            alignment_df.to_csv(
                OUTPUT_PATH
                + CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH]
                + "_%s_alignment_df.csv" % uniprot,
                index=True,
            )

            alignment_df2 = alignment_df[alignment_df.inconsistent == False]
            seq_mapping[uniprot] = {
                r["e_index"]: r["u_index_0"] for i, r in alignment_df2.iterrows()
            }

        if seq_mapping:
            return seq_mapping
        else:
            return None
