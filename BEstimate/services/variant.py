# Copyright (C) 2025 Genome Research Ltd.

import pandas
import re


class Variant:
    """
    A class for handling variant information and VEP (Variant Effect Predictor) analysis.

    This class processes genomic variants and their consequences using Ensembl's VEP API.
    It handles HGVS nomenclature, consequence prediction, and clinical significance assessment.

    :ivar hgvs: Original HGVS string
    :ivar hgvsc: HGVS coding sequence nomenclature
    :ivar hgvsp: HGVS protein nomenclature
    :ivar vep: VEP response data
    :ivar gene: Gene symbol
    :ivar strand: Strand orientation
    :ivar transcript: Transcript identifier
    :ivar variant_class: Type of variant (SNV, insertion, etc.)
    :ivar consequence_terms: List of consequence terms
    :ivar most_severe_consequence: Most severe predicted consequence
    :ivar protein_change: Protein-level change description
    :ivar old_aa: Original amino acid
    :ivar new_aa: Changed amino acid
    :ivar clinical: Clinical significance information
    """

    def __init__(self, hgvs: str, gene: str, strand: int, transcript: str) -> None:
        """
        Initialize Variant object with genomic variant information.

        :param hgvs: HGVS nomenclature string
        :param gene: Gene symbol where the variant is located
        :param strand: Strand orientation (1 or -1)
        :param transcript: Ensembl transcript identifier
        """
        self.hgvs, self.hgvsc, self.hgvsp = hgvs, None, None
        self.vep = None
        self.gene, self.strand = gene, strand
        self.transcript = transcript
        self.allele = None
        self.regulatory = None
        self.motif, self.motif_TFs = None, None
        self.variant_class, self.consequence_terms, self.biotype = None, None, None
        self.most_severe_consequence = None
        self.cdna_change, self.cds_position = None, None
        self.old_codon, self.new_codon = None, None
        self.protein_change, self.protein_position = None, None
        self.old_aa, self.new_aa = None, None
        self.old_aa_chem, self.new_aa_chem = None, None
        self.synonymous, self.stop, self.proline = None, None, None
        self.protein, self.swissprot = None, None
        self.polyphen_score, self.polyphen_prediction = None, None
        self.sift_score, self.sift_prediction = None, None
        self.cadd_phred, self.cadd_raw, self.lof = None, None, None
        self.impact, self.blosum62 = None, None
        self.clinical, self.clinical_id, self.clinical_sig = None, None, None
        self.clinvar_id = None
        self.cosmic, self.cosmic_id = None, None
        self.ancestral_populations = None

    def extract_vep_obj(self, vep_json: list) -> None:
        """
        Extract VEP response data for this variant from batch response.

        Searches through the VEP batch response JSON to find the data
        corresponding to this variant's HGVS string and stores it in
        the vep attribute.

        :param vep_json: VEP API response containing multiple variant results
        """
        for vep in vep_json:
            if vep["input"] == self.hgvs:
                self.vep = vep

    def extract_hgvsp(self, hgvsp: str, which: str) -> str | None:
        """
        Parse HGVS protein nomenclature to extract specific components.

        Parses protein-level HGVS nomenclature to extract amino acid information
        including original amino acid, new amino acid, and position.

        :param hgvsp: HGVS protein nomenclature string
        :param which: Component to extract ('old_aa', 'new_aa', or 'position')

        :return: Requested component from HGVS protein nomenclature
        """
        aa_3to1 = {
            "Ala": "A",
            "Arg": "R",
            "Asn": "N",
            "Asp": "D",
            "Cys": "C",
            "Glu": "E",
            "Gln": "Q",
            "Gly": "G",
            "His": "H",
            "Ile": "I",
            "Leu": "L",
            "Lys": "K",
            "Met": "M",
            "Phe": "F",
            "Pro": "P",
            "Ser": "S",
            "Thr": "T",
            "Trp": "W",
            "Tyr": "Y",
            "Val": "V",
            "Ter": "*",
        }
        if hgvsp is not None:
            protein_change = hgvsp.split("p.")[1]
            if len(protein_change.split("delins")) == 1:
                # SNP
                if len(protein_change.split("=")) == 1:
                    if len(protein_change.split("?")) == 1:
                        if len(protein_change.split("ext")) == 1:
                            if which == "old_aa":
                                return aa_3to1[protein_change[:3]]
                            if which == "new_aa":
                                return aa_3to1[protein_change[-3:]]
                            if which == "position":
                                return protein_change[3:-3]
                        else:
                            # Extension for termination or start Ter629GlnextTer1 | Met1ext-5
                            if protein_change[:3] == "Ter":
                                alteration = protein_change.split("ext")[0]
                                extension_amount = (
                                    int(protein_change.split("ext")[1][3:]) - 1
                                )
                                if which == "old_aa":
                                    return aa_3to1[alteration[:3]]
                                if which == "new_aa":
                                    return (
                                        aa_3to1[alteration[-3:]]
                                        + "X%s" % extension_amount
                                        + "*"
                                    )
                                if which == "position":
                                    return alteration[3:-3]
                            else:
                                if which == "old_aa":
                                    return aa_3to1[protein_change[:3]]
                                if which == "new_aa":
                                    extension_amount = (
                                        abs(int(protein_change.split("ext")[1])) - 1
                                    )
                                    return (
                                        aa_3to1[protein_change[:3]]
                                        + "X-%s" % extension_amount
                                        + aa_3to1[protein_change[:3]]
                                    )
                                if which == "position":
                                    return protein_change.split("ext")[0][3:]

                    else:
                        # Start codon lost - Met1? | MetAla1_?2
                        if which == "old_aa":
                            aa1 = list()
                            aa_string = re.match(
                                "([a-z]+)([0-9]+)", protein_change.split("?")[0], re.I
                            ).groups()[0]
                            for i in [
                                aa_string[x : x + 3]
                                for x in range(0, len(aa_string), 3)
                            ]:
                                aa1.append(aa_3to1[i])
                            return ";".join(aa1)

                        if which == "new_aa":
                            if protein_change[-1] == "?" or protein_change[-2] == "?":
                                return "-"
                        if which == "position":
                            return re.match(
                                "([a-z]+)([0-9]+)", protein_change.split("?")[0], re.I
                            ).groups()[1]

                else:
                    if which == "old_aa":
                        # Synonymous variant
                        aa1 = list()
                        aa_string = re.match(
                            "([a-z]+)([0-9]+)", protein_change.split("=")[0], re.I
                        ).groups()[0]
                        for i in [
                            aa_string[x : x + 3] for x in range(0, len(aa_string), 3)
                        ]:
                            aa1.append(aa_3to1[i])
                        return ";".join(aa1)
                    if which == "new_aa":
                        # Synonymous variant
                        aa1 = list()
                        aa_string = re.match(
                            "([a-z]+)([0-9]+)", protein_change.split("=")[0], re.I
                        ).groups()[0]
                        for i in [
                            aa_string[x : x + 3] for x in range(0, len(aa_string), 3)
                        ]:
                            aa1.append(aa_3to1[i])
                        return ";".join(aa1)
                    if which == "position":
                        return re.match(
                            "([a-z]+)([0-9]+)", protein_change.split("=")[0], re.I
                        ).groups()[1]

            elif len(protein_change.split("delins")) > 1:
                # Substitution
                if which == "old_aa":
                    aa1 = list()
                    for i in protein_change.split("delins")[0].split("_"):
                        aa1.append(aa_3to1[i[:3]])
                    return ";".join(aa1)

                if which == "new_aa":
                    aa2 = list()
                    for i in [
                        protein_change.split("delins")[1][x : x + 3]
                        for x in range(0, len(protein_change.split("delins")[1]), 3)
                    ]:
                        aa2.append(aa_3to1[i])
                    return ";".join(aa2)

                if which == "position":
                    pos = list()
                    for i in protein_change.split("delins")[0].split("_"):
                        pos.append(re.match("([a-z]+)([0-9]+)", i, re.I).groups()[1])
                    return ";".join(pos)

        else:
            return None

    def extract_consequences(self) -> None:
        """
        Extract and process variant consequences from VEP response data.

        Parses the VEP (Variant Effect Predictor) response to extract variant
        consequences, protein effects, clinical significance, and other annotations.
        Sets multiple object attributes with processed consequence information.

        .. note::
            This method processes regulatory features, motif features, transcript
            consequences, and clinical annotations from the VEP response. It also
            determines amino acid chemical property changes and clinical significance.
        """
        consequence_terms = list()
        ancestral_populations = list()
        # Dictionary to find the chemical properperty change due to the edit
        aa_chem = {
            "G": "Non-Polar",
            "A": "Non-Polar",
            "V": "Non-Polar",
            "C": "Polar",
            "P": "Non-Polar",
            "L": "Non-Polar",
            "I": "Non-Polar",
            "M": "Non-Polar",
            "W": "Non-Polar",
            "F": "Non-Polar",
            "S": "Polar",
            "T": "Polar",
            "Y": "Polar",
            "N": "Polar",
            "Q": "Polar",
            "K": "Charged",
            "R": "Charged",
            "H": "Charged",
            "D": "Charged",
            "E": "Charged",
            "*": "-",
        }

        if "allele_string" in self.vep.keys():
            self.allele = self.vep["allele_string"]

        if "most_severe_consequence" in self.vep.keys():
            self.most_severe_consequence = self.vep["most_severe_consequence"]

        if "variant_class" in self.vep.keys():
            self.variant_class = self.vep["variant_class"]

        if "regulatory_feature_consequences" in self.vep.keys():
            for r in self.vep["regulatory_feature_consequences"]:
                if "strand" in r.keys():
                    if r["strand"] == self.strand:
                        if "regulatory_feature_id" in r.keys():
                            self.regulatory = r["regulatory_feature_id"]
                        if "consequence_terms" in r.keys():
                            for cons_term in r["consequence_terms"]:
                                if cons_term not in consequence_terms:
                                    consequence_terms.append(cons_term)

        if "motif_feature_consequences" in self.vep.keys():
            for m in self.vep["motif_feature_consequences"]:
                if "motif_feature_id" in m.keys():
                    self.motif = m["motif_feature_id"]
                if "transcription_factors" in m.keys():
                    self.motif_TFs = ", ".join(
                        [tf for tf in m["transcription_factors"]]
                    )
                if "consequence_terms" in m.keys():
                    for cons_term in m["consequence_terms"]:
                        if cons_term not in consequence_terms:
                            consequence_terms.append(cons_term)

        if "transcript_consequences" in self.vep.keys():
            for t in self.vep["transcript_consequences"]:
                if (
                    t["gene_symbol"] == self.gene
                    and t["transcript_id"] == self.transcript
                ):
                    if "hgvsc" in t.keys():
                        self.hgvsc = t["hgvsc"]
                    if "biotype" in t.keys():
                        self.biotype = t["biotype"]
                    if "hgvsp" in t.keys():
                        self.hgvsp = t["hgvsp"]
                        self.protein_position = self.extract_hgvsp(
                            hgvsp=self.hgvsp, which="position"
                        )
                        self.old_aa = self.extract_hgvsp(
                            hgvsp=self.hgvsp, which="old_aa"
                        )
                        self.new_aa = self.extract_hgvsp(
                            hgvsp=self.hgvsp, which="new_aa"
                        )
                        self.old_aa_chem = (
                            aa_chem[self.old_aa]
                            if self.old_aa is not None
                            and self.old_aa in aa_chem.keys()
                            and len(self.old_aa) == 1
                            else (
                                ";".join(
                                    [
                                        aa_chem[i]
                                        for i in self.old_aa.split(";")
                                        if i in aa_chem.keys()
                                    ]
                                )
                                if self.old_aa is not None and len(self.old_aa) > 1
                                else None
                            )
                        )
                        self.new_aa_chem = (
                            aa_chem[self.new_aa]
                            if self.new_aa is not None
                            and self.new_aa in aa_chem.keys()
                            and len(self.new_aa) == 1
                            else (
                                ";".join(
                                    [
                                        aa_chem[i]
                                        for i in self.new_aa.split(";")
                                        if i in aa_chem.keys()
                                    ]
                                )
                                if self.new_aa is not None and len(self.new_aa) > 1
                                else None
                            )
                        )

                    if "protein_id" in t.keys():
                        self.protein = t["protein_id"]
                    if "amino_acids" in t.keys():
                        self.protein_change = t["amino_acids"]
                    if "codons" in t.keys():
                        self.cdna_change = t["codons"]
                        self.old_codon = (
                            self.cdna_change.split("/")[0]
                            if self.cdna_change is not None
                            and pandas.isna(self.cdna_change) is False
                            and type(self.cdna_change) != float
                            else None
                        )
                        self.new_codon = (
                            self.cdna_change.split("/")[1]
                            if self.cdna_change is not None
                            and pandas.isna(self.cdna_change) is False
                            and type(self.cdna_change) != float
                            else None
                        )
                    if "cds_start" in t.keys() and "cds_end" in t.keys():
                        self.cds_position = (
                            str(t["cds_start"]) + "-" + str(t["cds_end"])
                        )

                    if self.cdna_change and self.protein_change:
                        self.synonymous = (
                            True
                            if self.old_codon is not None
                            and self.new_codon is not None
                            and self.old_aa is not None
                            and self.new_aa is not None
                            and self.old_codon != self.new_codon
                            and self.old_aa == self.new_aa
                            else (
                                None
                                if self.old_codon is None
                                and self.new_codon is None
                                or self.old_aa is None
                                and self.new_aa is None
                                else False
                            )
                        )
                        self.proline = (
                            True
                            if self.synonymous is not None
                            and self.synonymous == False
                            and self.new_aa is not None
                            and "P" in self.new_aa.split(";")
                            else False
                        )
                        self.stop = (
                            True
                            if self.new_aa is not None
                            and self.new_aa == "*"
                            and len(self.new_aa) == 1
                            else (
                                True
                                if self.new_aa is not None
                                and "*" in self.new_aa
                                and len(self.new_aa) > 1
                                else False
                            )
                        )

                    if "swissprot" in t.keys():
                        self.swissprot = t["swissprot"][0].split(".")[0]
                    if "polyphen_score" in t.keys():
                        self.polyphen_score = t["polyphen_score"]
                    if "polyphen_prediction" in t.keys():
                        self.polyphen_prediction = t["polyphen_prediction"]
                    if "sift_score" in t.keys():
                        self.sift_score = t["sift_score"]
                    if "sift_prediction" in t.keys():
                        self.sift_prediction = t["sift_prediction"]
                    if "cadd_phred" in t.keys():
                        self.cadd_phred = t["cadd_phred"]
                    if "cadd_raw" in t.keys():
                        self.cadd_raw = t["cadd_raw"]
                    if "lof" in t.keys():
                        self.lof = t["lof"]
                    if "impact" in t.keys():
                        self.impact = t["impact"]
                    if "blosum62" in t.keys():
                        self.blosum62 = t["blosum62"]
                    if "consequence_terms" in t.keys():
                        for cons_term in t["consequence_terms"]:
                            if cons_term not in consequence_terms:
                                consequence_terms.append(cons_term)

        if consequence_terms:
            self.consequence_terms = ", ".join(consequence_terms)

        if "colocated_variants" in self.vep.keys():
            self.clinical = True
            cosmic_id = list()
            clinvar_id = list()
            for c in self.vep["colocated_variants"]:
                if "allele_string" in c.keys():
                    if c["allele_string"] == "COSMIC_MUTATION":
                        self.cosmic = True
                        if "id" in c.keys():
                            if c["id"] not in cosmic_id:
                                cosmic_id.append(c["id"])

                if "clin_sig" in c.keys():
                    self.clinical_sig = ", ".join([i for i in c["clin_sig"]])
                if "id" in c.keys():
                    self.clinical_id = c["id"]

                if "var_synonyms" in c.keys():
                    if type(c["var_synonyms"]) == str:
                        for clnv in c["var_synonyms"]["ClinVar"]:
                            for cl_id in clnv:
                                if cl_id not in clinvar_id:
                                    clinvar_id.append(cl_id)
                if "frequencies" in c.keys():
                    for alele, freq_dict in c["frequencies"].items():
                        for pop, val in freq_dict.items():
                            if val >= 0.01:
                                if pop not in ancestral_populations:
                                    ancestral_populations.append(pop)

            if cosmic_id:
                self.cosmic_id = ", ".join(cosmic_id)
            if clinvar_id:
                self.clinvar_id = ", ".join(clinvar_id)
            if ancestral_populations:
                self.ancestral_populations = ", ".join(ancestral_populations)