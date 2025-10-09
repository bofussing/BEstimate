# Copyright (C) 2025 Genome Research Ltd.

import re
import requests

class Uniprot:
    """
    A class for interacting with the UniProt API to retrieve protein information.

    This class handles UniProt protein data retrieval including sequence information,
    protein domains, post-translational modifications (PTMs), and mutagenesis data.
    It provides methods to extract and analyze protein features for base editor analysis.

    :param uniprotid: UniProt accession identifier for the protein of interest
    :type uniprotid: str

    :ivar uniprotid: UniProt accession identifier
    :ivar reviewed: Whether the UniProt entry is reviewed (SwissProt) or unreviewed (TrEMBL)
    :ivar sequence: Protein amino acid sequence
    :ivar domains: Dictionary of protein domains and their positions
    :ivar phosphorylation_sites: Dictionary of phosphorylation sites and descriptions
    :ivar ubiquitination_sites: Dictionary of ubiquitination sites and descriptions
    :ivar methylation_sites: Dictionary of methylation sites and descriptions
    :ivar acetylation_sites: Dictionary of acetylation sites and descriptions
    :ivar mutagenesis: Dictionary of mutagenesis data from UniProt
    :ivar server: Base URL for UniProt API
    """

    def __init__(self, uniprotid: str):
        """
        Initialize UniProt object with protein identifier.

        :param uniprotid: UniProt accession identifier
        """
        self.uniprotid, self.reviewed = uniprotid, None
        self.sequence = None
        self.domains = dict()
        self.phosphorylation_sites = dict()
        self.ubiquitination_sites = dict()
        self.methylation_sites = dict()
        self.acetylation_sites = dict()
        self.mutagenesis = dict()
        self.server = "https://www.ebi.ac.uk/proteins/api/"

    def extract_uniprot(self) -> str:
        """
        Extract protein information from UniProt API.

        Retrieves protein sequence, domains, and post-translational modifications
        from the UniProt database using the API. Processes features including
        domains, binding sites, and various PTMs (phosphorylation, methylation,
        ubiquitination, acetylation).

        :return: Status message indicating completion of API request

        .. note::
            This method populates the object's attributes with data from UniProt.
            If no data is found for certain categories, the corresponding attributes
            are set to None.
        """

        uniprot_api = "proteins?offset=0&size=-1&accession=%s" % self.uniprotid
        api_request = requests.get(
            self.server + uniprot_api, headers={"Accept": "application/json"}
        )

        # Check the response of the server for the request
        if api_request.status_code != 200:
            return "No response from UniProt!\n"

        else:
            for i in api_request.json():
                self.reviewed = False if i["info"]["type"] == "TrEMBL" else True
                self.sequence = i["sequence"]["sequence"]
                if "features" in i.keys() and i["features"] != []:
                    for ftr in i["features"]:
                        if ftr["type"] == "MOD_RES" and ftr["category"] == "PTM":
                            if "description" in ftr.keys():
                                pos, ptm = ftr["begin"], ftr["description"]
                                ptm = ptm.split(";")[0]
                                # Phosphorylation
                                phos = (
                                    ptm
                                    if re.search(r"Phospho", ptm)
                                    or re.search(r"phospho", ptm)
                                    else None
                                )
                                if phos is not None:
                                    self.phosphorylation_sites[pos] = phos
                                # Methylation
                                methy = (
                                    ptm
                                    if re.search(r"Methyl", ptm)
                                    or re.search(r"methyl", ptm)
                                    else None
                                )
                                if methy is not None:
                                    self.methylation_sites[pos] = methy
                                # Ubiquitination
                                ubi = (
                                    ptm
                                    if re.search(r"Ub", ptm) or re.search(r"ub", ptm)
                                    else None
                                )
                                if ubi is not None:
                                    self.ubiquitination_sites[pos] = ubi
                                # Acetylation
                                acety = (
                                    ptm
                                    if re.search(r"Ace", ptm)
                                    or re.search(r"acetyl", ptm)
                                    else None
                                )
                                if acety is not None:
                                    self.acetylation_sites[pos] = acety

                        if ftr["category"] == "DOMAINS_AND_SITES":
                            if ftr["type"] == "BINDING":
                                domain = ftr["ligand"]["name"] + " binding site"
                            else:
                                if "description" in ftr.keys():
                                    domain = ftr["description"]
                            domain_range = list(
                                range(int(ftr["begin"]), int(ftr["end"]))
                            )
                            if domain not in self.domains.keys():
                                self.domains[domain] = domain_range
                            else:
                                t = self.domains[domain]
                                for p in domain_range:
                                    if p not in t:
                                        t.append(p)
                                self.domains[domain] = t

        if self.phosphorylation_sites == dict():
            self.phosphorylation_sites = None
        if self.methylation_sites == dict():
            self.methylation_sites = None
        if self.acetylation_sites == dict():
            self.acetylation_sites = None
        if self.ubiquitination_sites == dict():
            self.ubiquitination_sites = None
        if self.domains == dict():
            self.domains = None

        return "UniProt API request is done."

    def find_domain(self, protein_edit_location: int, old_aa: str) -> str | None:
        """
        Find protein domain at a specific amino acid position.

        Checks if the given protein position falls within any known protein domains
        and validates that the amino acid at that position matches the expected residue.

        :param protein_edit_location: Position in the protein sequence (1-indexed)
        :param old_aa: Expected amino acid at the given position

        :return: Domain name if position is within a domain and amino acid matches, None otherwise
        """

        edit_domain = None
        if self.domains != {} and self.domains is not None:
            for domain, domain_range in self.domains.items():
                if int(protein_edit_location) in domain_range:
                    if old_aa == self.sequence[protein_edit_location - 1]:
                        edit_domain = domain
        return edit_domain

    def find_ptm_site(self, ptm_type, protein_edit_location, old_aa):
        """
        Find post-translational modification site at a specific position.

        Checks if the given protein position corresponds to a known PTM site
        of the specified type and validates the amino acid at that position.

        :param ptm_type: Type of PTM to search for ('phosphorylation', 'methylation', 'ubiquitination', 'acetylation')
        :type ptm_type: str
        :param protein_edit_location: Position in the protein sequence (1-indexed)
        :type protein_edit_location: int
        :param old_aa: Expected amino acid at the given position
        :type old_aa: str

        :return: PTM description if found at the position, None otherwise
        :rtype: str or None
        """

        edit_ptm_site = None
        if ptm_type == "phosphorylation":
            d = self.phosphorylation_sites
        if ptm_type == "methylation":
            d = self.methylation_sites
        if ptm_type == "ubiquitination":
            d = self.ubiquitination_sites
        if ptm_type == "acetylation":
            d = self.acetylation_sites

        if d != {} and d is not None:
            for ptm_pos, ptm in d.items():
                if int(ptm_pos) == int(protein_edit_location):
                    if old_aa == self.sequence[int(protein_edit_location) - 1]:
                        edit_ptm_site = ptm
        return edit_ptm_site

    def extract_mutagenesis(self) -> str | None:
        """
        Extract mutagenesis data from UniProt API.

        Retrieves experimental mutagenesis data for the protein, including
        amino acid substitutions and their phenotypic effects.

        :return: Status message indicating completion of mutagenesis data extraction

        .. note::
            This method populates the mutagenesis attribute with experimental
            data from UniProt's mutagenesis features.
        """

        mut_api = "features/%s?types=MUTAGEN" % self.uniprotid
        mut_api_request = requests.get(
            self.server + mut_api, headers={"Accept": "application/json"}
        )

        # Check the response of the server for the request
        if mut_api_request.status_code != 200:
            return "No response from UniProt Mutagen!\n"

        else:
            if "features" in mut_api_request.json().keys():
                if mut_api_request.json()["features"]:
                    for mut in mut_api_request.json()["features"]:
                        if mut["category"] == "MUTAGENESIS":
                            if mut["begin"] == mut["end"]:
                                self.mutagenesis[mut["begin"]] = {
                                    mut["alternativeSequence"]: [mut["description"]]
                                }
                                if (
                                    mut["alternativeSequence"]
                                    not in self.mutagenesis[mut["begin"]].keys()
                                ):
                                    self.mutagenesis[mut["begin"]][
                                        mut["alternativeSequence"]
                                    ] = [mut["description"]]
                                else:
                                    if (
                                        mut["description"]
                                        not in self.mutagenesis[mut["begin"]][
                                            mut["alternativeSequence"]
                                        ]
                                    ):
                                        self.mutagenesis[mut["begin"]][
                                            mut["alternativeSequence"]
                                        ].append(mut["description"])

    def find_mutagenesis(self, protein_edit_location: int, new_aa: str) -> str | None:
        """
        Find mutagenesis data for a specific amino acid substitution.

        Searches for experimental mutagenesis data at the given position with
        the specified amino acid substitution.

        :param protein_edit_location: Position in the protein sequence (1-indexed)
        :param new_aa: New amino acid after substitution

        :return: Mutagenesis description(s) if found, None otherwise
        """

        if str(protein_edit_location) in self.mutagenesis.keys():
            mut = self.mutagenesis[str(protein_edit_location)]
            if new_aa in mut.keys():
                mutagenesis = mut[new_aa]
                if len(mutagenesis) == 1:
                    return mutagenesis[0]
                else:
                    return ";".join(mutagenesis)

