# -----------------------------------------------------------------------------------------#
#                                                                                          #
#                                  B E s t i m a t e                                       #
#                        Author : Cansu Dincer cd7@sanger.ac.uk                            #
#                         Dr Matthew Coelho & Dr Mathew Garnett                            #
#                              Wellcome Sanger Institute                                   #
#                                                                                          #
# -----------------------------------------------------------------------------------------#


import typing as t
import os, sys, pandas, argparse, requests, json
import warnings
from Bio.pairwise2 import format_alignment

from BEstimate.datafiles import DataFiles
from BEstimate import constants
from BEstimate import x_crispranalyser
import BEstimate

from BEstimate.services.uniprot import Uniprot
from BEstimate.services.ensembl import Ensembl

import BEstimate.sequences.utils as seq_utils

# GLOBAL VARIABLES
OT_PATH: str = ""
OUTPUT_PATH: str = ""
CLI_ARGS: dict[str, t.Any] = {}


def take_input() -> dict[str, t.Any]:
    """
    Parse command line arguments for BEstimate base editor analysis.

    This function sets up an argument parser to handle all the input parameters
    required for base editor site analysis including gene information, PAM sequences,
    activity windows, and output options.

    :return: Dictionary containing all parsed command line arguments

    .. note::
        This function configures and parses command line arguments using argparse.
        Default values are provided for most optional parameters.
    """
    parser = argparse.ArgumentParser(
        prog=constants.PROGRAM_NAME,
        usage="%(prog)s [inputs]",
        description="""
                                     **********************************
                                     Find and Analyse Base Editor sites
                                     **********************************""",
    )

    for group in parser._action_groups:
        if group.title == "optional arguments":
            group.title = "Inputs"
        elif "positional arguments":
            group.title = "Mandatory Inputs"

    # BASIC INFORMATION
    version_str = f"{constants.PROGRAM_NAME} {BEstimate.__version__}"
    parser.add_argument(
        "--version",
        action="version",
        version=version_str,
        help="Show program's version number and exit.",
    )
    parser.add_argument(
        "-gene",
        dest="GENE",
        required=True,
        help="The hugo symbol of the interested gene!",
    )

    parser.add_argument(
        "-assembly",
        dest="ASSEMBLY",
        required=True,
        default="GRCh38",
        help="The genome assembly that will be used!",
    )

    parser.add_argument(
        "-transcript",
        dest="TRANSCRIPT",
        default=None,
        help="The interested ensembl transcript id",
    )

    parser.add_argument(
        "-uniprot", dest="UNIPROT", default=None, help="The interested Uniprot id"
    )

    # PAM AND PROTOSPACER INFORMATION

    # The NGG PAM will be used unless otherwise specified.
    parser.add_argument(
        "-pamseq",
        dest="PAMSEQ",
        default="NGG",
        help="The PAM sequence in which features used "
        "for searching activity window and editable nucleotide.",
    )
    parser.add_argument(
        "-pamwin",
        dest="PAMWINDOW",
        default="21-23",
        help="The index of the PAM sequence when starting "
        "from the first index of protospacer as 1.",
    )
    parser.add_argument(
        "-actwin",
        dest="ACTWINDOW",
        default="4-8",
        help="The index of the activity window when starting "
        "from the first index of protospacer as 1.",
    )
    parser.add_argument(
        "-protolen",
        dest="PROTOLEN",
        default="20",
        help="The total protospacer and PAM length.",
    )

    # VEP and PROTEIN LEVEL ANALYSIS
    parser.add_argument(
        "-vep",
        dest="VEP",
        action="store_true",
        help="The boolean option if user wants to analyse the edits through VEP and Uniprot.",
    )

    # MUTATION INFORMATION
    parser.add_argument(
        "-mutation_file",
        dest="MUTATION_FILE",
        default=None,
        type=argparse.FileType("r"),
        help="If you have more than one mutations, a file for the mutations on the "
        "interested gene that you need to integrate into guide and/or annotation analysis",
    )

    # gRNA FLANKING REGIONS
    parser.add_argument(
        "-flank",
        dest="FLAN",
        action="store_true",
        help="The boolean option if the user wants to add flanking sequences of the gRNAs",
    )
    parser.add_argument(
        "-flank3",
        dest="FLAN_3",
        default="7",
        help="The number of nucleotides in the 3' flanking region",
    )
    parser.add_argument(
        "-flank5",
        dest="FLAN_5",
        default="11",
        help="The number of nucleotides in the 5' flanking region",
    )

    # BE INFORMATION

    parser.add_argument(
        "-edit",
        dest="EDIT",
        choices=["A", "T", "G", "C"],
        help="The nucleotide which will be edited.",
    )

    parser.add_argument(
        "-edit_to",
        dest="EDIT_TO",
        choices=["A", "T", "G", "C"],
        help="The nucleotide after edition.",
    )

    # OUTPUT

    parser.add_argument(
        "-o",
        dest=constants.ARGS_KEY_OUTPUT_PATH,
        default=os.getcwd() + "/",
        help="The path for output. If not specified the current directory will be used!",
    )

    parser.add_argument(
        "-ofile",
        dest=constants.ARGS_KEY_OUTPUT_PATH,
        default="output",
        help='The output file name, if not specified "position" will be used!',
    )

    # OFF TARGETS
    parser.add_argument(
        "-ot",
        dest="OFF_TARGET",
        action="store_true",
        help="Whether off targets will be computed or not",
    )
    parser.add_argument(
        "-genome",
        dest="GENOME",
        default="Homo_sapiens.GRCh38.dna.chromosome",
        help="(If -ot provided) name of the genome file",
    )
    parser.add_argument(
        "-v_ensembl",
        dest="VERSION",
        default="113",
        help="The ensembl version in which genome will be retrieved "
        "(if the assembly is GRCh37 then please use <=75)",
    )
    parser.add_argument(
        "-ot_path",
        dest=constants.ARGS_KEY_OT_PATH,
        default=os.getcwd() + "/../offtargets/",
    )

    parsed_input = parser.parse_args()
    input_dict = vars(parsed_input)
    input_dict = _clean_and_globalize_ot_path(input_dict)
    input_dict = _clean_and_globalize_output_path(input_dict)
    return input_dict


def run_offtargets(genome: str, file_name: str, final_df: str) -> bool:
    """
    Perform off-target analysis for identified gRNA sequences.

    Analyzes potential off-target sites for the gRNAs identified in the BEstimate
    analysis using the CRISPRAnalyser tool. Generates detailed and summary reports
    of off-target predictions.

    :param genome: Genome file name/prefix for off-target analysis
    :param file_name: Base name for input and output files
    :param final_df: Suffix for the input dataframe file

    :return: True if off-targets were found and analyzed, False otherwise

    .. note::
        This function requires pre-built genome indices and databases for
        off-target analysis. The analysis uses binary index files and
        SQLite databases created by the x_genome.py script.
    """
    print(f"Summary Data Frame was read from {file_name}{final_df}\n")

    file_prefix = genome.replace(".dna.chromosome", "")
    has_off_targets = x_crispranalyser.get_off_targets(
        input_csv_file=os.path.join(OUTPUT_PATH, f"{file_name}{final_df}"),
        binary_index_file=f"{OT_PATH}grna_bin/{file_prefix}.bin",
        output_csv_file_base=f"{OUTPUT_PATH}{file_name}",
        db_file=f"{OT_PATH}crispr_db/{file_prefix}.db",
    )

    if has_off_targets:
        return True
    else:
        print("No alignment - off target")
        return False


def _clean_and_globalize_output_path(
    orginal_args: dict[str, t.Any],
) -> dict[str, t.Any]:
    arg_key = constants.ARGS_KEY_OUTPUT_PATH
    has_truthy_value = bool(orginal_args.get(arg_key, ""))
    if has_truthy_value:
        raw_value = orginal_args[arg_key]
        has_trailing_slash = raw_value[-1] == "/"
        clean_value = raw_value if has_trailing_slash else raw_value + "/"
    else:
        err_msg = f"An output path must be provided via the '{constants.ARGS_KEY_OUTPUT_PATH}' argument."
        raise RuntimeError(err_msg)
    global OUTPUT_PATH
    OUTPUT_PATH = clean_value
    return orginal_args


def _clean_and_globalize_ot_path(orginal_args: dict[str, t.Any]) -> dict[str, t.Any]:
    # TODO: Target for removal & refactor
    arg_key = constants.ARGS_KEY_OT_PATH
    has_truthy_value = bool(orginal_args.get(arg_key, ""))
    if has_truthy_value:
        raw_value = orginal_args[arg_key]
        has_trailing_slash = raw_value[-1] == "/"
        clean_value = raw_value if has_trailing_slash else raw_value + "/"
    else:
        # TODO: dangerous - relies on cwd being BEstimate/
        clean_value = os.getcwd() + "/../offtargets/"
    global OT_PATH
    OT_PATH = clean_value
    return orginal_args


###########################################################################################
# Execution


def main():
    """
    Execute the complete BEstimate analysis pipeline.

    This is the main execution function that orchestrates the entire base editor
    analysis workflow including gene sequence retrieval, gRNA site identification,
    variant effect prediction, and optional off-target analysis.

    The function processes command line arguments and runs the analysis pipeline
    with the following main steps:

    1. Extract gene information from Ensembl
    2. Identify potential gRNA target sites
    3. Find editable nucleotides within activity windows
    4. Perform VEP annotation (if requested)
    5. Add protein domain and PTM annotations
    6. Generate summary reports
    7. Perform off-target analysis (if requested)

    :raises SystemExit: If no corresponding Ensembl Gene ID is found

    .. note::
        This function uses global variables for arguments and data paths.
        It creates output files in the specified output directory and
        prints progress messages to stdout.
    """
    # Data w/out API opportunity
    yulab_df = DataFiles.get_homo_sapiens_interfaces_as_dataframe()
    global CLI_ARGS
    CLI_ARGS = take_input()

    print(
        """
--------------------------------------------------------------
           B E s t i m a t e

           Wellcome Sanger Institute

--------------------------------------------------------------
    """
    )
    if CLI_ARGS["VEP"]:
        vep = True
    else:
        vep = False

    if CLI_ARGS["OFF_TARGET"]:
        ot_analysis = True
    else:
        ot_analysis = False

    if CLI_ARGS["TRANSCRIPT"]:
        transcript = CLI_ARGS["TRANSCRIPT"]
    else:
        transcript = None

    if CLI_ARGS["MUTATION_FILE"]:
        mutations = list()
        for line in CLI_ARGS["MUTATION_FILE"].readlines():
            mutations.append(line.strip())
    else:
        mutations = None

    print(
        """
The given arguments are:\nGene: %s\nAssembl: %s\nEnsembl transcript ID: %s\nUniprot ID: %s\nPAM sequence: %s\nPAM window: %s
Protospacer length: %s\nActivity window: %s\nNucleotide change: %s>%s\nVEP and Uniprot analysis: %s\nMutation on genome: %s
Off target analysis: %s"""
        % (
            CLI_ARGS["GENE"],
            CLI_ARGS["ASSEMBLY"],
            CLI_ARGS["TRANSCRIPT"],
            CLI_ARGS["UNIPROT"],
            CLI_ARGS["PAMSEQ"],
            CLI_ARGS["PAMWINDOW"],
            CLI_ARGS["PROTOLEN"],
            CLI_ARGS["ACTWINDOW"],
            CLI_ARGS["EDIT"],
            CLI_ARGS["EDIT_TO"],
            vep,
            ", ".join("" if mutations is None else mutations),
            ot_analysis,
        )
    )

    print(
        """\n

--------------------------------------------------------------
        Ensembl Gene Information
--------------------------------------------------------------
    \n"""
    )

    ensembl_obj = Ensembl(hugo_symbol=CLI_ARGS["GENE"], assembly=CLI_ARGS["ASSEMBLY"])

    ensembl_obj.extract_gene_id()

    if ensembl_obj.gene_id == "":
        sys.exit("No corresponding Ensembl Gene ID could be found!")

    ensembl_obj.extract_sequence(ensembl_obj.gene_id, mutations=mutations)

    if ensembl_obj.gene_range[0] < ensembl_obj.gene_range[1]:
        ensembl_obj.extract_info(
            chromosome=ensembl_obj.chromosome,
            loc_start=ensembl_obj.gene_range[0],
            loc_end=ensembl_obj.gene_range[1],
            transcript=transcript,
        )
    else:
        ensembl_obj.extract_info(
            chromosome=ensembl_obj.chromosome,
            loc_start=ensembl_obj.gene_range[1],
            loc_end=ensembl_obj.gene_range[0],
            transcript=transcript,
        )

    print(
        """\n
--------------------------------------------------------------
        gRNAs - Targetable Sites
--------------------------------------------------------------
    \n"""
    )
    try:
        os.mkdir(OUTPUT_PATH)
    except FileExistsError:
        pass

    file_name = CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_edit_df.csv"

    file_name = CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_edit_df.csv"

    if CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_crispr_df.csv" not in os.listdir(
        OUTPUT_PATH
    ):
        crispr_df = seq_utils.extract_grna_sites(
            hugo_symbol=CLI_ARGS["GENE"],
            searched_nucleotide=CLI_ARGS["EDIT"],
            pam_window=[
                int(CLI_ARGS["PAMWINDOW"].split("-")[0]),
                int(CLI_ARGS["PAMWINDOW"].split("-")[1]),
            ],
            activity_window=[
                int(CLI_ARGS["ACTWINDOW"].split("-")[0]),
                int(CLI_ARGS["ACTWINDOW"].split("-")[1]),
            ],
            pam_sequence=CLI_ARGS["PAMSEQ"],
            protospacer_length=CLI_ARGS["PROTOLEN"],
            flan=CLI_ARGS["FLAN"],
            flan_3=CLI_ARGS["FLAN_3"],
            flan_5=CLI_ARGS["FLAN_5"],
            ensembl_object=ensembl_obj,
        )

        if len(crispr_df.index) != 0:
            print("CRISPR Data Frame was created!")
        crispr_df.to_csv(
            OUTPUT_PATH + CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_crispr_df.csv",
            index=False,
        )

        print(
            "CRISPR Data Frame was written in %s as %s\n"
            % (OUTPUT_PATH, CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_crispr_df.csv")
        )

    else:
        print(
            "CRISPR Data Frame was read from %s as %s\n\n"
            % (OUTPUT_PATH, CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_crispr_df.csv")
        )
        crispr_df = pandas.read_csv(
            OUTPUT_PATH + CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_crispr_df.csv"
        )
    print(
        """\n
--------------------------------------------------------------
        gRNAs - Editable Sites
--------------------------------------------------------------
    \n"""
    )
    if CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_edit_df.csv" not in os.listdir(
        OUTPUT_PATH
    ):
        edit_df = seq_utils.find_editable_nucleotide(
            crispr_df=crispr_df,
            searched_nucleotide=CLI_ARGS["EDIT"],
            activity_window=[
                int(CLI_ARGS["ACTWINDOW"].split("-")[0]),
                int(CLI_ARGS["ACTWINDOW"].split("-")[1]),
            ],
            pam_window=[
                int(CLI_ARGS["PAMWINDOW"].split("-")[0]),
                int(CLI_ARGS["PAMWINDOW"].split("-")[1]),
            ],
            ensembl_object=ensembl_obj,
            mutations=mutations,
        )

        if len(edit_df.index) != 0:
            print("Edit Data Frame was created!")

        edit_df.to_csv(
            OUTPUT_PATH + CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_edit_df.csv",
            index=False,
        )

        print(
            "Edit Data Frame was written in %s as %s"
            % (OUTPUT_PATH, CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_edit_df.csv\n")
        )

    else:
        print(
            "Edit Data Frame was read from %s as %s\n\n"
            % (OUTPUT_PATH, CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_edit_df.csv")
        )
        edit_df = pandas.read_csv(
            OUTPUT_PATH + CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_edit_df.csv"
        )

    if CLI_ARGS["VEP"]:
        print(
            """\n
--------------------------------------------------------------
        Annotation - VEP Annotation
--------------------------------------------------------------
        \n"""
        )
        file_name = CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_summary_df.csv"
        whole_vep_df = pandas.DataFrame()
        if CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_vep_df.csv" not in os.listdir(
            OUTPUT_PATH
        ):
            if CLI_ARGS[
                constants.ARGS_KEY_OUTPUT_PATH
            ] + "_hgvs_df.csv" not in os.listdir(OUTPUT_PATH):
                hgvs_df = seq_utils.extract_hgvs_df(
                    edit_df=edit_df,
                    ensembl_object=ensembl_obj,
                    transcript_id=CLI_ARGS["TRANSCRIPT"],
                    edited_nucleotide=CLI_ARGS["EDIT"],
                    new_nucleotide=CLI_ARGS["EDIT_TO"],
                    activity_window=[
                        int(CLI_ARGS["ACTWINDOW"].split("-")[0]),
                        int(CLI_ARGS["ACTWINDOW"].split("-")[1]),
                    ],
                    mutations=mutations,
                )
                if hgvs_df is not None and len(hgvs_df.index) != 0:
                    hgvs_df.to_csv(
                        OUTPUT_PATH
                        + CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH]
                        + "_hgvs_df.csv"
                    )
                    print("HGVS nomenclatures were collected.\n")
            else:
                hgvs_df = pandas.read_csv(
                    OUTPUT_PATH
                    + CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH]
                    + "_hgvs_df.csv"
                )
                print("HGVS nomenclatures were collected.\n")

            if hgvs_df is not None and len(hgvs_df.index) != 0:
                whole_vep_df = seq_utils.retrieve_vep_info(
                    hgvs_df=hgvs_df,
                    ensembl_object=ensembl_obj,
                    uniprot=CLI_ARGS["UNIPROT"],
                    transcript_id=CLI_ARGS["TRANSCRIPT"],
                )
                if len(whole_vep_df.index) != 0:
                    print("VEP Data Frame was created!")
                    whole_vep_df.to_csv(
                        OUTPUT_PATH
                        + CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH]
                        + "_vep_df.csv"
                    )
                    print(
                        "VEP Data Frame was written in %s as %s\n\n"
                        % (
                            OUTPUT_PATH,
                            CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_vep_df.csv",
                        )
                    )
                else:
                    print("VEP Data Frame cannot be created because it is empty!")
        else:
            print(
                "VEP Data Frame was read from %s as %s\n\n"
                % (
                    OUTPUT_PATH,
                    CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_vep_df.csv",
                )
            )
            whole_vep_df = pandas.read_csv(
                OUTPUT_PATH + CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_vep_df.csv"
            )

        print(
            """\n
--------------------------------------------------------------
        Annotation - Protein Annotation
--------------------------------------------------------------
        \n"""
        )
        protein_df = pandas.DataFrame()
        if CLI_ARGS[
            constants.ARGS_KEY_OUTPUT_PATH
        ] + "_protein_df.csv" not in os.listdir(OUTPUT_PATH):
            print("Adding Uniprot ID, corresponding Domain and PTM information..")
            if len(whole_vep_df.index) != 0:
                uniprot_df = seq_utils.annotate_edits(
                    ensembl_object=ensembl_obj,
                    vep_df=whole_vep_df,
                    uniprot_id=CLI_ARGS["UNIPROT"],
                    output_path=OUTPUT_PATH,
                )
                if uniprot_df is not None and len(uniprot_df.index) != 0:
                    print("Adding affected interface and interacting partners..")
                    protein_df = seq_utils.annotate_interface(
                        annotated_edit_df=uniprot_df,
                        uniprot_id=CLI_ARGS["UNIPROT"],
                        yulab_df=yulab_df,
                    )

                    if protein_df is not None and len(protein_df.index) != 0:
                        print("Protein Data Frame was created!")
                        protein_output_file = os.path.join(
                            OUTPUT_PATH,
                            CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH]
                            + "_protein_df.csv",
                        )
                        protein_df.to_csv(protein_output_file, index=False)
                        print(
                            f"Protein Data Frame was written as {protein_output_file}\n"
                        )
                    else:
                        print(
                            "Protein Data Frame cannot be created because it is empty."
                        )
                else:
                    print("Protein Data Frame cannot be created because it is empty.")
            else:
                print("Protein Data Frame cannot be created because it is empty.")
        else:
            print(
                "Protein Data Frame was read from %s as %s\n\n"
                % (
                    OUTPUT_PATH,
                    CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_protein_df.csv",
                )
            )
            protein_df = pandas.read_csv(
                OUTPUT_PATH
                + CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH]
                + "_protein_df.csv"
            )

        if len(protein_df.index) > 0:
            if CLI_ARGS[
                constants.ARGS_KEY_OUTPUT_PATH
            ] + "_summary_df.csv" not in os.listdir(OUTPUT_PATH):
                print("Summarising information..")
                summary_df = seq_utils.summarise_guides(last_df=protein_df)

                if summary_df is not None and len(summary_df.index) != 0:
                    print("Summary Data Frame was created!")
                    summary_df.to_csv(
                        OUTPUT_PATH
                        + CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH]
                        + "_summary_df.csv",
                        index=False,
                    )
                    final_df = "_summary_df.csv"
                    print(
                        "Summary Data Frame was written in %s as %s\n\n"
                        % (
                            OUTPUT_PATH,
                            CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH]
                            + "_summary_df.csv",
                        )
                    )
                else:
                    print("Summary Data Frame cannot be created because it is empty.")

            else:
                print(
                    "Summary Data Frame was read from %s as %s\n\n"
                    % (
                        OUTPUT_PATH,
                        CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH] + "_summary_df.csv",
                    )
                )
                summary_df = pandas.read_csv(
                    OUTPUT_PATH
                    + CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH]
                    + "_summary_df.csv"
                )
                final_df = "_summary_df.csv"
        else:
            print("Protein Data Frame cannot be created because it is empty.")

    else:
        final_df = "_edit_df.csv"

    if CLI_ARGS["OFF_TARGET"]:
        print(
            """\n
--------------------------------------------------------------
        Annotation - Off Target Annotation
--------------------------------------------------------------
                \n"""
        )
        try:
            os.mkdir(os.getcwd() + "/../offtargets")
        except FileExistsError:
            pass

        try:
            os.mkdir(os.getcwd() + "/../offtargets/output/")

        except FileExistsError:
            pass
        if CLI_ARGS["ASSEMBLY"] == "GRCh37":
            file_main_text = "Homo_sapiens.GRCh37.%s.%s" % (
                CLI_ARGS["VERSION"],
                CLI_ARGS["PAMSEQ"],
            )
        elif CLI_ARGS["ASSEMBLY"] == "GRCh38":
            file_main_text = "Homo_sapiens.GRCh38.%s" % CLI_ARGS["PAMSEQ"]
        if "%s.bin" % file_main_text in os.listdir(f"{OT_PATH}grna_bin/"):
            _ = run_offtargets(
                genome=file_main_text,
                file_name=CLI_ARGS[constants.ARGS_KEY_OUTPUT_PATH],
                final_df=final_df,
            )
        else:
            print("Please download and index your genome file\nRun x_genome.py first.")
    print(
        """\n
--------------------------------------------------------------
    The BEstimate analysis successfully finished!
--------------------------------------------------------------
    \n"""
    )
    return


if __name__ == "__main__":
    deprecation_msg = (
        f"You cannot run 'python {__file__}' directly, please use the CLI command "
        f"'{constants.PROGRAM_NAME}' instead."
    )
    warnings.warn(deprecation_msg, DeprecationWarning)
    sys.exit(1)
