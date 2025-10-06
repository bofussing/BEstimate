import subprocess
import filecmp
import shutil


# tests that off-target files are correctly generated
def test_output(tmp_path):
    d = tmp_path / "test"
    d.mkdir()
    shutil.copytree("data/offtargets", d / "offtargets")
    subprocess.run(
        [
            "BEstimate",
            "-gene",
            "SRY",
            "-assembly",
            "GRCh38",
            "-pamseq",
            "NGG",
            "-edit",
            "A",
            "-edit_to",
            "G",
            "-ot",
            "-o",
            "/tmp/output",
            "-ot_path",
            "offtargets",
            "-ofile",
            "SRY_TEST",
        ],
        cwd=d,
    )

    for file in ["crispr_df", "ot_annotated_df", "edit_df", "wge_return"]:
        assert filecmp.cmp(
            d / f"SRY_TEST/SRY_TEST_{file}.csv",
            d / f"offtargets/output/SRY_TEST_{file}.csv",
        )
