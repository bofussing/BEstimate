import subprocess
import filecmp
import os


# tests that off-target files are correctly generated
def test_output(tmp_path, request):
    d = tmp_path / "test"
    d.mkdir()
    test_dir = os.path.dirname(request.path)
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
            "-ot_path",
            f"{test_dir}/data/offtargets/",
            "-ofile",
            "SRY_TEST",
        ],
        cwd=d,
    )

    for file in ["crispr_df", "ot_annotated_df", "edit_df", "wge_return"]:
        assert filecmp.cmp(
            d / f"SRY_TEST/SRY_TEST_{file}.csv",
            d / f"{test_dir}/data/offtargets/output/SRY_TEST_{file}.csv",
        )
