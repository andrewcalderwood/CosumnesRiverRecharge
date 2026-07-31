"""
copy_model_outputs.py

Copies output_clean folders from groundwater model run directories
to a central output location, preserving the model run folder name.

Source:  F:/WRDAPP/GWFlowModel/Cosumnes/Economic/<model_run>/output_clean/
Dest:    C:/Users/andrew/Box/SESYNC_paper1/model_results/2026_3_20/<model_run>/output_clean/
"""

import shutil
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
loadpth = pd.read_csv('00_model_path.txt',header=None).iloc[0,0]
# SOURCE_BASE = Path(r"D:/WRDAPP/GWFlowModel/Cosumnes/Economic")
SOURCE_BASE = Path(loadpth)

import os
usr_dir = os.path.expanduser('~')
DEST_BASE   = Path(usr_dir+r"/Box/SESYNC_paper1/model_results/2026_7_27")
print(DEST_BASE)

# +
# MODEL_RUNS = [
#     "input_write_2014_2022_R200",
#     "input_write_2014_2022_R203",
#     "input_write_2014_2022_R204",
#     "input_write_2016_2022_R200",
#     "input_write_2018_2022_R200",
# ]
# -

MODEL_RUNS = [
    # "input_write_2014_2025_R200",
    # "input_write_2014_2025_R203",
    # "input_write_2014_2025_R204",
    "input_write_2014_2025_R300",
    "input_write_2014_2025_R303",
    "input_write_2014_2025_R304",
]

OUTPUT_SUBFOLDER = "output_clean"

# ---------------------------------------------------------------------------
# Copy logic
# ---------------------------------------------------------------------------
def copy_output_clean(dry_run: bool = False) -> None:
    """
    Iterate over each model run and copy its output_clean folder.

    Parameters
    ----------
    dry_run : bool
        If True, print what would be copied without actually copying anything.
    """
    print(f"{'[DRY RUN] ' if dry_run else ''}Starting copy operation.../n")

    for run in MODEL_RUNS:
        src = SOURCE_BASE / run / OUTPUT_SUBFOLDER
        dst = DEST_BASE   / run / OUTPUT_SUBFOLDER

        if not src.exists():
            print(f"  [SKIP]  Source not found: {src}")
            continue

        if dst.exists():
            print(f"  [WARN]  Destination already exists (will overwrite): {dst}")

        if dry_run:
            print(f"  [DRY]   Would copy:/n          {src}/n       -> {dst}")
        else:
            # Remove existing destination so copytree starts fresh
            if dst.exists():
                shutil.rmtree(dst)

            shutil.copytree(src, dst)
            print(f"  [OK]    Copied:/n          {src}/n       -> {dst}")

    print("/nDone.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Copy output_clean folders from model runs to a central location."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview what would be copied without making any changes.",
    )
    args = parser.parse_args()

    copy_output_clean(dry_run=args.dry_run)




