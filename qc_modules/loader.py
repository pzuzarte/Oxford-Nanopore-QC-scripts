import glob
import random
import pandas as pd

# Columns required for all analyses
CORE_COLUMNS = [
    "read_id",
    "channel",
    "mux",
    "start_time",
    "duration",
    "passes_filtering",
    "template_start",
    "num_events_template",
    "template_duration",
    "sequence_length_template",
    "mean_qscore_template",
]

# Optional columns present in some file versions
OPTIONAL_COLUMNS = [
    "median_template",
    "mad_template",
    "strand_score_template",
    "minknow_events",
    "pore_type",
    "experiment_id",
    "sample_id",
    "end_reason",
]

BARCODE_COLUMNS = [
    "barcode_arrangement",
    "barcode_score",
    "barcode_kit",
]

# ---------------------------------------------------------------------------
# Old-format column aliases
# ---------------------------------------------------------------------------
# Maps original file column names (Albacore / early Guppy) to the canonical
# names used throughout this suite.  Applied during loading so all downstream
# code only needs to know the canonical names.
#
# Albacore (v2.x and older) differences:
#   - No  passes_filtering column  → derived from Q-score threshold
#   - strand_score                 → strand_score_template
#   - num_called_template          → num_events_template
#   - mean_qscore                  → mean_qscore_template
#   - sequence_length              → sequence_length_template
#   - duration (also the read dur) → kept as-is; template_duration may be absent
#   - called_events                → num_events_template  (very old Albacore)
#   - start_time already canonical in most versions
#
# MinKNOW / early Guppy variants:
#   - channel_number               → channel
#   - read_number                  → (not mapped; not currently used)
#   - start_mux                    → mux
_COLUMN_ALIASES = {
    # Q-score ----------------------------------------------------------------
    'mean_qscore':              'mean_qscore_template',   # Albacore ≤ 2.x
    'qscore':                   'mean_qscore_template',   # occasional variant
    'q_score':                  'mean_qscore_template',   # occasional variant
    'quality_mean':             'mean_qscore_template',   # rare MinKNOW variant

    # Read length ------------------------------------------------------------
    'sequence_length':          'sequence_length_template',  # Albacore ≤ 2.x
    'read_length':              'sequence_length_template',  # occasional variant

    # Events -----------------------------------------------------------------
    'num_called_template':      'num_events_template',    # Albacore variant
    'called_events':            'num_events_template',    # very old Albacore
    'num_events':               'num_events_template',    # rare variant

    # Strand score -----------------------------------------------------------
    'strand_score':             'strand_score_template',  # Albacore (no suffix)

    # Channel / mux ----------------------------------------------------------
    'channel_number':           'channel',               # early MinKNOW
    'start_mux':                'mux',                   # early MinKNOW

    # Template duration ------------------------------------------------------
    'duration':                 'template_duration',     # Albacore: one duration col
}

# Columns that MUST be present after alias resolution for the suite to work.
# Used only for diagnostics — a clear error is raised if any are missing.
_REQUIRED_COLUMNS = [
    'sequence_length_template',
    'mean_qscore_template',
    'start_time',
    'channel',
]

# Q-score threshold used to derive passes_filtering when the column is absent
# (Albacore-era files predate the column; ONT's historical default was Q7).
_DEFAULT_QSCORE_THRESHOLD = 7.0


def find_summary_file(search_dir='.'):
    matches = glob.glob(f"{search_dir}/sequencing_summary*.txt")
    if not matches:
        raise FileNotFoundError(
            f"No sequencing_summary*.txt file found in '{search_dir}'. "
            "Use --file to specify the path explicitly."
        )
    if len(matches) > 1:
        print(f"[WARNING] Multiple summary files found: {matches}. Using {matches[0]}.")
    return matches[0]


def find_pore_activity_file(search_dir='.'):
    matches = glob.glob(f"{search_dir}/pore_activity*.csv")
    if not matches:
        return None
    return matches[0]


def find_throughput_file(search_dir='.'):
    matches = glob.glob(f"{search_dir}/throughput_*.csv")
    if not matches:
        return None
    return matches[0]


def load_summary(filepath, subsample=1.0):
    """
    Load a sequencing summary TSV file.

    Handles both current (Guppy / Dorado) and legacy (Albacore / early Guppy)
    column naming conventions.  Key differences handled automatically:

      passes_filtering        Missing in Albacore files → derived from
                              mean_qscore_template >= 7.0 with a printed warning.
      mean_qscore             Albacore name for mean_qscore_template → renamed.
      sequence_length         Albacore name for sequence_length_template → renamed.
      strand_score            Albacore name for strand_score_template → renamed.
      num_called_template /
        called_events         Albacore names for num_events_template → renamed.
      duration (single col)   Albacore uses one duration column; mapped to
                              template_duration when no template_duration exists.

    Parameters
    ----------
    filepath : str
        Path to the sequencing_summary*.txt file.
    subsample : float
        Fraction of rows to retain (0 < subsample <= 1).

    Returns
    -------
    df : pd.DataFrame
    is_barcoded : bool
        True if the file contains barcode_arrangement column.
    """
    # --- Peek at header -------------------------------------------------------
    header = pd.read_csv(filepath, sep='\t', nrows=0)
    file_cols = set(header.columns)

    # Build a map of original-file-column → canonical-name for aliased columns.
    # Special case: only map 'duration' → 'template_duration' when the file
    # has NO 'template_duration' column of its own (avoids clobbering Guppy
    # files that carry both duration and template_duration).
    col_rename = {}
    for orig, canon in _COLUMN_ALIASES.items():
        if orig not in file_cols:
            continue
        if orig == 'duration' and 'template_duration' in file_cols:
            continue   # modern file — keep both columns as-is
        col_rename[orig] = canon

    if col_rename:
        print(f'[loader] Old-format columns detected and renamed: {col_rename}')

    # All canonical column names we want
    wanted_canonical = set(CORE_COLUMNS + OPTIONAL_COLUMNS + BARCODE_COLUMNS)

    # Select original file column names that map to something we want
    usecols = [
        col for col in header.columns
        if col_rename.get(col, col) in wanted_canonical
    ]

    # --- Load -----------------------------------------------------------------
    if subsample < 1.0:
        df = pd.read_csv(
            filepath,
            sep='\t',
            usecols=usecols,
            skiprows=lambda i: i > 0 and random.random() > subsample,
        )
    else:
        df = pd.read_csv(filepath, sep='\t', usecols=usecols)

    # Apply column renames so all downstream code uses canonical names
    if col_rename:
        df.rename(columns=col_rename, inplace=True)

    # --- Diagnostics: check required columns are present ---------------------
    missing = [c for c in _REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        # Detect the common case: a MinKNOW acquisition summary (raw run index)
        # which has read_id / channel / start_time but NO basecall quality columns.
        # These are written by MinKNOW itself; the quality columns come from the
        # basecaller (Dorado / Guppy) which writes its own sequencing_summary file.
        _acquisition_indicators = {'filename_pod5', 'filename_fast5',
                                   'filename_fastq', 'parent_read_id', 'run_id'}
        _quality_indicators = {
            'sequence_length_template', 'sequence_length',
            'mean_qscore_template', 'mean_qscore', 'qscore', 'q_score',
        }
        looks_like_acquisition = (
            bool(_acquisition_indicators & file_cols)
            and not bool(_quality_indicators & file_cols)
        )

        if looks_like_acquisition:
            print(
                f'\n[loader] ERROR: This file appears to be a MinKNOW ACQUISITION '
                f'summary (a raw run index), not a basecaller summary.\n'
                f'\n'
                f'  It has file-path / run-ID columns but NO quality or length columns.\n'
                f'  Basecall quality data (Q-score, read length, etc.) are written by\n'
                f'  the basecaller (Dorado or Guppy) into a separate file, typically\n'
                f'  also called sequencing_summary.txt but located in the basecall\n'
                f'  output directory.\n'
                f'\n'
                f'  Columns found in this file: {sorted(file_cols)}\n'
                f'\n'
                f'  What to do:\n'
                f'    1. Locate the Dorado / Guppy basecall output directory.\n'
                f'    2. Pass that summary file explicitly:\n'
                f'         python ont_qc.py --file /path/to/basecall/sequencing_summary.txt\n'
            )
        else:
            print(
                f'\n[loader] ERROR: After alias resolution the following required '
                f'columns are still missing: {missing}\n'
                f'         Columns present in the file: {sorted(file_cols)}\n'
                f'         Please check that the file is a valid sequencing_summary TSV.\n'
            )
        raise KeyError(
            f"Required columns missing after loading: {missing}. "
            f"File columns were: {sorted(file_cols)}"
        )

    # --- passes_filtering -----------------------------------------------------
    if 'passes_filtering' not in df.columns:
        # Albacore-era files do not have this column.  Derive it from Q-score.
        print(
            f'[loader] passes_filtering column not found (likely an Albacore-era file). '
            f'Deriving from mean_qscore_template >= {_DEFAULT_QSCORE_THRESHOLD}.'
        )
        df['passes_filtering'] = (
            df['mean_qscore_template'] >= _DEFAULT_QSCORE_THRESHOLD
        )
    else:
        # Normalise to Python bool.
        # The column arrives as strings ('True'/'False'), Python bools stored as
        # object dtype (common in mid-run files), or already bool/int.
        pf = df['passes_filtering']
        if pf.dtype == object:
            df['passes_filtering'] = pf.map(
                lambda v: str(v).strip().upper() == 'TRUE'
            )
        elif pf.dtype != bool:
            df['passes_filtering'] = pf.astype(bool)

    # --- sequencing_speed -----------------------------------------------------
    if 'template_duration' in df.columns:
        df['sequencing_speed'] = (
            df['sequence_length_template']
            / df['template_duration'].replace(0, float('nan'))
        )
    else:
        df['sequencing_speed'] = float('nan')

    is_barcoded = 'barcode_arrangement' in df.columns

    return df, is_barcoded
