
import pandas as pd
import statistics
from typing import Optional

CONTROL_COLS = [
    "avg_cov_ACAD9",
    "avg_cov_ATR",
    "avg_cov_CYP11B1",
    "avg_cov_EDNRB",
    "avg_cov_FASTKD2",
    "avg_cov_FOXN1",
    "avg_cov_HEXB",
    "avg_cov_IQCB1",
    "avg_cov_ITGA6",
    "avg_cov_IVD",
    "avg_cov_LMNA",
    "avg_cov_LRPPRC",
    "avg_cov_NTRK1",
    "avg_cov_PTEN",
    "avg_cov_RAB3GAP1",
    "avg_cov_RAPSN",
    "avg_cov_SIL1",
    "avg_cov_SLC22A5",
    "avg_cov_SLC35D1",
    "avg_cov_STIM1",
]

def parse_duplication_markers(row):
    """
    Parse SMAca markers for silent carriers outputs:

    - g.27134T>G:
        Extract the token before the space, take its *last nucleotide*.
        If last base == 'G' → ALT (duplication marker present)
        If last base == 'T' → REF
        Else → treat as REF (conservative)

    - g.27706_27707delAT:
        Extract the token before the space.
        If it contains substring 'AT' → REF
        Else → ALT (marker present)
    """

    def _get_token(value):
        """Return the token before the first whitespace."""
        if pd.isna(value):
            return ""
        s = str(value).strip()
        if not s:
            return ""
        return s.split()[0]

    # Extract tokens
    tok_27134 = _get_token(row.get("g.27134T>G", None))
    tok_27706 = _get_token(row.get("g.27706_27707delAT", None))

    # ----- Marker 1: g.27134T>G -----
    # Last nucleotide determines the allele
    last_base = tok_27134[-1] if tok_27134 else ""

    if last_base == "G":
        has_27134_alt = True
    else:
        has_27134_alt = False  # 'T' or anything else -> treat as reference

    # ----- Marker 2: g.27706_27707delAT -----
    # If token contains 'AT' → reference (no deletion)
    if tok_27706 and ("AT" in tok_27706):
        has_27706_delAT = False
    else:
        has_27706_delAT = True

    return has_27134_alt, has_27706_delAT

def classify_smaca_sample(df_smaca, smaca_cv_fail_threshold, smaca_cv_warn_threshold, smaca_low_cov_abs, smaca_low_cov_rel):
    '''

    Apply the SMA carrier algorithm to a single sample's metrics.

       QC logic (updated):
       -------------------
       - CV_control = std_control / mean_control
       - If CV_control <= cv_warn_threshold    → qc_flag = "QC_OK"
       - If cv_warn_threshold < CV_control <= cv_fail_threshold → "QC_WARN"
       - If CV_control > cv_fail_threshold     → "QC_FAIL"
       - If CV_control cannot be computed      → "QC_UNKNOWN"

       - SMN_coverage_flag (soft warning, does NOT stop classification):
           * compute mean_cov_SMN1, mean_cov_SMN2
           * "LOW_COVERAGE_WARNING" if both:
                - mean_cov_SMN1 < low_cov_abs AND mean_cov_SMN2 < low_cov_abs
              OR
                - mean_cov_SMN1 and mean_cov_SMN2 < low_cov_rel * mean_control
           * otherwise "OK"

       Carrier logic (unchanged in structure, thresholds refined):
       -----------------------------------------------------------
       - Classical carriers:
           any Pi_j < 1/3  → "LIKELY_SMA_CARRIER (1-copy SMN1)"

       - Inconclusive:
           no Pi_j < 1/3  AND  Pi_mean < 0.40 → "Inconclusive"

       - Otherwise:
           provisional "LIKELY_NON_CARRIER"

       - Silent carriers:
           if provisional "LIKELY_NON_CARRIER" AND both duplication markers
           (27134 non-ref AND 27706 non-ref) → "PUTATIVE_SILENT_SMA_CARRIER"

    :param df_smaca: Data frame with SMAca output
    :param smaca_cv_fail_threshold: Max CV_control for QC warning
    :param smaca_cv_warn_threshold: Max CV_control for QC OK
    :param smaca_low_cov_abs: Absolute mean coverage (×) below which both SMN1 and SMN2 trigger a low-coverage warning.
    :param smaca_low_cov_rel: raction of mean control-gene coverage; if both SMN1 and SMN2 fall below this fraction, low-coverage warning.
    :return:
    '''


    row = df_smaca.iloc[0]

    avg_cov_control_genes = [float(row[c]) for c in CONTROL_COLS]
    has_27134_nonref, has_27706_nonref = parse_duplication_markers(row)

    Pi_a = float(row["Pi_a"])
    Pi_b = float(row["Pi_b"])
    Pi_c = float(row["Pi_c"])
    scale_factor = float(row["scale_factor"])
    cov_SMN1_a = float(row["cov_SMN1_a"])
    cov_SMN1_b = float(row["cov_SMN1_b"])
    cov_SMN1_c = float(row["cov_SMN1_c"])
    cov_SMN2_a = float(row["cov_SMN2_a"])
    cov_SMN2_b = float(row["cov_SMN2_b"])
    cov_SMN2_c = float(row["cov_SMN2_c"])
    std_control = float(row["std_control"])
    avg_cov_control_genes = avg_cov_control_genes
    has_27134_nonref = has_27134_nonref
    has_27706_nonref = has_27706_nonref
    cv_fail_threshold = smaca_cv_fail_threshold
    cv_warn_threshold = smaca_cv_warn_threshold
    low_cov_abs = smaca_low_cov_abs
    low_cov_rel = smaca_low_cov_rel

    # --- Derived Pi statistics ----------------------------------------
    Pi_values = [Pi_a, Pi_b, Pi_c]
    Pi_mean = statistics.mean(Pi_values)
    Pi_min = min(Pi_values)
    Pi_max = max(Pi_values)

    # --- Coverage at the SMN locus -----------------------------------
    # Per-gene means (for SMN coverage flag)
    mean_cov_SMN1 = statistics.mean([cov_SMN1_a, cov_SMN1_b, cov_SMN1_c])
    mean_cov_SMN2 = statistics.mean([cov_SMN2_a, cov_SMN2_b, cov_SMN2_c])
    # Overall mean (for reporting)
    mean_cov_SMN = statistics.mean(
        [cov_SMN1_a, cov_SMN1_b, cov_SMN1_c,
         cov_SMN2_a, cov_SMN2_b, cov_SMN2_c]
    )

    # --- Control-gene coverage & CV_control ---------------------------
    if avg_cov_control_genes:
        mean_control = statistics.mean(avg_cov_control_genes)
    else:
        mean_control = 0.0

    if mean_control > 0:
        CV_control: Optional[float] = std_control / mean_control
    else:
        CV_control = None

    # --- Total copies and dup markers --------------------------------
    total_copies = round(4 * scale_factor)
    has_dup_marker = has_27134_nonref or has_27706_nonref


    # ==============================================================    # STEP 1 – Quality control
    # ==============================================================

    # QC flag from CV_control
    if CV_control is None:
        qc_flag = "QC_UNKNOWN"
    else:
        if CV_control <= cv_warn_threshold:
            qc_flag = "QC_OK"
        elif CV_control <= cv_fail_threshold:
            qc_flag = "QC_WARN"
        else:
            qc_flag = "QC_FAIL"

    # SMN coverage soft flag
    def _is_low_abs(val: float) -> bool:
        return val < low_cov_abs

    def _is_low_rel(val: float, mean_ctrl: float) -> bool:
        return (mean_ctrl > 0) and (val < low_cov_rel * mean_ctrl)

    both_abs_low = _is_low_abs(mean_cov_SMN1) and _is_low_abs(mean_cov_SMN2)
    both_rel_low = _is_low_rel(mean_cov_SMN1, mean_control) and _is_low_rel(
        mean_cov_SMN2, mean_control
    )

    if both_abs_low or both_rel_low:
        SMN_coverage_flag = "LOW_COVERAGE_WARNING"
    else:
        SMN_coverage_flag = "OK"

    # Classification goes ahead; qc_flag and SMN_coverage_flag are just annotations.

    # ==============================================================    # STEP 2 – Classic carriers (1-copy SMN1)
    #   Rule from SMAca: Pi_a < 1/3, Pi_b < 1/3 or Pi_c < 1/3.
    # ==============================================================
    # SMN2 gene is not included in RR catalogue, so no matter the estimated number of SMN2 copies (scale factor is not used)

    call = "None"

    if any(p < (1.0 / 3.0) for p in Pi_values):
        call = "LIKELY_SMA_CARRIER (1-copy SMN1)"

    # ==============================================================    # STEP 3 – Inconclusive/Check manually (any p > 0.33 but Pi_mean < 0.40
    # ==============================================================
    elif Pi_mean < 0.40:
        call = "Inconclusive"
    # ==============================================================    # STEP 4 – Likely non-carriers (any p > 0.33 and Pi_mean >= 0.40)
    # ==============================================================
    else:
        call = "LIKELY_NON_CARRIER"
    # ==============================================================    # STEP 5 – Silent carriers
    # ==============================================================

    if call == "LIKELY_NON_CARRIER" and has_dup_marker:
        call = "PUTATIVE_SILENT_SMA_CARRIER"

    # Final result dict
    result_dict= {
        "call": call,
        "qc_flag": qc_flag,
        "SMN_coverage_flag": SMN_coverage_flag,
        "Pi_mean": Pi_mean,
        "Pi_range": (Pi_min, Pi_max),
        "scale_factor": scale_factor,
        "total_copies_est": total_copies,
        "mean_cov_SMN": mean_cov_SMN,
        "mean_cov_SMN1": mean_cov_SMN1,
        "mean_cov_SMN2": mean_cov_SMN2,
        "mean_control": mean_control,
        "CV_control": CV_control,
        "related_HPOs_for_sample": 'NA'
    }

    return result_dict


def read_smaca_file(SMAca_output_file):
    """
    Read SMAca output where:
      - The header is comma-separated and starts with '#'
      - Data rows are '|' separated
      - Lines starting with '#' (after the header) are comments and must be ignored
      - All columns are converted to float except:
          * id        (string)
          * g.27134T>G
          * g.27706_27707delAT
    """
    with open(SMAca_output_file, "r") as fh:
        # --- Read and parse header (comma-separated, starts with '#') ---
        header_line = fh.readline().rstrip("\n")
        # Remove leading '#' and spaces, then split by comma
        header_line = header_line.lstrip("#").strip()
        columns = [c.strip() for c in header_line.split(",")]

        for line in fh:
            line = line.strip()
            if not line.startswith("#"):
                row = line.split("|")

    # Build initial DataFrame as strings
    df = pd.DataFrame([row], columns=columns)

    print(df)
    # Replace empty strings with NaN
    df = df.replace({"": pd.NA})

    # Columns that must stay as strings
    marker_cols = ["g.27134T>G", "g.27706_27707delAT"]
    string_cols = ["id"] + marker_cols

    # Convert all other columns to float
    for col in df.columns:
        if col not in string_cols:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def SMN1_collection(SMAca_output_file, SMAca_thresholds):
    '''

    :param SMAca_output_file: Results generated by SMAca
    :param SMAca_thresholds: data structure with specific thresholds from SMAca software
    :return:
    '''

    # 1. Read SMAca output file. The first line (header) is split by comma. Each entry is aplit by |

    df_smaca = read_smaca_file(SMAca_output_file)


    smaca_cv_fail_threshold = SMAca_thresholds.cv_fail # Max CV_control for QC warning
    smaca_cv_warn_threshold = SMAca_thresholds.cv_warn # Max CV_control for QC OK
    smaca_low_cov_abs = SMAca_thresholds.low_cov_absolute# Absolute mean coverage (×) below which both SMN1 and SMN2 trigger a low-coverage warning.
    smaca_low_cov_rel = SMAca_thresholds.low_cov_relative # fraction of mean control-gene coverage; if both SMN1 and SMN2 fall below this fraction, low-coverage warning.

    # PENDING: refine coverage threshold (if any)
    return (classify_smaca_sample(df_smaca, smaca_cv_fail_threshold, smaca_cv_warn_threshold, smaca_low_cov_abs, smaca_low_cov_rel ))