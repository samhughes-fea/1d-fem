"""
Generate a LaTeX table and PDF from gci_richardson_abaqus_deflection_rotation.csv.

Outputs:
- gci_richardson_abaqus_deflection_rotation_table.tex (standalone compilable document)
- gci_richardson_abaqus_deflection_rotation_table.pdf (after running pdflatex)

Run from repo root or from this directory:
  python post_processing/validation_visualisers/deflection_tables/csv_to_latex_table.py
  cd post_processing/validation_visualisers/output && pdflatex gci_richardson_abaqus_deflection_rotation_table.tex
"""
from pathlib import Path
import csv


def float_fmt(x: str, decimals: int = 4) -> str:
    """Format a numeric CSV value for LaTeX; handle nan/inf."""
    x = x.strip().lower()
    if x in ("nan", "") or (x.startswith("-") and x[1:].strip().lower() == "nan"):
        return "---"
    try:
        v = float(x)
        if abs(v) >= 1e4 or (abs(v) < 1e-3 and v != 0):
            return f"{v:.4e}"
        return f"{v:.{decimals}f}"
    except ValueError:
        return x


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    output_dir = script_dir.parent / "output"
    csv_path = output_dir / "gci_richardson_abaqus_deflection_rotation.csv"
    tex_path = output_dir / "gci_richardson_abaqus_deflection_rotation_table.tex"

    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    rows: list[list[str]] = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if not row or all(cell.strip() == "" for cell in row):
                continue
            rows.append(row)

    # Column indices
    col = {h: i for i, h in enumerate(header)}
    job_id_i = col["job_id"]
    load_type_i = col["load_type"]
    qoi_i = col["QoI"]
    n_fine_i = col["n_fine"]
    phi_fine_i = col["phi_fine"]
    phi_ext_i = col["phi_ext"]
    abaqus_i = col["Abaqus"]
    err_ext_i = col["error_ext_Abaqus_pct"]
    gci_12_i = col["GCI_12_pct"]
    gci_23_i = col["GCI_23_pct"]
    p_obs_i = col["p_obs"]

    # Build table body
    lines: list[str] = []
    for r in rows:
        job_id = r[job_id_i]
        load_type = r[load_type_i]
        qoi = r[qoi_i].replace("_", " ")
        n_fine = r[n_fine_i]
        phi_fine = float_fmt(r[phi_fine_i], 4)
        phi_ext = float_fmt(r[phi_ext_i], 4)
        abaqus = float_fmt(r[abaqus_i], 4)
        err_ext = float_fmt(r[err_ext_i], 3)
        gci_12 = float_fmt(r[gci_12_i], 4)
        gci_23 = float_fmt(r[gci_23_i], 4)
        p_obs = float_fmt(r[p_obs_i], 2)
        lines.append(
            f"    {job_id} & {load_type} & {qoi} & {n_fine} & {p_obs} & "
            f"{phi_fine} & {phi_ext} & {abaqus} & {err_ext} & {gci_12} & {gci_23} \\\\"
        )
    table_body = "\n".join(lines)

    tex_content = r"""\documentclass[11pt]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{booktabs}
\usepackage{array}
\usepackage{graphicx}
\usepackage{geometry}
\geometry{margin=1in}

\title{GCI--Richardson vs Abaqus: Deflection and Rotation}
\date{}

\begin{document}
\maketitle

\begin{table}[htbp]
  \centering
  \small
  \resizebox{\textwidth}{!}{%
  \begin{tabular}{@{}llllccccccc@{}}
    \toprule
    Job & Load type & QoI & $n_{\mathrm{fine}}$ & $p_{\mathrm{obs}}$ &
    $\phi_{\mathrm{fine}}$ & $\phi_{\mathrm{ext}}$ & Abaqus &
    Error ext.\ (\%) & GCI$_{12}$ (\%) & GCI$_{23}$ (\%) \\
    \midrule
""" + table_body + r"""
    \bottomrule
  \end{tabular}%
  }
  \caption{GCI--Richardson extrapolation vs Abaqus reference (tip deflection in mm, tip rotation in deg).}
  \label{tab:gci-richardson-abaqus}
\end{table}

\noindent\textbf{Notes:}
\begin{itemize}
  \item \textbf{Tip rotation:} Abaqus reference (tip rotation) is zero when rotational DOF (UR) was not written to the ODB. Request field output \texttt{U} and \texttt{UR} in the Abaqus step and re-run extraction to obtain rotation comparison.
  \item \textbf{Distributed loads (Triangular, Parabolic):} Rows with Abaqus deflection zero indicate missing or non-matching Abaqus reference (e.g.\ job not run, or tip node/ordering differs). Error vs.\ Abaqus is then not meaningful for those rows; only point-load and UDL rows have a non-zero Abaqus reference in the current data.
\end{itemize}

\end{document}
"""

    tex_path.write_text(tex_content, encoding="utf-8")
    print(f"Wrote {tex_path}")


if __name__ == "__main__":
    main()
