"""
Generate a LaTeX table and PDF from gci_richardson_abaqus_deflection_rotation.csv.

Outputs (in grid_convergence_study/gci_tables/):
- gci_richardson_abaqus_deflection_rotation_table.tex (standalone compilable document)
- gci_richardson_abaqus_deflection_rotation_table.pdf (built automatically if pdflatex is available)

Run from repo root:
  python post_processing/validation_visualisers/grid_convergence_study/csv_to_latex_table.py
"""
from pathlib import Path
import csv
import subprocess
import sys

# For ABAQUS_REFERENCE_N in caption/notes
_SCRIPT_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPT_DIR.parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
from post_processing.validation_visualisers.abaqus.config import ABAQUS_REFERENCE_N


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
    validation_dir = script_dir.parent
    output_dir = validation_dir / "grid_convergence_study" / "gci_tables"
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

    # QoI display: match Roark style ($u_y$ (mm), $\theta_z$ (deg))
    qoi_tex = {
        "tip_deflection_mm": r"$u_y$ (mm)",
        "tip_rotation_deg": r"$\theta_z$ (deg)",
    }

    # Build table body
    lines: list[str] = []
    for r in rows:
        job_id = r[job_id_i]
        load_type = r[load_type_i]
        qoi_raw = r[qoi_i].strip()
        qoi = qoi_tex.get(qoi_raw, qoi_raw.replace("_", " "))
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

    n_ref = ABAQUS_REFERENCE_N
    tex_content = r"""\documentclass[11pt]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{booktabs}
\usepackage{amsmath}
\usepackage{array}
\usepackage{graphicx}
\usepackage[margin=1in]{geometry}

\begin{document}

\begin{table}[htbp]
\centering
\caption{GCI--Richardson extrapolation vs.\ Abaqus ($n=""" + str(n_ref) + r"""$): deflection and rotation.}
\label{tab:gci-richardson-abaqus}
\resizebox{\textwidth}{!}{%
\begin{tabular}{@{}cll r rr r r rr rr@{}}
\toprule
Job & Load type & QoI & $n_{\mathrm{fine}}$ & $p_{\mathrm{obs}}$ &
$\phi_{\mathrm{fine}}$ & $\phi_{\mathrm{ext}}$ & Abaqus ($n=""" + str(n_ref) + r"""$) &
Error ext.\ (\%) & GCI$_{12}$ (\%) & GCI$_{23}$ (\%) \\
\midrule
""" + table_body + r"""
\bottomrule
\end{tabular}%
}
\medskip
\raggedright
\small
\textbf{Note.} FEM grids: $n_1=128$, $n_2=64$, $n_3=32$ elements (refinement ratio $r=2$); safety factor $F_s=1.25$. Abaqus reference: """ + str(n_ref) + r"""-element mesh (converged benchmark). $\phi_{\mathrm{fine}}$ = fine-grid FEM solution; $\phi_{\mathrm{ext}}$ = Richardson extrapolate. Error ext.\ = $100(\phi_{\mathrm{ext}}-\mathrm{Abaqus})/|\mathrm{Abaqus}|$. GCI$_{12}$, GCI$_{23}$ from three-grid formula. Tip deflection in mm, tip rotation in deg. Abaqus tip rotation is zero if UR was not in the ODB.
\end{table}

\end{document}
"""

    tex_path.write_text(tex_content, encoding="utf-8")
    print(f"Wrote {tex_path}")

    # Build PDF in gci_tables/ if pdflatex is available (run twice for refs)
    pdf_path = output_dir / "gci_richardson_abaqus_deflection_rotation_table.pdf"
    try:
        for _ in range(2):
            result = subprocess.run(
                ["pdflatex", "-interaction=nonstopmode", tex_path.name],
                cwd=str(output_dir),
                capture_output=True,
                timeout=60,
                text=True,
            )
            if result.returncode != 0:
                print("pdflatex failed:", result.stderr[-500:] if result.stderr else result.returncode)
                break
        else:
            if pdf_path.is_file():
                print(f"Wrote {pdf_path}")
            else:
                print("pdflatex succeeded but PDF not found.")
    except FileNotFoundError:
        print("pdflatex not found; PDF not built. Run manually: cd gci_tables && pdflatex gci_richardson_abaqus_deflection_rotation_table.tex")
    except subprocess.TimeoutExpired:
        print("pdflatex timed out; PDF not built.")


if __name__ == "__main__":
    main()
    sys.exit(0)
