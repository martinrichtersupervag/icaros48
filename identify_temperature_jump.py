"""
identify_temperature_jump.py
============================
Samostatný detektor a analyzátor teplotního skoku, který využívá výpočetní
knihovnu `signal_analyzer.py`, načte libovolný datový soubor (bez jakéhokoliv
porovnávání s ostatními soubory v adresáři `data/`) a čistě z absolutních
hodnot po sobě jdoucích rozdílů |y_{i+1} - y_i| v intervalech:
    10 s, 20 s, 60 s a 120 s od detekovaného počátku děje t0
vypočítá velikost teplotního skoku |Delta T| = |T - 24 °C| a úroveň
detekovaného signálu.

Výsledek vypíše v přehledné tabulce do konzole a automaticky uloží do souboru:
    <input filename>-analysis.txt

Použití:
  python identify_temperature_jump.py "data/0 C.txt"
  python identify_temperature_jump.py "data/2 C.txt" "data/24 C.txt" "data/38 C.txt"
  python identify_temperature_jump.py --all
"""

from __future__ import annotations

import os
import re
import sys
import glob
import argparse
from typing import Dict, List, Tuple, Any

import numpy as np

from signal_analyzer import SignalAnalyzer, SystemConfig

WINDOWS_SEC: Tuple[int, ...] = (10, 20, 60, 120)
REF_TEMP_C: float = 24.0


class AbsoluteSignalJumpDetector:
    """
    Samostatný detektor velikosti teplotního skoku |Delta T| = |T - 24 °C|
    založený výhradně na absolutních změnách signálu S_w = sum(|y_{i+1} - y_i|)
    a rychlosti změn v časových oknech 10 s, 20 s, 60 s a 120 s od počátku děje t0.
    Neporovnává soubor s žádnou databází ostatních souborů.
    """

    # Analytické kalibrační koeficienty fyzikálního zákona pro jednotlivá okna w in (10, 20, 60, 120):
    #   |Delta T|_w = clip(c0 + c1 * sqrt(max(0, S_w - S0_w)) + c2 * max(0, S_w - S0_w), 0.0, 24.0)
    # kde S_w [uW] je součet absolutních rozdílů po sobě jdoucích hodnot v okně [t0, t0 + w]
    # a S0_w [uW] je základní klidový práh šumu/driftu při |Delta T| = 0 °C (24 °C).
    WINDOW_FORMULAS: Dict[int, Tuple[float, float, float, float]] = {
        10:  (6.0,  0.0, 2.750, -0.0040),
        20:  (11.0, 0.0, 2.150, -0.0030),
        60:  (24.0, 0.0, 1.550, -0.0018),
        120: (41.0, 0.0, 1.250, -0.0010),
    }

    # Váhy jednotlivých časových oken pro výpočet celkového konsenzu |Delta T|
    WINDOW_WEIGHTS: Dict[int, float] = {
        10:  0.15,
        20:  0.20,
        60:  0.30,
        120: 0.35,
    }

    def __init__(self, fs: float = 16.0):
        self.config = SystemConfig(resample_fs=fs, reference_temp_c=REF_TEMP_C)
        self.analyzer = SignalAnalyzer(self.config)

    @classmethod
    def compute_jump_from_abs_diff(cls, window_s: int, sum_abs_dv_uW: float) -> float:
        """
        Vypočítá velikost teplotního skoku |Delta T| [°C] přímo ze součtu
        absolutních hodnot rozdílů po sobě jdoucích bodů S_w = sum(|y_{i+1} - y_i|) [uW].
        """
        s0, c0, c1, c2 = cls.WINDOW_FORMULAS[window_s]
        eff_s = max(0.0, sum_abs_dv_uW - s0)
        if eff_s <= 0.0:
            return 0.0
        val = c0 + c1 * np.sqrt(eff_s) + c2 * eff_s
        return float(np.clip(val, 0.0, 24.0))

    @staticmethod
    def classify_signal_level(abs_dt: float) -> str:
        """Klasifikace detekovaného signálu a teplotního skoku podle |Delta T|."""
        if abs_dt < 1.5:
            return "Úroveň 0 – Klidový signál / bez teplotního skoku (T ≈ 24 °C)"
        elif abs_dt <= 6.5:
            return "Úroveň 1 – Slabý podnět / malý teplotní skok (|ΔT| ≈ 2–6 °C, pomalý dojezd)"
        elif abs_dt <= 14.5:
            return "Úroveň 2 – Střední teplotní skok (|ΔT| ≈ 7–14 °C, tlumené kmitání)"
        elif abs_dt <= 19.5:
            return "Úroveň 3 – Silný teplotní skok (|ΔT| ≈ 15–19 °C, rychlé vícenásobné kmity)"
        else:
            return "Úroveň 4 – Extrémní teplotní skok (|ΔT| ≥ 20 °C, vysokofrekvenční náběh)"

    def analyze_file(self, filepath: str) -> Dict[str, Any]:
        """
        Načte jeden konkrétní soubor pomocí knihovny `signal_analyzer.py`,
        detekuje počátek t0 a v oknech 10, 20, 60 a 120 s vyhodnotí teplotní skok
        čistě z absolutních změn po sobě jdoucích hodnot.
        """
        analysis = self.analyzer.analyze_temperature_jump_windows(filepath, WINDOWS_SEC)
        finfo = analysis["file_info"]
        win_data = analysis["windows"]

        parsed = self.analyzer.parse_starlab_file(filepath)
        t_raw = parsed["t_raw"]
        v_clean_uW = parsed["v_clean"] * 1e6
        t0 = analysis["onset_t0_s"]

        table_rows: List[Dict[str, Any]] = []
        window_estimates: List[float] = []

        for w in WINDOWS_SEC:
            wm = win_data[w]
            m_raw = (t_raw >= t0) & (t_raw <= t0 + float(w))
            v_win_raw = v_clean_uW[m_raw]
            raw_consecutive_diffs = np.abs(np.diff(v_win_raw)) if len(v_win_raw) > 1 else np.array([0.0])
            sum_raw_uW = float(np.sum(raw_consecutive_diffs))
            max_step_uW = float(np.max(raw_consecutive_diffs))
            mean_step_uW = float(np.mean(raw_consecutive_diffs))

            s_norm_uW = wm["sum_abs_dv_norm_uW"]
            est_abs_dt = self.compute_jump_from_abs_diff(w, s_norm_uW)
            window_estimates.append(est_abs_dt)

            table_rows.append({
                "window_s": w,
                "interval_label": f"0 – {w} s",
                "t_abs_range": f"{t0:.1f} – {t0 + w:.1f} s",
                "n_samples": int(len(v_win_raw)),
                "sum_abs_diff_raw_uW": sum_raw_uW,
                "sum_abs_diff_norm_uW": s_norm_uW,
                "mean_step_diff_uW": mean_step_uW,
                "max_step_diff_uW": max_step_uW,
                "mean_rate_uW_s": wm["mean_rate_uW_s"],
                "max_speed_uW_s": wm["max_speed_uW_s"],
                "fringe_sweeps": wm["fringe_sweeps"],
                "est_abs_delta_T_C": est_abs_dt,
            })

        w_arr = np.array([self.WINDOW_WEIGHTS[w] for w in WINDOWS_SEC])
        consensus_dt = float(np.sum(w_arr * np.array(window_estimates)))
        std_dt = float(np.std(window_estimates))

        return {
            "file_info": finfo,
            "onset_t0_s": t0,
            "initial_y0_uW": analysis["initial_y0_uW"],
            "window_table": table_rows,
            "consensus": {
                "abs_delta_T_C": consensus_dt,
                "uncertainty_std_C": std_dt,
                "signal_level": self.classify_signal_level(consensus_dt),
            },
        }


def format_estimation_table(result: Dict[str, Any]) -> str:
    """Naformátuje tabulku detekce teplotního skoku z absolutních změn pro okna 10, 20, 60 a 120 s."""
    fi = result["file_info"]
    cons = result["consensus"]
    t0 = result["onset_t0_s"]
    y0 = result["initial_y0_uW"]

    lines: List[str] = []
    lines.append("=" * 112)
    lines.append(f"DETEKTOR TEPLOTNÍHO SKOKU Z ABSOLUTNÍCH ZMĚN SIGNÁLU | Soubor: {fi['filename']}")
    lines.append("=" * 112)
    lines.append(f"  Cesta k souboru    : {fi['filepath']}")
    lines.append(
        f"  Datum měření       : {fi['logged_timestamp'] or 'neuvedeno'} | Rozsah senzoru: {fi['sensor_range']} "
        f"(faktor citlivosti: {fi['range_scale_factor']:.2f}x)"
    )
    lines.append(f"  Detekovaný počátek : t0 = {t0:.3f} s | Výchozí klidová hladina y0 = {y0:.2f} uW | Referenční klid: {REF_TEMP_C:.1f} °C")

    if len(fi["glitches"]) > 0:
        g_str = ", ".join(f"{g['time_s']:.2f}..{g['end_time_s']:.2f}s" for g in fi["glitches"])
        lines.append(f"  [!] OPRAVENÝ VÝPADEK SIGNÁLU: Detekován a interpolován výpadek na 0 W v t = {g_str}")

    lines.append("-" * 112)
    header = (
        f"{'Interval':<10} | {'Čas děje':<15} | {'Σ|y_i+1 - y_i|':>14} | {'Σ|Δy| (norm.)':>13} | "
        f"{'Prům. krok |Δy|':>15} | {'Max krok |Δy|':>13} | {'Rychlost |Δy|/s':>15} | "
        f"{'Teplotní skok |ΔT|':>18}"
    )
    units = (
        f"{'od t0 [s]':<10} | {'v záznamu [s]':<15} | {'surový [uW]':>14} | {'[uW_300]':>13} | "
        f"{'[uW/vzorek]':>15} | {'[uW/vzorek]':>13} | {'[uW/s]':>15} | "
        f"{'od 24 °C [°C]':>18}"
    )
    lines.append(header)
    lines.append(units)
    lines.append("-" * 112)

    for row in result["window_table"]:
        lines.append(
            f"{row['interval_label']:<10} | {row['t_abs_range']:<15} | "
            f"{row['sum_abs_diff_raw_uW']:14.2f} | {row['sum_abs_diff_norm_uW']:13.2f} | "
            f"{row['mean_step_diff_uW']:15.3f} | {row['max_step_diff_uW']:13.2f} | "
            f"{row['mean_rate_uW_s']:15.3f} | "
            f"{'±' + format(row['est_abs_delta_T_C'], '4.1f') + ' °C':>18}"
        )

    lines.append("=" * 112)
    return "\n".join(lines)


def get_analysis_output_path(input_filepath: str) -> str:
    """Sestaví výstupní cestu '<input filename>-analysis.txt'."""
    base, ext = os.path.splitext(input_filepath)
    if ext.lower() in (".txt", ".csv", ".tsv"):
        return f"{base}-analysis.txt"
    return f"{input_filepath}-analysis.txt"


def save_and_print_estimation(result: Dict[str, Any], print_to_console: bool = True) -> str:
    """Vypíše tabulku do konzole a uloží ji do souboru '<input filename>-analysis.txt'."""
    report_text = format_estimation_table(result)
    out_path = get_analysis_output_path(result["file_info"]["filepath"])
    with open(out_path, "w", encoding="utf-8") as fp:
        fp.write(report_text + "\n")

    if print_to_console:
        print("\n" + report_text)
        print(f"  -> Výsledek uložen do souboru: {out_path}")
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Detektor teplotního skoku |T - 24 °C| čistě z absolutních změn |y_{i+1} - y_i| v oknech 10, 20, 60 a 120 s."
    )
    parser.add_argument(
        "files",
        nargs="*",
        help="Cesta k jednomu nebo více datovým souborům (např. 'data/0 C.txt').",
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="data",
        help="Výchozí adresář pro hledání souborů (výchozí: 'data').",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Spustit detektor samostatně pro každý soubor v adresáři data/ a uložit <filename>-analysis.txt.",
    )
    args = parser.parse_args()

    detector = AbsoluteSignalJumpDetector()

    target_files: List[str] = []
    if args.all:
        target_files = sorted(
            [
                f for f in glob.glob(os.path.join(args.data_dir, "*.txt"))
                if not f.lower().endswith("-analysis.txt")
            ],
            key=lambda x: int(re.search(r"(\d+)", os.path.basename(x)).group(1))
            if re.search(r"(\d+)", os.path.basename(x)) else 999,
        )
    elif args.files:
        for f in args.files:
            if os.path.exists(f):
                target_files.append(f)
            elif os.path.exists(os.path.join(args.data_dir, f)):
                target_files.append(os.path.join(args.data_dir, f))
            else:
                print(f"[CHYBA] Soubor nebyl nalezen: {f}", file=sys.stderr)
    else:
        default_samples = ["0 C.txt", "1 C.txt", "2 C.txt", "24 C.txt", "38 C.txt", "47 C.txt", "48 C.txt"]
        for s in default_samples:
            p = os.path.join(args.data_dir, s)
            if os.path.exists(p):
                target_files.append(p)

    if args.all:
        print("\n" + "=" * 75)
        print("POROVNÁNÍ ODHADU TEPLOTNÍHO SKOKU |ΔT| V OKNECH 10 s, 20 s, 60 s A 120 s")
        print("=" * 75)
        print(
            f"{'Soubor':<9} | {'t0 [s]':>7} | {'Σ|Δy| 10s':>9} | {'Σ|Δy| 20s':>9} | {'Σ|Δy| 60s':>9} | "
            f"{'Σ|Δy| 120s':>10} | {'|ΔT| 10s':>8} | {'|ΔT| 20s':>8} | {'|ΔT| 60s':>8} | {'|ΔT| 120s':>9}"
        )
        print("-" * 75)
        saved_count = 0
        for fpath in target_files:
            res = detector.analyze_file(fpath)
            save_and_print_estimation(res, print_to_console=False)
            saved_count += 1
            wt = {r["window_s"]: r for r in res["window_table"]}
            print(
                f"{res['file_info']['filename']:<9} | {res['onset_t0_s']:7.2f} | "
                f"{wt[10]['sum_abs_diff_norm_uW']:9.2f} | {wt[20]['sum_abs_diff_norm_uW']:9.2f} | "
                f"{wt[60]['sum_abs_diff_norm_uW']:9.2f} | {wt[120]['sum_abs_diff_norm_uW']:10.2f} | "
                f"{wt[10]['est_abs_delta_T_C']:6.1f}°C | {wt[20]['est_abs_delta_T_C']:6.1f}°C | "
                f"{wt[60]['est_abs_delta_T_C']:6.1f}°C | {wt[120]['est_abs_delta_T_C']:7.1f}°C"
            )
        print("=" * 75)
        print(f"Všechny detailní výstupy ({saved_count}x '<filename>-analysis.txt') byly uloženy.")
    else:
        for fpath in target_files:
            res = detector.analyze_file(fpath)
            save_and_print_estimation(res, print_to_console=True)


if __name__ == "__main__":
    main()
