"""
CorrectTimeStampOphirVega.py

Opravuje kvantovanou časovou osu v log souborech Ophir StarLab / Vega PD300-IR.
- Vstupní soubory v adresáři `data` zůstávají beze změny.
- Opravené soubory jsou uloženy do adresáře `data_correctedTimeStamp`.
- Skutečná vzorkovací frekvence je 15 Hz (dt = 1/15 s), Timestamp na 6 desetinných míst.
- Sloupec Channel A i zbytek souboru zůstávají byte-for-byte beze změny,
  pouze do řádku `;Notes:` se připíše poznámka.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np


DEFAULT_INPUT_DIR = Path("data")
DEFAULT_OUTPUT_DIR = Path("data_correctedTimeStamp")
DEFAULT_FS_HZ = 15.0
NOTE_TEXT = "21112026a Corrected timestamps, sample frequency 15 Hz"
MAX_ALLOWED_RESIDUAL_S = 0.025  # 25 ms pojistka proti skutečnému výpadku vzorků


def detect_eol(line: str) -> str:
    if line.endswith("\r\n"):
        return "\r\n"
    if line.endswith("\n"):
        return "\n"
    return ""


def check_for_dropouts(raw_timestamps: list[float], file_name: str) -> float:
    """Ověří, že v původní časové ose není skutečný výpadek vzorků."""
    t = np.asarray(raw_timestamps, dtype=np.float64)
    idx = np.arange(len(t), dtype=np.float64)
    slope, intercept = np.polyfit(idx, t, 1)
    max_res = float(np.max(np.abs(t - (slope * idx + intercept))))
    if max_res > MAX_ALLOWED_RESIDUAL_S:
        print(
            f"[VAROVÁNÍ] {file_name}: maximální reziduum časové osy = "
            f"{max_res * 1000:.2f} ms (> {MAX_ALLOWED_RESIDUAL_S * 1000:.0f} ms) – "
            "možný výpadek vzorku!"
        )
    return max_res


def correct_file(src_path: Path, dst_path: Path, fs_hz: float = DEFAULT_FS_HZ) -> tuple[int, float]:
    """Opraví časovou řadu v jednom souboru a uloží výsledek do `dst_path`."""
    raw_bytes = src_path.read_bytes()
    text = raw_bytes.decode("utf-8", errors="surrogateescape")
    lines = text.splitlines(keepends=True)

    out_lines: list[str] = []
    raw_timestamps: list[float] = []
    in_data = False
    sample_idx = 0

    for line in lines:
        if not in_data:
            if line.startswith(";Notes:"):
                eol = detect_eol(line)
                existing_note = line[len(";Notes:") :].strip()
                if existing_note:
                    new_line = f";Notes: {existing_note}; {NOTE_TEXT}{eol}"
                else:
                    new_line = f";Notes: {NOTE_TEXT}{eol}"
                out_lines.append(new_line)
                continue

            out_lines.append(line)
            if "Timestamp" in line and "Channel A" in line:
                in_data = True
            continue

        # Datová část: nahradíme pouze první sloupec (Timestamp před prvním tabulátorem)
        if "\t" not in line or not line.strip():
            out_lines.append(line)
            continue

        ts_field, rest = line.split("\t", 1)
        raw_timestamps.append(float(ts_field.strip()))

        corrected_t = sample_idx / fs_hz
        width = max(len(ts_field), 14)
        new_ts_field = f"{corrected_t:{width}.6f}"
        out_lines.append(f"{new_ts_field}\t{rest}")
        sample_idx += 1

    if sample_idx == 0:
        raise ValueError(f"V souboru {src_path} nebyla nalezena žádná datová řada.")

    max_res = check_for_dropouts(raw_timestamps, src_path.name)

    dst_path.parent.mkdir(parents=True, exist_ok=True)
    dst_path.write_bytes("".join(out_lines).encode("utf-8", errors="surrogateescape"))
    return sample_idx, max_res


def process_directory(
    input_dir: Path = DEFAULT_INPUT_DIR,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    fs_hz: float = DEFAULT_FS_HZ,
) -> None:
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Vstupní adresář neexistuje: {input_dir}")

    files = sorted(input_dir.glob("*.txt"))
    if not files:
        print(f"Ve vstupním adresáři {input_dir} nebyly nalezeny žádné soubory .txt.")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    max_res_all = 0.0
    total_samples = 0
    for src in files:
        dst = output_dir / src.name
        n_samples, max_res = correct_file(src, dst, fs_hz=fs_hz)
        total_samples += n_samples
        max_res_all = max(max_res_all, max_res)

    print(
        f"Hotovo: zpracováno {len(files)} souborů ({total_samples} vzorků celkem) -> {output_dir}\n"
        f"Vzorkovací frekvence: {fs_hz:.6f} Hz (dt = {1.0 / fs_hz:.6f} s)\n"
        f"Max. reziduum původní časové osy napříč soubory: {max_res_all * 1000:.2f} ms"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Oprava časových značek v log souborech Ophir StarLab / Vega PD300-IR."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Vstupní adresář se zdrojovými soubory (výchozí: data)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Výstupní adresář pro opravené soubory (výchozí: data_correctedTimeStamp)",
    )
    parser.add_argument(
        "--fs",
        type=float,
        default=DEFAULT_FS_HZ,
        help="Vzorkovací frekvence v Hz (výchozí: 15.0)",
    )
    args = parser.parse_args()
    process_directory(input_dir=args.input_dir, output_dir=args.output_dir, fs_hz=args.fs)


if __name__ == "__main__":
    main()
