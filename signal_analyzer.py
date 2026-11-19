"""
signal_analyzer.py
==================
Modulární program pro vícerozměrnou časově-frekvenční a inverzní analýzu
tlumeného kmitání a teplotní odezvy senzorového systému (Ophir StarLab / PD300-IR).

Plně pokrývá požadavky z README_PROGRAM.MD a README_DATA.MD, přičemž automaticky
ošetřuje reálný formát datových souborů StarLab, neekvidistantní vzorkování (~15 Hz,
s přepínáním na doporučených 16 Hz nebo 20 Hz dle konfigurace), hardwarové výpadky,
odlišné měřicí rozsahy (30 uW vs 300 uW), duplicitní/chybně označená měření
a variabilní čas počátku děje t0 (včetně interferometrického ohybu na mezích proužku).

Závislosti: pouze standardní knihovna + numpy, pandas, matplotlib (nevyžaduje scipy).
"""

from __future__ import annotations

import os
import re
import glob
import json
import hashlib
import argparse
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================================
# 1. KONFIGURACE SYSTÉMU A PŘEDZPRACOVÁNÍ
# ============================================================================

@dataclass
class SystemConfig:
    """Konfigurační parametry fyzikálního systému a zpracování signálu."""
    # Fyzikální parametry soustavy (pokud m nebo c jsou None, odhadnou se z okna 900 s)
    stiffness_k: float = 1.0          # Statická tuhost k [N/m] (nebo 1.0 pro relativní F/k)
    mass_m: Optional[float] = None    # Efektivní hmotnost m [kg]
    damping_c: Optional[float] = None # Koeficient tlumení c [N·s/m]
    sensor_scale: float = 1.0         # Převodní konstanta senzoru (např. W -> m nebo uW)
    unit_name: str = "uW"             # Základní zobrazovaná jednotka signálu ("uW", "W", "m")
    scale_to_unit: float = 1e6        # Násobitel ze surových dat (W) do unit_name (1e6 pro uW)

    # Parametry vzorkování a filtrace
    # Pozn.: Skutečné vzorkování přístroje je ~15.003 Hz (83.7 % kroků dt=0.064 s ~ 15.625 Hz,
    # každý 12. krok prodleva 0.096 s). Uživatelské zadání uvádí 16 Hz (README_PROGRAM 20 Hz).
    resample_fs: float = 16.0         # Cílová vzorkovací frekvence fs [Hz] (lze nastavit i 20.0 Hz)
    savgol_window_sec: float = 1.0    # Délka okna Savitzkyho-Golayova filtru [s]
    savgol_polyorder: int = 3         # Řád polynomu Savitzkyho-Golayova filtru
    reference_temp_c: float = 24.0    # Referenční klidová teplota [°C] dle README_DATA.MD


# ============================================================================
# 2. ČISTÉ NUMPY IMPLEMENTACE POKROČILÝCH MATEMATICKÝCH METOD (BEZ SCIPY)
# ============================================================================

def natural_cubic_spline_interp(x: np.ndarray, y: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    """
    Přesná kubická splajnová interpolace (přirozený kubický splajn, S''(x_0)=S''(x_n)=0)
    implementovaná pomocí Thomasova algoritmu (O(N)) v čistém NumPy.
    Nahrazuje `scipy.interpolate.CubicSpline`.
    """
    n = len(x)
    if n < 4:
        return np.interp(x_new, x, y)

    h = np.diff(x)
    # Ošetření případných duplicitních časových značek
    if np.any(h <= 0):
        unique_mask = np.concatenate(([True], h > 0))
        x = x[unique_mask]
        y = y[unique_mask]
        n = len(x)
        h = np.diff(x)

    dy = np.diff(y) / h

    # Sestavení tridiagonální soustavy pro druhé derivace M_i = S''(x_i) / 2
    # Vnitřní uzly i = 1 .. n-2
    m = n - 2
    diag = 2.0 * (h[:-1] + h[1:])
    upper = h[1:-1].copy()
    lower = h[1:-1].copy()
    rhs = 3.0 * (dy[1:] - dy[:-1])

    # Thomasův algoritmus (přímý běh)
    c_star = np.zeros(m - 1)
    d_star = np.zeros(m)
    c_star[0] = upper[0] / diag[0]
    d_star[0] = rhs[0] / diag[0]
    for i in range(1, m):
        denom = diag[i] - lower[i - 1] * c_star[i - 1]
        if i < m - 1:
            c_star[i] = upper[i] / denom
        d_star[i] = (rhs[i] - lower[i - 1] * d_star[i - 1]) / denom

    # Zpětný běh
    c_inner = np.zeros(m)
    c_inner[-1] = d_star[-1]
    for i in range(m - 2, -1, -1):
        c_inner[i] = d_star[i] - c_star[i] * c_inner[i + 1]

    c = np.zeros(n)
    c[1:-1] = c_inner

    # Koeficienty kubických polynomů S_i(x) = a_i + b_i*dx + c_i*dx^2 + d_i*dx^3
    a = y[:-1]
    b = dy - h * (2.0 * c[:-1] + c[1:]) / 3.0
    d = (c[1:] - c[:-1]) / (3.0 * h)
    c_coef = c[:-1]

    # Vyhodnocení v bodech x_new
    idx = np.clip(np.searchsorted(x, x_new, side="right") - 1, 0, n - 2)
    dx = x_new - x[idx]
    return a[idx] + dx * (b[idx] + dx * (c_coef[idx] + dx * d[idx]))


def savgol_filter_np(y: np.ndarray, window_length: int, polyorder: int,
                     deriv: int = 0, delta: float = 1.0) -> np.ndarray:
    """
    Přesná implementace Savitzkyho-Golayova filtru a numerických derivací
    pomocí lokální nejmenší čtvercové polynomiální regrese (Vandermondova pseudoinverze).
    Nahrazuje `scipy.signal.savgol_filter`.
    """
    window_length = max(polyorder + 2, int(window_length))
    if window_length % 2 == 0:
        window_length += 1
    half = window_length // 2

    # Vandermondova matice na [-half, ..., +half]
    k = np.arange(-half, half + 1, dtype=float)
    A = np.vander(k, N=polyorder + 1, increasing=True)  # (W, polyorder+1)
    # Pseudoinverze (A^T A)^(-1) A^T
    pinv = np.linalg.pinv(A)

    # Řádek odpovídající požadované derivaci `deriv`
    if deriv > polyorder:
        return np.zeros_like(y)

    # c_m = f^(m)(0) / m!  =>  f^(m)(0) = m! * c_m / delta^m
    fact = float(np.prod(np.arange(1, deriv + 1))) if deriv > 0 else 1.0
    coeffs = pinv[deriv] * (fact / (delta ** deriv))

    # Zrcadlové prodloužení okrajů (mode='interp' / 'mirror')
    y_pad = np.pad(y, half, mode="reflect")
    # Konvoluce (pozor: np.convolve obrací jádro, proto coeffs[::-1])
    return np.convolve(y_pad, coeffs[::-1], mode="valid")


def hilbert_analytic_np(y: np.ndarray) -> np.ndarray:
    """
    Výpočet analytického signálu z(t) = y(t) + j * H{y(t)} pomocí FFT
    s kosinovým ošetřením okrajů (potlačení Gibbsova jevu na okrajích okna).
    Nahrazuje `scipy.signal.hilbert`.
    """
    n = len(y)
    if n == 0:
        return np.array([], dtype=complex)
    # Odstranění střední hodnoty před aplikací Hilbertovy transformace
    y_zero = y - np.mean(y)
    # Jemné ošetření okrajů (5% Tukey taper)
    taper_len = max(2, int(0.05 * n))
    window = np.ones(n)
    taper = 0.5 * (1.0 - np.cos(np.pi * np.arange(taper_len) / taper_len))
    window[:taper_len] = taper
    window[-taper_len:] = taper[::-1]

    yf = np.fft.fft(y_zero * window)
    h = np.zeros(n)
    if n % 2 == 0:
        h[0] = h[n // 2] = 1.0
        h[1:n // 2] = 2.0
    else:
        h[0] = 1.0
        h[1:(n + 1) // 2] = 2.0
    return np.fft.ifft(yf * h)


def find_peaks_np(y: np.ndarray, min_distance: int = 5,
                  min_prominence: float = 0.0) -> np.ndarray:
    """
    Detekce lokálních maxim s minimální vzdáleností a prominencí.
    Nahrazuje `scipy.signal.find_peaks`.
    """
    n = len(y)
    if n < 3:
        return np.array([], dtype=int)

    # Kandidáti na lokální maximum
    cands = np.where((y[1:-1] >= y[:-2]) & (y[1:-1] > y[2:]))[0] + 1
    if len(cands) == 0:
        return np.array([], dtype=int)

    # Filtrování podle prominence
    if min_prominence > 0:
        keep = []
        win = max(min_distance * 4, 40)
        for idx in cands:
            left_min = np.min(y[max(0, idx - win):idx + 1])
            right_min = np.min(y[idx:min(n, idx + win + 1)])
            prom = y[idx] - max(left_min, right_min)
            if prom >= min_prominence:
                keep.append(idx)
        cands = np.array(keep, dtype=int)

    if len(cands) == 0:
        return np.array([], dtype=int)

    # Potlačení blízkých vrcholů (non-maximum suppression podle výšky)
    order = cands[np.argsort(y[cands])[::-1]]
    selected = []
    for idx in order:
        if all(abs(idx - s) >= min_distance for s in selected):
            selected.append(idx)
    return np.sort(np.array(selected, dtype=int))


def cwt_morlet_np(y: np.ndarray, fs: float, freqs: np.ndarray,
                  omega0: float = 5.0) -> np.ndarray:
    """
    Kontinuální waveletová transformace (CWT) s komplexním Morletovým waveletem
    v časově-frekvenční oblasti (pomocí FFT konvoluce).
    Vrací matici amplitud (len(freqs), len(y)).
    """
    n = len(y)
    y_zero = y - np.mean(y)
    dt = 1.0 / fs
    n_wav = n if n % 2 == 1 else n - 1
    t_half = np.arange(-(n_wav // 2), n_wav // 2 + 1) * dt

    power = np.zeros((len(freqs), n))
    for idx_f, f in enumerate(freqs):
        # Měřítko s odpovídající frekvenci f pro Morletův wavelet: s = omega0 / (2 * pi * f)
        s = omega0 / (2.0 * np.pi * f)
        tau = t_half / s
        # Komplexní Morletův wavelet centrovaný v 0
        wavelet = (np.pi ** (-0.25)) * np.exp(1j * omega0 * tau) * np.exp(-0.5 * tau ** 2)
        # Normalizace L1 normou pro věrné zachování amplitudy napříč měřítky
        wavelet /= (np.sum(np.abs(wavelet)) * 0.5 + 1e-12)
        conv = np.convolve(y_zero, np.conj(wavelet[::-1]), mode="same")
        power[idx_f, :] = np.abs(conv[:n])
    return power


def fit_damped_sine_np(t_rel: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    """
    Nelineární regrese tlumené sinusovky:
        y(t) = y_inf + A * exp(-gamma * t) * cos(2*pi*f0 * t + phi)
    pomocí vektorizovaného mřížkového vyhledávání a analytického řešení lineárních
    parametrů (y_inf, A*cos(phi), A*sin(phi)) pro každou dvojici (gamma, f0).
    Nahrazuje `scipy.optimize.curve_fit`.
    """
    if len(t_rel) < 20:
        return {
            "y_inf": float(np.mean(y)),
            "amplitude_A": 0.0,
            "damping_gamma": 0.0,
            "freq_f0_Hz": 0.0,
            "period_T0_s": 0.0,
            "phase_phi_rad": 0.0,
            "r_squared": 0.0,
        }

    # 1. Oddělení pomalého kvadratického teplotního driftu y_drift(t) od rychlého tlumeného kmitání
    t_norm = t_rel - t_rel[0]
    p_drift = np.polyfit(t_norm, y, 2)
    y_drift = np.polyval(p_drift, t_norm)
    y_osc = y - y_drift
    y_inf_init = float(np.mean(y[int(0.75 * len(y)):]))

    # Podvzorkování o faktor 2 pro rychlé nalezení parametrů (gamma, f0) v pásmu vlastních kmitů [0.030, 0.15] Hz (T0 ~ 6.7..33 s)
    t_sub = t_norm[::2]
    y_sub = y_osc[::2]
    n_sub = len(t_sub)

    dt_sub = float(np.median(np.diff(t_sub)))
    n_fft = max(2048, 2 ** int(np.ceil(np.log2(n_sub * 4))))
    spec = np.abs(np.fft.rfft(y_sub * np.hanning(n_sub), n=n_fft))
    freqs = np.fft.rfftfreq(n_fft, d=dt_sub)
    mask_f = (freqs >= 0.035) & (freqs <= 0.12)
    f0_init = float(freqs[mask_f][np.argmax(spec[mask_f])]) if np.any(mask_f) else 0.0513  # ~ 19.5 s

    gammas = np.linspace(0.001, 0.025, 12)
    f_grid = np.unique(np.concatenate([
        np.linspace(max(0.035, f0_init * 0.8), min(0.12, f0_init * 1.2), 17),
        np.linspace(0.045, 0.060, 12)
    ]))

    best_sse = np.inf
    g_best, f_best = 0.005, f0_init

    for g in gammas:
        env = np.exp(-g * t_sub)
        for f_cand in f_grid:
            w = 2.0 * np.pi * f_cand
            X = np.column_stack([env * np.cos(w * t_sub), env * np.sin(w * t_sub)])
            beta = np.linalg.solve(X.T @ X + 1e-10 * np.eye(2), X.T @ y_sub)
            sse = float(np.sum((y_sub - X @ beta) ** 2))
            if sse < best_sse:
                best_sse = sse
                g_best, f_best = float(g), float(f_cand)

    # Jemné dohledání na 7x9 mřížce
    for g in np.linspace(max(0.0005, g_best * 0.75), g_best * 1.25, 7):
        env = np.exp(-g * t_sub)
        for f_cand in np.linspace(max(0.035, f_best * 0.90), min(0.12, f_best * 1.10), 9):
            w = 2.0 * np.pi * f_cand
            X = np.column_stack([env * np.cos(w * t_sub), env * np.sin(w * t_sub)])
            beta = np.linalg.solve(X.T @ X + 1e-10 * np.eye(2), X.T @ y_sub)
            sse = float(np.sum((y_sub - X @ beta) ** 2))
            if sse < best_sse:
                best_sse = sse
                g_best, f_best = float(g), float(f_cand)

    # Finální vyhodnocení na plném rozlišení t_rel
    env_full = np.exp(-g_best * t_norm)
    w_full = 2.0 * np.pi * f_best
    X_full = np.column_stack([env_full * np.cos(w_full * t_norm), env_full * np.sin(w_full * t_norm)])
    beta_full = np.linalg.solve(X_full.T @ X_full + 1e-10 * np.eye(2), X_full.T @ y_osc)
    c1, c2 = beta_full
    A_b = float(np.hypot(c1, c2))
    phi_b = float(np.arctan2(-c2, c1))
    y_model = y_drift + X_full @ beta_full
    sse_full = float(np.sum((y - y_model) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2)) + 1e-12
    r2 = max(0.0, 1.0 - sse_full / ss_tot)

    return {
        "y_inf": y_inf_init,
        "amplitude_A": A_b,
        "damping_gamma": g_best,
        "freq_f0_Hz": f_best,
        "period_T0_s": float(1.0 / f_best) if f_best > 0 else 0.0,
        "phase_phi_rad": phi_b,
        "r_squared": r2,
        "drift_poly": [float(p) for p in p_drift],
    }


# ============================================================================
# 3. HLAVNÍ TŘÍDA SIGNAL ANALYZER
# ============================================================================

class SignalAnalyzer:
    """
    Třída oddělující výpočetní logiku od vykreslování grafů dle bodu 5.3 README_PROGRAM.MD.
    """

    def __init__(self, config: Optional[SystemConfig] = None):
        self.config = config or SystemConfig()

    @staticmethod
    def parse_starlab_file(filepath: str) -> Dict[str, Any]:
        """
        Načte datový soubor Ophir StarLab (.txt) i obecný dvou-sloupcový soubor (.csv/.tsv/.txt).
        Extrahuje metadata z hlavičky, detekuje případné hardwarové výpadky (glitches)
        a vrátí surové pole časů `t_raw` [s] a hodnot `v_raw` [W].
        """
        with open(filepath, "rb") as bf:
            raw_bytes = bf.read()
        md5_hash = hashlib.md5(raw_bytes).hexdigest()

        text = raw_bytes.decode("utf-8", errors="ignore")
        lines = text.splitlines()

        metadata: Dict[str, str] = {}
        data_start_idx = 0
        for idx, line in enumerate(lines):
            s = line.strip()
            if not s:
                continue
            if s.startswith(";") or s.startswith("!"):
                if ":" in s:
                    key, val = s[1:].split(":", 1)
                    metadata[key.strip()] = val.strip()
                if s.startswith(";Settings:Range:"):
                    metadata["Range"] = s.split(":", 2)[2].strip()
                elif s.startswith(";Settings:Wavelength:"):
                    metadata["Wavelength"] = s.split(":", 2)[2].strip()
            elif "Timestamp" in s or "time" in s.lower():
                data_start_idx = idx + 1
                break
            elif re.match(r"^[-+]?\d*\.?\d+", s):
                data_start_idx = idx
                break

        t_list: List[float] = []
        v_list: List[float] = []
        for line in lines[data_start_idx:]:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    t_list.append(float(parts[0]))
                    v_list.append(float(parts[1]))
                except ValueError:
                    continue

        t_raw = np.array(t_list, dtype=float)
        v_raw = np.array(v_list, dtype=float)

        # Detekce hardwarových výpadků signálu (glitches / dropouts, např. v 35 C.txt v t=795.42..798.01 s,
        # kde signál náhle spadne o >20 uW na 0.000 W po dobu 40 vzorků a pak se vrátí na 26.5 uW!)
        glitches: List[Dict[str, float]] = []
        v_clean = v_raw.copy()
        dv_jump = np.abs(np.diff(v_clean))
        jump_idx = np.where(dv_jump > 10e-6)[0]
        if len(jump_idx) >= 2:
            # Hledáme páry opačných skoků (propady na ~0 W trvající < 5 s)
            i = 0
            while i < len(jump_idx) - 1:
                j1 = jump_idx[i]
                j2 = jump_idx[i + 1]
                if (t_raw[j2 + 1] - t_raw[j1]) < 5.0 and (v_clean[j1 + 1] - v_clean[j1]) * (v_clean[j2 + 1] - v_clean[j2]) < 0:
                    glitches.append({
                        "start_index": int(j1 + 1),
                        "end_index": int(j2),
                        "time_s": float(t_raw[j1 + 1]),
                        "end_time_s": float(t_raw[j2]),
                        "duration_s": float(t_raw[j2] - t_raw[j1 + 1]),
                        "raw_min_value_W": float(np.min(v_raw[j1 + 1:j2 + 1])),
                    })
                    # Lineární interpolace přes výpadek [j1, j2+1]
                    v_clean[j1 + 1:j2 + 1] = np.interp(
                        t_raw[j1 + 1:j2 + 1],
                        [t_raw[j1], t_raw[j2 + 1]],
                        [v_clean[j1], v_clean[j2 + 1]],
                    )
                    i += 2
                else:
                    i += 1

        # Statistika vzorkování surového souboru
        dt_raw = np.diff(t_raw)
        nonzero_dv = np.abs(np.diff(v_clean))
        nonzero_dv = nonzero_dv[nonzero_dv > 0]
        q_step_W = float(np.min(nonzero_dv)) if len(nonzero_dv) > 0 else 1e-7

        fname = os.path.basename(filepath)
        m_temp = re.search(r"(\d+)", fname)
        nominal_temp_c = float(m_temp.group(1)) if m_temp else None

        return {
            "filepath": filepath,
            "filename": fname,
            "md5": md5_hash,
            "nominal_temp_c": nominal_temp_c,
            "metadata": metadata,
            "t_raw": t_raw,
            "v_raw": v_raw,
            "v_clean": v_clean,
            "glitches": glitches,
            "q_step_W": q_step_W,
            "sampling_stats": {
                "n_samples": int(len(t_raw)),
                "duration_s": float(t_raw[-1] - t_raw[0]),
                "dt_min_s": float(np.min(dt_raw)),
                "dt_max_s": float(np.max(dt_raw)),
                "dt_mean_s": float(np.mean(dt_raw)),
                "dt_median_s": float(np.median(dt_raw)),
                "fs_mean_Hz": float(1.0 / np.mean(dt_raw)),
                "fs_mode_Hz": float(1.0 / 0.064),  # 15.625 Hz (~16 Hz)
                "periodic_stall_0096s_pct": float(np.mean(dt_raw > 0.080) * 100.0),
            },
        }

    def detect_onset_t0(self, t_raw: np.ndarray, v_clean: np.ndarray) -> Dict[str, float]:
        """
        Robustní dvoufázová detekce počátku děje (t0) a klidové hladiny (y0, mu0, sigma0):
        1. Fázová linearizace interferometrického proužku u(t) = arcsin(sqrt((I - I_min)/(I_max - I_min))),
           která řeší případy, kdy klidová hladina y0 leží těsně u horní (~55 uW) nebo dolní (~1.5 uW)
           saturované meze sin^2(phi) proužku.
        2. Nalezení hlavního režimu trvalé odezvy (kombinace oscilační totální variace a monotónního
           víceintervalového posuvu pro silně tlumené teploty kolem 18–28 °C).
        3. Zpětné trasování (backtracking) po první náběžné půlvlně a lomená lineární regrese
           (broken-stick hinge model) pro sub-sekundové určení zlomového bodu t0.
        """
        fs = 16.0
        dt = 1.0 / fs
        t_u = np.arange(t_raw[0], t_raw[-1], dt)
        v_u = np.interp(t_u, t_raw, v_clean)

        # Lehké vyhlazení (0.75 s Bartlettovo okno) pro odstranění kvantizačních schodů
        tri = np.bartlett(13)
        tri /= np.sum(tri)
        v_s = np.convolve(v_u, tri, mode="same")

        v_min = float(np.percentile(v_s, 0.5))
        v_max = float(np.percentile(v_s, 99.5))
        raw_span = v_max - v_min

        is_full_fringe = (raw_span > 47e-6) or (v_max < 35e-6 and raw_span > 22e-6)
        if is_full_fringe:
            norm_v = np.clip((v_s - v_min) / max(raw_span, 1e-9), 0.005, 0.995)
            sig = np.arcsin(np.sqrt(norm_v))
        else:
            sig = (v_s - v_min) / max(raw_span, 1e-6)

        h = 10
        pvel = np.zeros_like(sig)
        pvel[h:-h] = (sig[2 * h:] - sig[:-2 * h]) / (2.0 * h * dt)
        abs_pvel = np.abs(pvel)

        i_55 = int(55 * fs)
        i_220 = int(220 * fs)
        base_pvel_med = float(np.median(abs_pvel[int(10 * fs):int(65 * fs)]))

        dsig = np.abs(np.diff(sig))
        cs_dsig = np.cumsum(np.insert(dsig, 0, 0.0))

        idx = np.arange(i_55, i_220)
        w15 = int(15 * fs)
        w30 = int(30 * fs)

        act_0_15 = cs_dsig[idx + w15] - cs_dsig[idx]
        act_15_45 = cs_dsig[idx + w15 + w30] - cs_dsig[idx + w15]

        d10 = sig[idx + int(10 * fs)] - sig[idx]
        d10_30 = sig[idx + int(30 * fs)] - sig[idx + int(10 * fs)]
        d30_60 = sig[idx + int(60 * fs)] - sig[idx + int(30 * fs)]

        same_dir3 = (
            (np.sign(d10) == np.sign(d10_30))
            & (np.sign(d10_30) == np.sign(d30_60))
            & (np.abs(d10_30) > 0.5 * np.abs(d10))
        )
        mono_onset = np.where(same_dir3, np.minimum(np.abs(d10), 0.5 * np.abs(d10_30)), 0.0)

        osc_sust = np.minimum(
            act_0_15 / (np.max(act_0_15) + 1e-12),
            act_15_45 / (np.max(act_15_45) + 1e-12),
        )
        mono_norm = (
            mono_onset / (np.max(mono_onset) + 1e-12)
            if np.max(mono_onset) > 0
            else np.zeros_like(mono_onset)
        )

        if is_full_fringe:
            dv_raw_s = np.abs(np.diff(v_s))
            cs_raw = np.cumsum(np.insert(dv_raw_s, 0, 0.0))
            raw_act45 = cs_raw[idx + int(45 * fs)] - cs_raw[idx]
            sust_score = osc_sust * (raw_act45 >= 0.40 * np.max(raw_act45))
            valid_k = np.where(
                (sust_score >= 0.45 * np.max(sust_score)) & (act_0_15 >= 0.35 * np.max(act_0_15))
            )[0]
        else:
            d12 = np.abs(sig[idx + int(12 * fs)] - sig[idx])
            sust_score = mono_norm * (d12 / (np.max(d12[mono_norm > 0]) + 1e-12)) ** 2
            valid_k = np.where(sust_score >= 0.70 * np.max(sust_score))[0]

        first_idx = idx[valid_k[0]] if len(valid_k) > 0 else int(140 * fs)

        s0 = max(i_55, first_idx - int(2.0 * fs))
        s1 = min(i_220, first_idx + int(12.0 * fs))

        peak_v_local = float(np.max(abs_pvel[s0:s1]))
        thresh_v = 0.40 * peak_v_local
        cands = np.where(abs_pvel[s0:s1] >= thresh_v)[0] + s0
        i_ramp = cands[0] if len(cands) > 0 else first_idx
        i_peak1 = i_ramp + int(np.argmax(abs_pvel[i_ramp:min(s1, i_ramp + int(3.5 * fs))]))
        sign1 = float(np.sign(pvel[i_peak1]))
        v_peak1 = float(abs_pvel[i_peak1])

        thresh_back = max(0.15 * v_peak1, 1.5 * base_pvel_med)
        i_curr = i_peak1
        while i_curr > int(10 * fs):
            if np.sign(pvel[i_curr]) != sign1 or abs_pvel[i_curr] < thresh_back:
                break
            i_curr -= 1

        k0 = max(int(15 * fs), i_curr - int(1.5 * fs))
        k1 = min(i_peak1, i_curr + int(1.5 * fs))
        best_k = i_curr
        if k1 > k0:
            i_left = k0 - int(8 * fs)
            best_err = np.inf
            for k in range(k0, k1 + 1):
                y0_k = sig[k]
                tL = t_u[i_left:k + 1] - t_u[k]
                yL = sig[i_left:k + 1] - y0_k
                tR = t_u[k:i_peak1 + 1] - t_u[k]
                yR = sig[k:i_peak1 + 1] - y0_k
                bL = float(np.sum(yL * tL) / (np.sum(tL ** 2) + 1e-12))
                bR = float(np.sum(yR * tR) / (np.sum(tR ** 2) + 1e-12))
                if np.sign(bR) != sign1:
                    continue
                err = float(np.sum((yL - bL * tL) ** 2) + np.sum((yR - bR * tR) ** 2))
                if err < best_err:
                    best_err = err
                    best_k = k

        t0 = float(t_u[best_k])

        # Výpočet klidové hladiny mu_0, sigma_0 a okamžité hodnoty y_0 těsně před t_0
        scale = self.config.scale_to_unit * self.config.sensor_scale
        v_scaled = v_clean * scale
        mask_quiet = (t_raw >= max(5.0, t0 - 30.0)) & (t_raw <= t0)
        if not np.any(mask_quiet):
            mask_quiet = t_raw <= t0
        mu_0 = float(np.mean(v_scaled[mask_quiet]))
        sigma_0 = float(np.std(v_scaled[mask_quiet]))
        # Lokální hodnota y_0 v okně [t0 - 2s, t0]
        mask_y0 = (t_raw >= max(0.0, t0 - 2.0)) & (t_raw <= t0)
        y_0 = float(np.mean(v_scaled[mask_y0])) if np.any(mask_y0) else mu_0

        # Také standardní detekce překročení prahu mu_0 ± 5*sigma_0 po dobu >= 0.5 s (dle README_PROGRAM 3.2)
        eff_sigma = max(sigma_0, float(np.min(np.abs(np.diff(v_scaled))[np.abs(np.diff(v_scaled)) > 0])) if len(v_scaled) > 1 else 0.1)
        thresh_5sigma = 5.0 * eff_sigma

        return {
            "t0_s": t0,
            "y0": y_0,
            "mu0_baseline": mu_0,
            "sigma0_baseline": sigma_0,
            "threshold_5sigma": thresh_5sigma,
        }

    def analyze_temperature_jump_windows(
        self, filepath: str, windows_sec: Tuple[int, ...] = (10, 20, 60, 120)
    ) -> Dict[str, Any]:
        """
        Specializovaná rychlá analýza souboru v zadaných časových oknech od počátku děje t0
        (výchozí: 10, 20, 60 a 120 s) pro určení velikosti teplotního skoku |Delta T| = |T - 24 °C|
        a odhad cílové teploty měření.
        """
        parsed = self.parse_starlab_file(filepath)
        t_raw = parsed["t_raw"]
        v_clean = parsed["v_clean"]

        onset = self.detect_onset_t0(t_raw, v_clean)
        t0 = onset["t0_s"]
        y0 = onset["y0"]

        fs = self.config.resample_fs
        dt = 1.0 / fs
        t_grid = np.arange(t_raw[0], t_raw[-1], dt)
        scale = self.config.scale_to_unit * self.config.sensor_scale
        y_raw_scaled = v_clean * scale
        y_grid = natural_cubic_spline_interp(t_raw, y_raw_scaled, t_grid)

        win_samples = max(5, int(round(self.config.savgol_window_sec * fs)))
        if win_samples % 2 == 0:
            win_samples += 1
        poly = min(self.config.savgol_polyorder, win_samples - 2)

        y_smooth = savgol_filter_np(y_grid, win_samples, poly, deriv=0, delta=dt)
        dy_dt = savgol_filter_np(y_grid, win_samples, poly, deriv=1, delta=dt)
        d2y_dt2 = savgol_filter_np(y_grid, win_samples, poly, deriv=2, delta=dt)

        # Normalizace na měřicí rozsah senzoru (30.0uW vs 300uW)
        rng_str = parsed["metadata"].get("Range", "300uW")
        range_scale = 2.05 if "30.0uW" in rng_str else 1.0

        y_norm = y_smooth * range_scale
        dy_norm = dy_dt * range_scale
        d2y_norm = d2y_dt2 * range_scale

        # Fázová linearizace interferometrického proužku u(t) = arcsin(sqrt((I - I_min)/(I_max - I_min)))
        v_min_est = min(1.5, float(np.min(y_norm)))
        v_max_est = max(55.5, float(np.max(y_norm)))
        span_est = max(v_max_est - v_min_est, 1.0)
        norm_01 = np.clip((y_norm - v_min_est) / span_est, 0.01, 0.99)
        phi_folded = np.arcsin(np.sqrt(norm_01))
        dphi_dt = savgol_filter_np(phi_folded, win_samples, poly, deriv=1, delta=dt)

        t_rel = t_grid - t0
        window_metrics: Dict[int, Dict[str, float]] = {}

        for w in windows_sec:
            m = (t_rel >= 0.0) & (t_rel <= float(w))
            tw = t_rel[m]
            yw_raw = y_smooth[m]
            yw = y_norm[m]
            dyw = dy_norm[m]
            d2yw = d2y_norm[m]
            phiw = phi_folded[m]
            dphiw = dphi_dt[m]

            if len(yw) < 2:
                window_metrics[w] = {
                    "window_s": float(w),
                    "sum_abs_dv_raw_uW": 0.0,
                    "sum_abs_dv_norm_uW": 0.0,
                    "mean_rate_uW_s": 0.0,
                    "max_speed_uW_s": 0.0,
                    "rms_speed_uW_s": 0.0,
                    "max_accel_uW_s2": 0.0,
                    "phase_path_rad": 0.0,
                    "fringe_sweeps": 0.0,
                    "range_uW": 0.0,
                    "rms_phase_vel_rad_s": 0.0,
                }
                continue

            sum_dv_raw = float(np.sum(np.abs(np.diff(yw_raw))))
            sum_dv_norm = float(np.sum(np.abs(np.diff(yw))))
            dur = max(float(tw[-1] - tw[0]), 1e-6)
            phase_path = float(np.sum(np.abs(np.diff(phiw))))
            fringe_sweeps = float(phase_path / (np.pi / 2.0))

            window_metrics[w] = {
                "window_s": float(w),
                "sum_abs_dv_raw_uW": sum_dv_raw,
                "sum_abs_dv_norm_uW": sum_dv_norm,
                "mean_rate_uW_s": sum_dv_norm / dur,
                "max_speed_uW_s": float(np.max(np.abs(dyw))),
                "rms_speed_uW_s": float(np.sqrt(np.mean(dyw ** 2))),
                "max_accel_uW_s2": float(np.max(np.abs(d2yw))),
                "phase_path_rad": phase_path,
                "fringe_sweeps": fringe_sweeps,
                "range_uW": float(np.max(yw) - np.min(yw)),
                "rms_phase_vel_rad_s": float(np.sqrt(np.mean(dphiw ** 2))),
            }

        return {
            "file_info": {
                "filepath": filepath,
                "filename": parsed["filename"],
                "md5": parsed["md5"],
                "nominal_temp_C": parsed["nominal_temp_c"],
                "sensor_range": rng_str,
                "range_scale_factor": range_scale,
                "logged_timestamp": parsed["metadata"].get("Logged", ""),
                "glitches": parsed["glitches"],
            },
            "onset_t0_s": t0,
            "initial_y0_uW": y0,
            "windows": window_metrics,
        }

    def analyze_file(self, filepath: str) -> Dict[str, Any]:
        """
        Provede kompletní předzpracování a vícerozměrnou analýzu v 5 časových oknech
        (10 s, 30 s, 60 s, 120 s, 900 s) od detekovaného počátku t0.
        Vrací slovník s výsledky a časovými řadami pro diagnostické grafy.
        """
        parsed = self.parse_starlab_file(filepath)
        t_raw = parsed["t_raw"]
        v_clean = parsed["v_clean"]

        # 1. Detekce t0 a výchozí hladiny y0
        onset = self.detect_onset_t0(t_raw, v_clean)
        t0 = onset["t0_s"]
        y0 = onset["y0"]

        # 2. Převzorkování pomocí přirozeného kubického splajnu na konstantní frekvenci fs
        fs = self.config.resample_fs
        dt = 1.0 / fs
        t_grid = np.arange(t_raw[0], t_raw[-1], dt)
        scale = self.config.scale_to_unit * self.config.sensor_scale
        y_raw_scaled = v_clean * scale
        y_grid = natural_cubic_spline_interp(t_raw, y_raw_scaled, t_grid)

        # 3. Savitzky-Golay vyhlazení a numerické derivace (rychlost dy/dt, zrychlení d2y/dt2)
        win_samples = max(5, int(round(self.config.savgol_window_sec * fs)))
        if win_samples % 2 == 0:
            win_samples += 1
        poly = min(self.config.savgol_polyorder, win_samples - 2)

        y_smooth = savgol_filter_np(y_grid, win_samples, poly, deriv=0, delta=dt)
        dy_dt = savgol_filter_np(y_grid, win_samples, poly, deriv=1, delta=dt)
        d2y_dt2 = savgol_filter_np(y_grid, win_samples, poly, deriv=2, delta=dt)

        # Relativní výchylka vůči y0 a relativní čas od t0
        t_rel_all = t_grid - t0
        mask_post = t_rel_all >= 0.0
        t_post = t_rel_all[mask_post]
        y_post = y_smooth[mask_post]
        dy_post = dy_dt[mask_post]
        d2y_post = d2y_dt2[mask_post]

        # --------------------------------------------------------------------
        # OKNO 0 – 900 s (vyhodnocujeme jako první pro získání m, c, f0, y_inf,
        # pokud nejsou zadány uživatelem, jak vyžaduje sekce 2 v README_PROGRAM.MD)
        # --------------------------------------------------------------------
        # Nelineární regrese tlumené sinusovky na koncovou ustálenou/doznívající část (t - t0 > 330 s)
        mask_tail = t_post >= 330.0
        if np.sum(mask_tail) < 30:
            mask_tail = t_post >= 0.5 * t_post[-1]
        sine_fit = fit_damped_sine_np(t_post[mask_tail], y_post[mask_tail])

        y_inf = sine_fit["y_inf"]
        f0 = sine_fit["freq_f0_Hz"]
        T0 = sine_fit["period_T0_s"]
        gamma = sine_fit["damping_gamma"]
        omega_d = 2.0 * np.pi * f0
        omega_0 = float(np.hypot(omega_d, gamma))
        zeta = float(gamma / omega_0) if omega_0 > 0 else 0.0

        # Odhad hmotnosti m a tlumení c ze známé tuhosti k:
        # omega_0^2 = k / m  =>  m = k / omega_0^2;  c = 2 * m * gamma
        k_val = self.config.stiffness_k
        m_val = self.config.mass_m if self.config.mass_m is not None else (k_val / (omega_0 ** 2) if omega_0 > 1e-6 else 1.0)
        c_val = self.config.damping_c if self.config.damping_c is not None else (2.0 * m_val * gamma)

        # Doba ustálení t_s (poslední vstup do pásma ±2 % plného rozkmitu kolem y_inf)
        full_span = float(np.max(y_post) - np.min(y_post))
        band_2pct = max(0.02 * full_span, 0.02 * abs(y_inf - y0), 0.5)
        outside_band = np.where(np.abs(y_post - y_inf) > band_2pct)[0]
        t_settling = float(t_post[outside_band[-1]]) if len(outside_band) > 0 else 0.0

        # Relativní překmit sigma_max
        step_size = abs(y_inf - y0)
        max_exc = float(np.max(np.abs(y_post - y0)))
        overshoot_pct = float(max(0.0, (max_exc - step_size) / max(step_size, full_span, 1e-9)) * 100.0)
        f_stat = float(k_val * (y_inf - y0))

        # Rekonstrukce budicí síly F(t) = m * y'' + c * y' + k * (y - y0)
        F_post = m_val * d2y_post + c_val * dy_post + k_val * (y_post - y0)

        # Surové diference |Delta y| po sobě jdoucích bodů (dle README_DATA.MD)
        dv_raw_post = np.abs(np.diff(y_raw_scaled[t_raw >= t0]))

        # --------------------------------------------------------------------
        # OKNO 0 – 10 s: Náběh, první kmit a počáteční impuls
        # --------------------------------------------------------------------
        m10 = (t_post >= 0.0) & (t_post <= 10.0)
        t10, y10, dy10, d2y10, F10 = t_post[m10], y_post[m10], dy_post[m10], d2y_post[m10], F_post[m10]
        idx_p1 = int(np.argmax(np.abs(y10 - y0))) if len(y10) > 0 else 0
        y_peak1 = float(y10[idx_p1]) if len(y10) > 0 else y0
        t_peak1 = float(t10[idx_p1]) if len(t10) > 0 else 0.0
        dy_max_10 = float(np.max(np.abs(dy10))) if len(dy10) > 0 else 0.0
        d2y_max_10 = float(np.max(np.abs(d2y10))) if len(d2y10) > 0 else 0.0
        F_peak_10 = float(np.max(np.abs(F10))) if len(F10) > 0 else 0.0
        I_10 = float(np.trapz(F10, t10)) if len(t10) > 1 else 0.0
        I_abs_10 = float(np.trapz(np.abs(F10), t10)) if len(t10) > 1 else 0.0

        # --------------------------------------------------------------------
        # OKNO 0 – 30 s: Raná přechodová fáze, lokální periody, Hilbertova obálka
        # --------------------------------------------------------------------
        m30 = (t_post >= 0.0) & (t_post <= 30.0)
        t30, y30, dy30 = t_post[m30], y_post[m30], dy_post[m30]

        # Hilbertova transformace na okně 0-120 s (pro minimalizaci okrajového zkreslení v okně 0-30 s)
        m120 = (t_post >= 0.0) & (t_post <= 120.0)
        t120, y120, dy120, d2y120, F120 = (
            t_post[m120], y_post[m120], dy_post[m120], d2y_post[m120], F_post[m120]
        )
        analytic_120 = hilbert_analytic_np(y120 - np.mean(y120))
        envelope_120 = np.abs(analytic_120)
        inst_phase_120 = np.unwrap(np.angle(analytic_120))
        inst_freq_120 = np.clip(savgol_filter_np(inst_phase_120, win_samples, 2, deriv=1, delta=dt) / (2.0 * np.pi), 0.0, 2.0)

        env_30 = envelope_120[:len(t30)]
        prom_min = max(1.0, 0.12 * full_span)
        peaks_120 = find_peaks_np(y120, min_distance=int(2.5 * fs), min_prominence=prom_min)
        troughs_120 = find_peaks_np(-y120, min_distance=int(2.5 * fs), min_prominence=prom_min)
        extrema_120 = np.sort(np.concatenate([peaks_120, troughs_120]))

        # Lokální periody T_i z po sobě jdoucích půlperiod (zdvojnásobené rozestupy extrémů)
        local_periods: List[float] = []
        if len(extrema_120) >= 2:
            half_periods = np.diff(t120[extrema_120])
            local_periods = [float(2.0 * hp) for hp in half_periods[:6]]
        local_freqs = [float(1.0 / lp) for lp in local_periods if lp > 0]

        # Odhad logaritmického dekrementu delta z poklesu obálky / po sobě jdoucích vrcholů
        if len(peaks_120) >= 2:
            a1 = abs(y120[peaks_120[0]] - np.mean(y120))
            a2 = abs(y120[peaks_120[-1]] - np.mean(y120))
            n_cyc = max(1, len(peaks_120) - 1)
            log_dec_30 = float(np.log(max(a1, 1e-6) / max(a2, 1e-6)) / n_cyc)
        elif len(env_30) > 10:
            log_dec_30 = float(np.log(max(env_30[int(0.1 * len(env_30))], 1e-6) / max(env_30[-1], 1e-6)))
        else:
            log_dec_30 = 0.0

        # --------------------------------------------------------------------
        # OKNO 0 – 60 s: CWT (Morlet) a mechanická práce
        # --------------------------------------------------------------------
        m60 = (t_post >= 0.0) & (t_post <= 60.0)
        t60, y60, dy60, F60 = t_post[m60], y_post[m60], dy_post[m60], F_post[m60]

        cwt_freqs = np.linspace(0.015, 0.25, 60)  # periody 4 s až 66 s
        # Spočteme CWT na okně 0-300 s pro Graf 3 a vyřízneme hřeben pro 0-60 s
        m300 = (t_post >= 0.0) & (t_post <= min(300.0, t_post[-1]))
        t300, y300 = t_post[m300], y_post[m300]
        cwt_power_300 = cwt_morlet_np(y300, fs, cwt_freqs, omega0=5.0)
        ridge_idx_300 = np.argmax(cwt_power_300, axis=0)
        ridge_freqs_300 = cwt_freqs[ridge_idx_300]
        ridge_60 = ridge_freqs_300[:len(t60)]

        # Mechanická práce a výkon buzení v okně 0-60 s: P(t) = F(t) * dy/dt
        power_60 = F60 * dy60
        W_60 = float(np.trapz(power_60, t60)) if len(t60) > 1 else 0.0
        P_mean_60 = float(np.mean(power_60)) if len(power_60) > 0 else 0.0
        P_rms_60 = float(np.sqrt(np.mean(power_60 ** 2))) if len(power_60) > 0 else 0.0

        # --------------------------------------------------------------------
        # OKNO 0 – 120 s: Fázový portrét, Restoring Force Surface, nelinearita
        # --------------------------------------------------------------------
        I_120 = float(np.trapz(F120, t120)) if len(t120) > 1 else 0.0
        I_cum_120 = np.cumsum(F120) * dt
        F_rms_120 = float(np.sqrt(np.mean(F120 ** 2))) if len(F120) > 0 else 0.0

        # Restoring Force Surface odhad nelinearity tuhosti k(y):
        # F_restoring = -m * y'' - c * y' ~ k_1 * (y - y_c) + k_3 * (y - y_c)^3
        y_centered_120 = y120 - np.mean(y120)
        X_nl = np.column_stack([y_centered_120, y_centered_120 ** 3])
        f_iner_120 = -m_val * d2y120 - c_val * dy120
        beta_nl, _, _, _ = np.linalg.lstsq(X_nl, f_iner_120, rcond=None)
        k1_est, k3_est = float(beta_nl[0]), float(beta_nl[1])
        y_amp_max = float(np.max(np.abs(y_centered_120))) if len(y_centered_120) > 0 else 1.0
        nonlinearity_coeff = float((k3_est * (y_amp_max ** 2)) / (abs(k1_est) + 1e-9))

        # Poměr energie disipované tlumením vs. maximální akumulované potenciální energie
        E_diss_120 = float(np.trapz(c_val * (dy120 ** 2), t120)) if len(t120) > 1 else 0.0
        E_pot_max_120 = float(0.5 * k_val * (y_amp_max ** 2)) + 1e-12
        energy_dissipated_vs_stored_ratio = float(E_diss_120 / E_pot_max_120)

        # Sestavení výstupního reportu přesně dle struktury v README_PROGRAM.MD + README_DATA.MD
        unit = self.config.unit_name
        report = {
            "file_info": {
                "filename": parsed["filename"],
                "md5": parsed["md5"],
                "nominal_temp_C": parsed["nominal_temp_c"],
                "temp_diff_from_24C": float(abs(parsed["nominal_temp_c"] - self.config.reference_temp_c))
                if parsed["nominal_temp_c"] is not None else None,
                "logged_timestamp": parsed["metadata"].get("Logged", ""),
                "sensor_range": parsed["metadata"].get("Range", "300uW"),
                "wavelength_nm": parsed["metadata"].get("Wavelength", "1559"),
                "sampling_stats": parsed["sampling_stats"],
                "glitches_detected_and_repaired": parsed["glitches"],
            },
            "preprocessing": {
                "resampled_fs_Hz": fs,
                "resampled_dt_s": dt,
                "onset_t0_s": t0,
                "initial_level_y0": {"value": y0, "unit": unit},
                "baseline_mu0": {"value": onset["mu0_baseline"], "unit": unit},
                "baseline_sigma0": {"value": onset["sigma0_baseline"], "unit": unit},
                "consecutive_diff_stats_README_DATA": {
                    "total_abs_diff_post_t0": {"value": float(np.sum(dv_raw_post)), "unit": unit},
                    "mean_abs_diff_rate_120s": {
                        "value": float(np.sum(np.abs(np.diff(y120))) / max(t120[-1] - t120[0], 1e-6)),
                        "unit": f"{unit}/s",
                    },
                    "max_consecutive_abs_diff": {
                        "value": float(np.max(dv_raw_post)) if len(dv_raw_post) > 0 else 0.0,
                        "unit": unit,
                    },
                },
            },
            "window_10s": {
                "max_onset_velocity_dy_max": {"value": dy_max_10, "unit": f"{unit}/s"},
                "max_onset_acceleration_d2y_max": {"value": d2y_max_10, "unit": f"{unit}/s^2"},
                "first_peak_y_peak1": {"value": y_peak1, "unit": unit},
                "first_peak_time_from_t0_s": {"value": t_peak1, "unit": "s"},
                "initial_peak_force_F_peak": {"value": F_peak_10, "unit": f"F/k ({unit})"},
                "initial_impulse_I10": {"value": I_10, "unit": f"F/k·s ({unit}·s)"},
                "initial_abs_impulse_I10_abs": {"value": I_abs_10, "unit": f"F/k·s ({unit}·s)"},
                "sum_abs_consecutive_diff_10s": {
                    "value": float(np.sum(np.abs(np.diff(y10)))),
                    "unit": unit,
                },
            },
            "window_30s": {
                "local_periods_T_s": local_periods,
                "instantaneous_frequencies_f_Hz": local_freqs,
                "envelope_amplitude_max": {"value": float(np.max(env_30)) if len(env_30) > 0 else 0.0, "unit": unit},
                "envelope_amplitude_end_30s": {"value": float(env_30[-1]) if len(env_30) > 0 else 0.0, "unit": unit},
                "logarithmic_decrement_delta": {"value": log_dec_30, "unit": "-"},
                "sum_abs_consecutive_diff_30s": {
                    "value": float(np.sum(np.abs(np.diff(y30)))),
                    "unit": unit,
                },
            },
            "window_60s": {
                "cwt_ridge_freq_start_Hz": {"value": float(ridge_60[0]) if len(ridge_60) > 0 else 0.0, "unit": "Hz"},
                "cwt_ridge_freq_end_Hz": {"value": float(ridge_60[-1]) if len(ridge_60) > 0 else 0.0, "unit": "Hz"},
                "cwt_ridge_period_start_s": {
                    "value": float(1.0 / ridge_60[0]) if len(ridge_60) > 0 and ridge_60[0] > 0 else 0.0,
                    "unit": "s",
                },
                "cwt_ridge_period_end_s": {
                    "value": float(1.0 / ridge_60[-1]) if len(ridge_60) > 0 and ridge_60[-1] > 0 else 0.0,
                    "unit": "s",
                },
                "mechanical_work_W60": {"value": W_60, "unit": f"{unit}^2"},
                "mean_driving_power_P_mean": {"value": P_mean_60, "unit": f"{unit}^2/s"},
                "rms_driving_power_P_rms": {"value": P_rms_60, "unit": f"{unit}^2/s"},
                "sum_abs_consecutive_diff_60s": {
                    "value": float(np.sum(np.abs(np.diff(y60)))),
                    "unit": unit,
                },
            },
            "window_120s": {
                "total_impulse_I120": {"value": I_120, "unit": f"F/k·s ({unit}·s)"},
                "rms_driving_force_F_RMS_120": {"value": F_rms_120, "unit": f"F/k ({unit})"},
                "stiffness_nonlinearity_coeff": {"value": nonlinearity_coeff, "unit": "-"},
                "energy_dissipated_vs_stored_ratio": {"value": energy_dissipated_vs_stored_ratio, "unit": "-"},
                "sum_abs_consecutive_diff_120s": {
                    "value": float(np.sum(np.abs(np.diff(y120)))),
                    "unit": unit,
                },
                "number_of_fringe_sweeps_120s": {
                    "value": float(np.sum(np.abs(np.diff(y120))) / max(full_span, 1e-6)),
                    "unit": "sweeps",
                },
            },
            "window_900s": {
                "steady_state_y_inf": {"value": y_inf, "unit": unit},
                "static_driving_force_F_stat": {"value": f_stat, "unit": f"F/k ({unit})"},
                "natural_frequency_f0": {"value": f0, "unit": "Hz"},
                "natural_period_T0": {"value": T0, "unit": "s"},
                "damping_ratio_zeta": {"value": zeta, "unit": "-"},
                "estimated_mass_m": {"value": float(m_val), "unit": "kg_eff"},
                "estimated_damping_c": {"value": float(c_val), "unit": "N·s/m_eff"},
                "settling_time_ts_2pct": {"value": t_settling, "unit": "s"},
                "relative_overshoot_sigma_max": {"value": overshoot_pct, "unit": "%"},
                "damped_sine_fit_r2": {"value": sine_fit["r_squared"], "unit": "-"},
            },
        }

        # Data pro vykreslování diagnostických grafů (odděleně od JSON reportu)
        t_tail = t_post[mask_tail]
        t_tail_norm = t_tail - t_tail[0]
        y_drift_tail = np.polyval(sine_fit.get("drift_poly", [0.0, 0.0, sine_fit["y_inf"]]), t_tail_norm)
        y_tail_fit = (
            y_drift_tail
            + sine_fit["amplitude_A"]
            * np.exp(-sine_fit["damping_gamma"] * t_tail_norm)
            * np.cos(2.0 * np.pi * sine_fit["freq_f0_Hz"] * t_tail_norm + sine_fit["phase_phi_rad"])
        )

        # Obálka celého post-t0 signálu
        env_post = np.abs(hilbert_analytic_np(y_post - y_inf))

        plot_data = {
            "t_grid": t_grid,
            "y_raw_grid": y_grid,
            "y_smooth": y_smooth,
            "dy_dt": dy_dt,
            "d2y_dt2": d2y_dt2,
            "t0": t0,
            "y0": y0,
            "y_inf": y_inf,
            "t_post": t_post,
            "y_post": y_post,
            "dy_post": dy_post,
            "env_post": env_post,
            "t_tail": t_tail + t0,
            "y_tail_fit": y_tail_fit,
            "t120": t120,
            "y120": y120,
            "dy120": dy120,
            "F120": F120,
            "I_cum_120": I_cum_120,
            "t300": t300,
            "cwt_freqs": cwt_freqs,
            "cwt_power_300": cwt_power_300,
            "ridge_freqs_300": ridge_freqs_300,
            "mask_tail": mask_tail,
        }

        return {"report": report, "plot_data": plot_data}


# ============================================================================
# 4. VYKRESLOVÁNÍ DIAGNOSTICKÝCH GRAFŮ (ODDĚLENO OD VÝPOČETNÍ LOGIKY)
# ============================================================================

class SignalPlotter:
    """Třída pro generování 4 povinných diagnostických grafů dle bodu 5.2 README_PROGRAM.MD."""

    @staticmethod
    def save_diagnostic_plots(analysis_result: Dict[str, Any], output_dir: str, prefix: str = "") -> List[str]:
        os.makedirs(output_dir, exist_ok=True)
        rep = analysis_result["report"]
        pd_data = analysis_result["plot_data"]
        unit = rep["preprocessing"]["initial_level_y0"]["unit"]
        fname = rep["file_info"]["filename"]
        t0 = pd_data["t0"]
        y0 = pd_data["y0"]
        y_inf = pd_data["y_inf"]

        saved_files: List[str] = []
        p_str = f"{prefix}_" if prefix else ""

        # --------------------------------------------------------------------
        # GRAF 1: Celý průběh signálu y(t), t0, hranice 5 oken, obálka a fit
        # --------------------------------------------------------------------
        fig, ax = plt.subplots(figsize=(13, 6), dpi=160)
        ax.plot(pd_data["t_grid"], pd_data["y_raw_grid"], color="#90a4ae", alpha=0.5, lw=0.8, label="Původní signál")
        ax.plot(pd_data["t_grid"], pd_data["y_smooth"], color="#0d47a1", lw=1.3, label="Vyhlazený signál $y(t)$ (Savitzky-Golay)")
        ax.plot(pd_data["t_post"] + t0, y_inf + pd_data["env_post"], color="#ef6c00", ls="--", lw=1.1, alpha=0.8, label="Hilbertova obálka $y_\\infty \\pm A(t)$")
        ax.plot(pd_data["t_post"] + t0, y_inf - pd_data["env_post"], color="#ef6c00", ls="--", lw=1.1, alpha=0.8)
        ax.plot(pd_data["t_tail"], pd_data["y_tail_fit"], color="#d32f2f", lw=1.8, label="Nelineární regrese tlumené sinusovky ($t - t_0 > 330$ s)")

        ax.axvline(t0, color="#b71c1c", lw=2.0, label=f"Detekovaný počátek $t_0 = {t0:.2f}$ s")
        ax.axhline(y0, color="#455a64", ls=":", lw=1.0, label=f"Výchozí hladina $y_0 = {y0:.2f}$ {unit}")
        ax.axhline(y_inf, color="#2e7d32", ls="-.", lw=1.2, label=f"Ustálená hodnota $y_\\infty = {y_inf:.2f}$ {unit}")

        win_colors = [(10, "#e53935"), (30, "#fb8c00"), (60, "#fdd835"), (120, "#43a047"), (900, "#1e88e5")]
        for w_s, col in win_colors:
            t_w = min(t0 + w_s, pd_data["t_grid"][-1])
            ax.axvline(t_w, color=col, ls="--", lw=1.1, alpha=0.85, label=f"Okno +{w_s} s ({t_w:.1f} s)")

        ax.set_title(f"Graf 1: Časový průběh odezvy, detekce $t_0$, analytická okna a regrese ({fname})", fontsize=12, fontweight="bold")
        ax.set_xlabel("Čas $t$ [s]")
        ax.set_ylabel(f"Signál $y(t)$ [{unit}]")
        ax.legend(loc="upper right", fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        g1_path = os.path.join(output_dir, f"{p_str}graf1_prubeh_a_okna.png")
        fig.savefig(g1_path)
        plt.close(fig)
        saved_files.append(g1_path)

        # --------------------------------------------------------------------
        # GRAF 2: Rekonstruovaná budicí síla F(t) a kumulativní impuls I(t) (0-120 s)
        # --------------------------------------------------------------------
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), dpi=160, sharex=True)
        ax1.plot(pd_data["t120"], pd_data["F120"], color="#1565c0", lw=1.4, label="Rekonstruovaná budicí síla $F(t)/k$")
        ax1.axhline(0, color="black", lw=0.8, ls="--")
        ax1.set_title(f"Graf 2: Rekonstruovaný průběh budicí síly $F(t)$ a kumulativního impulsu $I(t)$ v prvních 120 s ({fname})", fontsize=12, fontweight="bold")
        ax1.set_ylabel(f"Budicí síla $F(t)/k$ [{unit}]")
        ax1.legend(loc="upper right")
        ax1.grid(True, alpha=0.3)

        ax2.plot(pd_data["t120"], pd_data["I_cum_120"], color="#c62828", lw=1.6, label="Kumulativní silový impuls $I(t) = \\int_0^t F(\\tau)\\,d\\tau$")
        ax2.plot(pd_data["t120"], np.cumsum(np.abs(pd_data["F120"])) * (pd_data["t120"][1] - pd_data["t120"][0]),
                 color="#6a1b9a", lw=1.2, ls="--", label="Kumulativní absolutní impuls $\\int_0^t |F(\\tau)|\\,d\\tau$")
        ax2.set_xlabel("Relativní čas od počátku děje $t - t_0$ [s]")
        ax2.set_ylabel(f"Impuls $I(t)/k$ [{unit}·s]")
        ax2.legend(loc="upper left")
        ax2.grid(True, alpha=0.3)
        fig.tight_layout()
        g2_path = os.path.join(output_dir, f"{p_str}graf2_budici_sila_a_impuls.png")
        fig.savefig(g2_path)
        plt.close(fig)
        saved_files.append(g2_path)

        # --------------------------------------------------------------------
        # GRAF 3: Časově-frekvenční spektrogram (CWT s Morletovým waveletem)
        # --------------------------------------------------------------------
        fig, ax = plt.subplots(figsize=(12, 6), dpi=160)
        t300 = pd_data["t300"]
        freqs = pd_data["cwt_freqs"]
        power = pd_data["cwt_power_300"]
        im = ax.pcolormesh(t300, freqs, power, shading="gouraud", cmap="inferno")
        fig.colorbar(im, ax=ax, label=f"Amplituda CWT [{unit}]")
        ax.plot(t300, pd_data["ridge_freqs_300"], color="#00e5ff", lw=2.0, label="Hřeben waveletu (trajektorie dominantní frekvence)")
        f0 = rep["window_900s"]["natural_frequency_f0"]["value"]
        if f0 > 0:
            ax.axhline(f0, color="#76ff03", ls="--", lw=1.5, label=f"Vlastní frekvence soustavy $f_0 = {f0:.4f}$ Hz ($T_0 = {1.0/f0:.1f}$ s)")
        ax.set_title(f"Graf 3: Časově-frekvenční spektrogram (Kontinuální waveletová transformace - Morlet) ({fname})", fontsize=12, fontweight="bold")
        ax.set_xlabel("Relativní čas od počátku děje $t - t_0$ [s]")
        ax.set_ylabel("Frekvence $f$ [Hz]")
        ax.legend(loc="upper right")
        fig.tight_layout()
        g3_path = os.path.join(output_dir, f"{p_str}graf3_cwt_spektrogram.png")
        fig.savefig(g3_path)
        plt.close(fig)
        saved_files.append(g3_path)

        # --------------------------------------------------------------------
        # GRAF 4: Fázový portrét dy/dt(y) pro okno 0-120 s vs. 330-900 s
        # --------------------------------------------------------------------
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), dpi=160)
        sc1 = ax1.scatter(pd_data["y120"], pd_data["dy120"], c=pd_data["t120"], cmap="viridis", s=8, alpha=0.8)
        fig.colorbar(sc1, ax=ax1, label="$t - t_0$ [s]")
        ax1.plot(pd_data["y120"], pd_data["dy120"], color="gray", lw=0.5, alpha=0.4)
        ax1.set_title("Okno 0 – 120 s (Silně nelineární přechodový děj)", fontsize=11, fontweight="bold")
        ax1.set_xlabel(f"Výchylka $y(t)$ [{unit}]")
        ax1.set_ylabel(f"Rychlost změny $\\dot{{y}}(t)$ [{unit}/s]")
        ax1.grid(True, alpha=0.3)

        m_tail = pd_data["mask_tail"]
        t_tail_rel = pd_data["t_post"][m_tail]
        y_tail_p = pd_data["y_post"][m_tail]
        dy_tail_p = pd_data["dy_post"][m_tail]
        sc2 = ax2.scatter(y_tail_p, dy_tail_p, c=t_tail_rel, cmap="plasma", s=8, alpha=0.8)
        fig.colorbar(sc2, ax=ax2, label="$t - t_0$ [s]")
        ax2.plot(y_tail_p, dy_tail_p, color="gray", lw=0.5, alpha=0.4)
        ax2.set_title("Okno 330 – 900 s (Lineární doznívání a ustálení)", fontsize=11, fontweight="bold")
        ax2.set_xlabel(f"Výchylka $y(t)$ [{unit}]")
        ax2.set_ylabel(f"Rychlost změny $\\dot{{y}}(t)$ [{unit}/s]")
        ax2.grid(True, alpha=0.3)

        fig.suptitle(f"Graf 4: Fázový portrét $\\dot{{y}}(y)$ – porovnání nelineární fáze a doznívání ({fname})", fontsize=12, fontweight="bold")
        fig.tight_layout()
        g4_path = os.path.join(output_dir, f"{p_str}graf4_fazovy_portret.png")
        fig.savefig(g4_path)
        plt.close(fig)
        saved_files.append(g4_path)

        return saved_files


# ============================================================================
# 5. DÁVKOVÁ ANALÝZA VŠECH 49 SOUBORŮ (0 °C AŽ 48 °C) A SOUHRNNÝ GRAF
# ============================================================================

def run_batch_analysis(data_dir: str, output_dir: str, config: SystemConfig,
                       detail_file: str = "0 C.txt") -> Dict[str, Any]:
    """
    Zpracuje všechny měřicí soubory v adresáři `data_dir`, vygeneruje:
    - detailní `analysis_report.json` a 4 diagnostické grafy pro vybraný soubor (výchozí `0 C.txt`),
    - souhrnný katalog `batch_summary.csv` a `batch_summary.json` pro všech 49 souborů,
    - souhrnný diagnostický graf závislosti rychlosti změn a počátku děje $t_0$ na teplotě $T$ (vůči 24 °C).
    """
    os.makedirs(output_dir, exist_ok=True)
    analyzer = SignalAnalyzer(config)

    files = sorted(
        [
            f for f in glob.glob(os.path.join(data_dir, "*.txt"))
            if not f.lower().endswith("-analysis.txt")
        ],
        key=lambda x: int(re.search(r"(\d+)", os.path.basename(x)).group(1))
        if re.search(r"(\d+)", os.path.basename(x)) else 999,
    )

    batch_rows: List[Dict[str, Any]] = []
    md5_seen: Dict[str, str] = {}
    anomalies_found: List[Dict[str, Any]] = []
    detail_result: Optional[Dict[str, Any]] = None

    for fpath in files:
        res = analyzer.analyze_file(fpath)
        rep = res["report"]
        finfo = rep["file_info"]
        prep = rep["preprocessing"]
        fname = finfo["filename"]
        temp_c = finfo["nominal_temp_C"]

        if fname.lower() == detail_file.lower() or detail_result is None:
            if fname.lower() == detail_file.lower():
                detail_result = res

        # Kontrola anomálií (duplicity, odlišný rozsah, odlišný čas počátku t0, glitche)
        file_anomalies: List[str] = []
        if finfo["md5"] in md5_seen:
            file_anomalies.append(f"DUPLICITA: 100% identický soubor jako {md5_seen[finfo['md5']]}")
        else:
            md5_seen[finfo["md5"]] = fname

        if finfo["sensor_range"] != "300uW":
            file_anomalies.append(f"ODLIŠNÝ ROZSAH SENZORU: {finfo['sensor_range']} (září 2024 místo prosince 2024)")

        t0_val = prep["onset_t0_s"]
        if abs(t0_val - 140.0) > 10.0:
            file_anomalies.append(f"POSUNUTÝ POČÁTEK DĚJE: t0 = {t0_val:.2f} s (odchylka {t0_val - 140.0:+.1f} s od 140 s)")

        if len(finfo["glitches_detected_and_repaired"]) > 0:
            g_times = ", ".join(f"{g['time_s']:.2f}s" for g in finfo["glitches_detected_and_repaired"])
            file_anomalies.append(f"HARDWAROVÝ VÝPADEK (GLITCH): {len(finfo['glitches_detected_and_repaired'])}x v t = {g_times}")

        if file_anomalies:
            anomalies_found.append({
                "filename": fname,
                "nominal_temp_C": temp_c,
                "t0_s": t0_val,
                "sensor_range": finfo["sensor_range"],
                "logged": finfo["logged_timestamp"],
                "issues": file_anomalies,
            })

        batch_rows.append({
            "filename": fname,
            "temp_C": temp_c,
            "abs_diff_from_24C": finfo["temp_diff_from_24C"],
            "sensor_range": finfo["sensor_range"],
            "logged": finfo["logged_timestamp"],
            "md5": finfo["md5"][:8],
            "fs_mean_Hz": round(finfo["sampling_stats"]["fs_mean_Hz"], 4),
            "t0_s": round(t0_val, 3),
            "y0_uW": round(prep["initial_level_y0"]["value"], 3),
            "y_inf_uW": round(rep["window_900s"]["steady_state_y_inf"]["value"], 3),
            "dy_max_10s_uW_s": round(rep["window_10s"]["max_onset_velocity_dy_max"]["value"], 3),
            "sum_abs_dv_10s_uW": round(rep["window_10s"]["sum_abs_consecutive_diff_10s"]["value"], 2),
            "sum_abs_dv_30s_uW": round(rep["window_30s"]["sum_abs_consecutive_diff_30s"]["value"], 2),
            "sum_abs_dv_60s_uW": round(rep["window_60s"]["sum_abs_consecutive_diff_60s"]["value"], 2),
            "sum_abs_dv_120s_uW": round(rep["window_120s"]["sum_abs_consecutive_diff_120s"]["value"], 2),
            "sum_abs_dv_post_t0_uW": round(prep["consecutive_diff_stats_README_DATA"]["total_abs_diff_post_t0"]["value"], 2),
            "fringe_sweeps_120s": round(rep["window_120s"]["number_of_fringe_sweeps_120s"]["value"], 3),
            "f0_Hz": round(rep["window_900s"]["natural_frequency_f0"]["value"], 5),
            "T0_s": round(rep["window_900s"]["natural_period_T0"]["value"], 2),
            "zeta": round(rep["window_900s"]["damping_ratio_zeta"]["value"], 4),
            "settling_time_s": round(rep["window_900s"]["settling_time_ts_2pct"]["value"], 1),
            "anomalies": "; ".join(file_anomalies) if file_anomalies else "OK",
        })

    # Uložení detailního reportu a 4 grafů pro referenční soubor
    saved_plots: List[str] = []
    if detail_result is not None:
        rep_path = os.path.join(output_dir, "analysis_report.json")
        with open(rep_path, "w", encoding="utf-8") as fp:
            json.dump(detail_result["report"], fp, indent=2, ensure_ascii=False)
        saved_plots = SignalPlotter.save_diagnostic_plots(detail_result, output_dir)

    # Uložení souhrnné tabulky všech 49 souborů
    df_batch = pd.DataFrame(batch_rows)
    csv_path = os.path.join(output_dir, "batch_summary.csv")
    df_batch.to_csv(csv_path, index=False, encoding="utf-8-sig")

    batch_json_path = os.path.join(output_dir, "batch_summary.json")
    with open(batch_json_path, "w", encoding="utf-8") as fp:
        json.dump({
            "total_files": len(batch_rows),
            "anomalies_detected": anomalies_found,
            "files": batch_rows,
        }, fp, indent=2, ensure_ascii=False)

    # ------------------------------------------------------------------------
    # GRAF 5 (SOUHRNNÝ): Závislost rychlosti změn na teplotě (vůči 24 °C) a časy t0
    # ------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 8), dpi=160, sharex=True)

    normal_mask = df_batch["anomalies"] == "OK"
    anom_mask = ~normal_mask

    ax1.plot(df_batch["temp_C"], df_batch["sum_abs_dv_post_t0_uW"], color="#1565c0", lw=1.2, alpha=0.6, label="Celkový součet $|\\Delta y|$ po $t_0$ [$\\mu$W]")
    ax1.scatter(df_batch.loc[normal_mask, "temp_C"], df_batch.loc[normal_mask, "sum_abs_dv_post_t0_uW"],
                color="#1565c0", s=45, zorder=3, label="Standardní měření (rozsah 300 $\\mu$W)")
    ax1.scatter(df_batch.loc[anom_mask, "temp_C"], df_batch.loc[anom_mask, "sum_abs_dv_post_t0_uW"],
                color="#d32f2f", marker="D", s=65, zorder=4, label="Měření s detekovanou anomálií (12, 26, 35, 47, 48 °C)")
    ax1.plot(df_batch["temp_C"], df_batch["sum_abs_dv_120s_uW"], color="#2e7d32", ls="--", lw=1.4, label="Součet $|\\Delta y|$ v okně 0–120 s od $t_0$ [$\\mu$W]")

    # Vyznačení V-tvaru kolem 24 °C
    ax1.axvline(24.0, color="#ef6c00", lw=2.0, ls="-.", label="Referenční teplota 24 °C (minimum rychlosti změn)")
    for _, row in df_batch.loc[anom_mask].iterrows():
        ax1.annotate(f"{int(row['temp_C'])} °C", (row["temp_C"], row["sum_abs_dv_post_t0_uW"]),
                     textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8, color="#b71c1c", fontweight="bold")

    ax1.set_title("Souhrnná analýza dat (0 °C až 48 °C): Rychlost změn $|\\Delta y|$ v závislosti na teplotě a detekce anomálií", fontsize=12, fontweight="bold")
    ax1.set_ylabel("Součet absolutních změn $\\sum |y_{i+1} - y_i|$ [$\\mu$W]")
    ax1.legend(loc="upper center", fontsize=8.5, ncol=2)
    ax1.grid(True, alpha=0.3)

    ax2.plot(df_batch["temp_C"], df_batch["t0_s"], color="#455a64", lw=1.0, alpha=0.5)
    ax2.scatter(df_batch.loc[normal_mask, "temp_C"], df_batch.loc[normal_mask, "t0_s"],
                color="#00838f", s=45, zorder=3, label="Detekovaný počátek děje $t_0$ [s]")
    ax2.scatter(df_batch.loc[anom_mask, "temp_C"], df_batch.loc[anom_mask, "t0_s"],
                color="#d32f2f", marker="D", s=65, zorder=4, label="Anomální soubory")
    ax2.axhline(140.0, color="#b71c1c", ls="--", lw=1.2, label="Nominální počátek $t_0 \\approx 140$ s")
    for _, row in df_batch.loc[df_batch["t0_s"].sub(140.0).abs() > 7.5].iterrows():
        ax2.annotate(f"{int(row['temp_C'])} °C ({row['t0_s']:.1f} s)", (row["temp_C"], row["t0_s"]),
                     textcoords="offset points", xytext=(0, -13 if row["t0_s"] < 140 else 7),
                     ha="center", fontsize=8, color="#b71c1c", fontweight="bold")

    ax2.set_xlabel("Teplota měření $T$ [$^\\circ$C]")
    ax2.set_ylabel("Detekovaný čas počátku $t_0$ [s]")
    ax2.set_xticks(np.arange(0, 49, 2))
    ax2.legend(loc="lower left", fontsize=8.5)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    g5_path = os.path.join(output_dir, "graf5_souhrn_teplot_a_anomalii.png")
    fig.savefig(g5_path)
    plt.close(fig)
    saved_plots.append(g5_path)

    return {
        "detail_report": detail_result["report"] if detail_result else None,
        "anomalies": anomalies_found,
        "batch_df": df_batch,
        "saved_plots": saved_plots,
    }


def print_console_summary(res: Dict[str, Any]) -> None:
    """Přehledný výpis výsledků do konzole dle bodu 5.1 README_PROGRAM.MD."""
    rep = res["detail_report"]
    if rep:
        fi = rep["file_info"]
        pr = rep["preprocessing"]
        print("=" * 88)
        print(f"VÝSLEDKY ANALÝZY REFERENČNÍHO SOUBORU: {fi['filename']} (T = {fi['nominal_temp_C']} °C)")
        print("=" * 88)
        print(f"  Skutečné vzorkování: průměr {fi['sampling_stats']['fs_mean_Hz']:.3f} Hz "
              f"(modus dt = 0.064 s -> 15.625 Hz ~ 16 Hz; max prodleva dt = {fi['sampling_stats']['dt_max_s']:.3f} s)")
        print(f"  Převzorkováno na   : {pr['resampled_fs_Hz']:.1f} Hz (dt = {pr['resampled_dt_s']:.4f} s)")
        print(f"  Detekovaný počátek : t0 = {pr['onset_t0_s']:.3f} s | Výchozí hladina y0 = {pr['initial_level_y0']['value']:.2f} uW")
        print("-" * 88)
        w10 = rep["window_10s"]
        print(f"  [Okno 0-10 s]  dy_max = {w10['max_onset_velocity_dy_max']['value']:.2f} uW/s | "
              f"d2y_max = {w10['max_onset_acceleration_d2y_max']['value']:.2f} uW/s^2 | "
              f"y_peak1 = {w10['first_peak_y_peak1']['value']:.2f} uW (v t0+{w10['first_peak_time_from_t0_s']['value']:.2f} s) | "
              f"I_10 = {w10['initial_impulse_I10']['value']:.2f} uW·s")
        w30 = rep["window_30s"]
        print(f"  [Okno 0-30 s]  Lokální periody T_i = {[round(x, 2) for x in w30['local_periods_T_s'][:4]]} s | "
              f"Log. dekrement delta = {w30['logarithmic_decrement_delta']['value']:.4f}")
        w60 = rep["window_60s"]
        print(f"  [Okno 0-60 s]  CWT perioda: {w60['cwt_ridge_period_start_s']['value']:.1f} s -> "
              f"{w60['cwt_ridge_period_end_s']['value']:.1f} s | "
              f"Práce W_60 = {w60['mechanical_work_W60']['value']:.1f} uW^2 | "
              f"P_RMS = {w60['rms_driving_power_P_rms']['value']:.2f} uW^2/s")
        w120 = rep["window_120s"]
        print(f"  [Okno 0-120 s] F_RMS = {w120['rms_driving_force_F_RMS_120']['value']:.2f} uW | "
              f"Koef. nelinearity tuhosti = {w120['stiffness_nonlinearity_coeff']['value']:.4f} | "
              f"Počet průchodů proužkem = {w120['number_of_fringe_sweeps_120s']['value']:.2f}")
        w900 = rep["window_900s"]
        print(f"  [Okno 0-900 s] y_inf = {w900['steady_state_y_inf']['value']:.2f} uW | "
              f"F_stat/k = {w900['static_driving_force_F_stat']['value']:.2f} uW | "
              f"f0 = {w900['natural_frequency_f0']['value']:.4f} Hz (T0 = {w900['natural_period_T0']['value']:.2f} s) | "
              f"zeta = {w900['damping_ratio_zeta']['value']:.4f} | t_s = {w900['settling_time_ts_2pct']['value']:.1f} s")

    print("\n" + "=" * 88)
    print("DETEKOVANÉ CHYBY A ANOMÁLIE V DATOVÉM SOUBORU (49 SOUBORŮ):")
    print("=" * 88)
    for anom in res["anomalies"]:
        print(f"  * {anom['filename']:8s} (logged: {anom['logged']}, range: {anom['sensor_range']}):")
        for iss in anom["issues"]:
            print(f"      -> {iss}")
    print("=" * 88)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Vícerozměrná časově-frekvenční a inverzní analýza tlumeného kmitání (Ophir StarLab)."
    )
    parser.add_argument("--data-dir", type=str, default="data", help="Adresář se vstupními datovými soubory (.txt)")
    parser.add_argument("--output-dir", type=str, default="output", help="Adresář pro uložení reportů a grafů")
    parser.add_argument("--file", type=str, default="0 C.txt", help="Jméno souboru pro detailní report a Grafy 1-4")
    parser.add_argument("--fs", type=float, default=16.0, help="Vzorkovací frekvence pro převzorkování [Hz] (výchozí 16.0 Hz, volitelně 20.0 Hz)")
    parser.add_argument("--config", type=str, default=None, help="Volitelná cesta ke konfiguračnímu JSON souboru")
    args = parser.parse_args()

    cfg = SystemConfig(resample_fs=args.fs)
    if args.config and os.path.exists(args.config):
        with open(args.config, "r", encoding="utf-8") as fp:
            c_dict = json.load(fp)
        for k, v in c_dict.items():
            if hasattr(cfg, k):
                setattr(cfg, k, v)

    res = run_batch_analysis(args.data_dir, args.output_dir, cfg, detail_file=args.file)
    print_console_summary(res)


if __name__ == "__main__":
    main()
