# Zadání: Program pro vícerozměrnou časově-frekvenční a inverzní analýzu tlumeného kmitání

## 1. Cíl programu

Vytvořit modulární skript v jazyce Python, který načte časový záznam odezvy fyzikálního systému, 
automaticky detekuje okamžik vzniku poruchy ($t_0$) 
a provede specializovanou analýzu v pěti vnořených časových oknech po vzniku změny: **10 s, 30 s, 60 s, 120 s a 900 s**. 
Cílem je kvantifikovat dynamiku odezvy a zrekonstruovat velikost a průběh poruchové veličiny $F(t)$.

## 2. Vstupní data a konfigurace
* Vstupní soubor: Textový soubor (.txt) je log optického měřicího systému Ophir StarLab Version

## 3. Povinné předzpracování signálu (Pre-processing)
* Určení výchozí hladiny ($y_0$) 
* detekce počátku ($t_0$) jevu 
* většinou kolem 140 vteřiny 
* délka pokusu 900 vteřin

* Výpočet klouzavého průměru $\mu_0$ a směrodatné odchylky $\sigma_0$ z klidové fáze před změnou.
* Automatická detekce okamžiku $t_0$ (v dodaném vzorku $t_0 \approx 141{,}12\text{ s}$), 
  kdy signál překročí práh $\mu_0 \pm 5\sigma_0$ a trend trvá alespoň $0{,}5\text{ s}$.
* Vyhlazení pro numerické derivace: Aplikace Savitzkyho–Golayova filtru (`scipy.signal.savgol_filter`) 
   pro výpočet rychlosti $\dot{y}(t)$ a zrychlení $\ddot{y}(t)$ bez zesílení vysokofrekvenčního šumu.

* vykreslení těchto průběhů

### 5. Požadované výstupy programu

1. **Strukturovaný report (`analysis_report.json` a přehledný výpis do konzole):**
* Slovník rozdělený podle klíčů `window_10s`, `window_30s`, `window_60s`, `window_120s`, `window_900s` obsahující všechny vypočtené metriky včetně fyzikálních jednotek.


2. **Sada diagnostických grafů (uložených jako `.png` ve vysokém rozlišení):**
* *Graf 1:* Celý průběh signálu $y(t)$ s vyznačeným okamžikem $t_0$, hranicemi všech 5 oken (10, 30, 60, 120, 900 s), obálkou a proloženou tlumenou sinusovkou v koncové části.
* *Graf 2:* Rekonstruovaný časový průběh budicí síly $F(t)$ a jejího kumulativního impulsu $I(t)$ v prvních 120 s.
* *Graf 3:* Časově-frekvenční spektrogram (CWT) ukazující pokles frekvence v první fázi a její ustálení na vlastní frekvenci $f_0$.
* *Graf 4:* Fázový portrét $\dot{y}(y)$ pro okno 0–120 s vs. 330–900 s (vizualizace nelinearity a útlumu).


3. **Kvalita kódu:**
* Čistý Python 3.10+, využití standardních vědeckých knihoven (`numpy`, `scipy`, `pandas`, `matplotlib`, volitelně `PyWavelets`).
* Oddělení výpočetní logiky (v samostatných modulech) 
  * třída `MathAnalyzer` - metody analýzy signálu
  * třída `SignalPlotter` - vykresluje grafy
  * třída `SignalAnalyzer` - řízení volání matematických metod (volba, úprava dat)
  * třída `MainProgramm`  - práce nad jedním nebo více souboiry, ukládání výsledků, řízení workflow