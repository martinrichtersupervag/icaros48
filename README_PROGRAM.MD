
## Zadání: Program pro vícerozměrnou časově-frekvenční a inverzní analýzu tlumeného kmitání

### 1. Cíl programu

Vytvořit modulární skript v jazyce Python, který načte časový záznam odezvy mechanického/fyzikálního systému, automaticky detekuje okamžik vzniku poruchy ($t_0$) a provede specializovanou analýzu v pěti vnořených časových oknech po vzniku změny: **10 s, 30 s, 60 s, 120 s a 900 s**. Cílem je kvantifikovat dynamiku odezvy a zrekonstruovat velikost a průběh budicí síly $F(t)$.

---

### 2. Vstupní data a konfigurace

* **Vstupní soubor:** Textový soubor (`.txt` / `.csv` / `.tsv`) se dvěma sloupci bez hlavičky (nebo s volitelnou hlavičkou):
1. `time` – čas v sekundách $[\text{s}]$ (pozor: vzorkovací krok kolísá mezi $\approx 0{,}063\text{ s}$ a $0{,}096\text{ s}$).
2. `value` – měřená výchylka / deformace / napětí v základní jednotce.


* **Konfigurační parametry systému (zadávané uživatelem v hlavičce skriptu nebo přes `config.json`):**
* Statická tuhost soustavy $k$, efektivní hmotnost $m$, koeficient tlumení $c$ (pokud nejsou $m, c$ známy, program je odhadne z okna 900 s a sílu vyjádří v relativních jednotkách $F/k$).
* Převodní konstanta senzoru (přepočet jednotek signálu na fyzikální posuv $[\text{m}]$).



---

### 3. Povinné předzpracování signálu (Pre-processing)

1. **Převzorkování (Resampling):** Převod neekvidistantního časového záznamu na přesně konstantní vzorkovací frekvenci $f_s = 20\text{ Hz}$ ($\Delta t = 0{,}05\text{ s}$) pomocí kubické splajnové interpolace (`scipy.interpolate.CubicSpline`).
2. **Určení výchozí hladiny ($y_0$) a detekce počátku ($t_0$):**
* Výpočet klouzavého průměru $\mu_0$ a směrodatné odchylky $\sigma_0$ z klidové fáze před změnou.
* Automatická detekce okamžiku $t_0$ (v dodaném vzorku $t_0 \approx 141{,}12\text{ s}$), kdy signál překročí práh $\mu_0 \pm 5\sigma_0$ a trend trvá alespoň $0{,}5\text{ s}$.


3. **Vyhlazení pro numerické derivace:** Aplikace Savitzkyho–Golayova filtru (`scipy.signal.savgol_filter`) pro výpočet rychlosti $\dot{y}(t)$ a zrychlení $\ddot{y}(t)$ bez zesílení vysokofrekvenčního šumu.

---

### 4. Specifikace analýzy pro jednotlivá časová okna ($\Delta t$ od $t_0$)

Každé časové okno zachycuje jinou fyzikální fázi děje, proto vyžaduje odlišný matematický aparát:

| Okno od $t_0$ | Fyzikální charakteristika fáze (dle vzorku) | Doporučené matematické metody | Požadované výstupní veličiny |
| --- | --- | --- | --- |
| **0 – 10 s** | Okamžitý náběh a první kmit ($T \approx 9{,}8\text{ s}$). Spektrální metody zde selhávají (málo period). | • Numerická derivace (Savitzky-Golay)<br>

<br>• Analýza náběžné hrany (Step/Impulse onset)<br>

<br>• Inverzní rovnice pohybu $F(t) = m\ddot{y} + c\dot{y} + ky$ | • Maximální rychlost náběhu $\dot{y}_{\max}$ a zrychlení $\ddot{y}_{\max}$<br>

<br>• První maximum $y_{\text{peak1}}$ a čas dosažení $t_{\text{peak1}}$<br>

<br>• Odhad počáteční špičkové síly $F_{\text{peak}}$ a počátečního impulsu $I_{10} = \int_0^{10} F(t)\,\mathrm{d}t$ |
| **0 – 30 s** | Raná přechodová fáze (2–3 cykly), rychlé prodlužování periody ($9{,}8\text{ s} \to 18\text{ s}$). | • Detekce lokálních extrémů (`find_peaks`) a průchodů nulou<br>

<br>• Hilbertova transformace (`scipy.signal.hilbert`) s ošetřením okrajů | • Posloupnost lokálních period $T_1, T_2, \dots$ a okamžitých frekvencí $f_i$<br>

<br>• Okamžitá obálka amplitudy $A(t)$<br>

<br>• První odhad logaritmického dekrementu $\delta$ |
| **0 – 60 s** | Rozvinutý nelineární přechodový děj ($T$ roste až na $32{,}4\text{ s}$). | • **Kontinuální waveletová transformace (CWT)** s Morletovým waveletem<br>

<br>• Výpočet mechanické práce | • Časově-frekvenční mapa (hřeben waveletu – trajektorie frekvence v čase)<br>

<br>• Práce vykonaná budicí silou $W_{60} = \int_0^{60} F(t)\dot{y}(t)\,\mathrm{d}t$<br>

<br>• Průměrný a efektivní (RMS) výkon buzení |
| **0 – 120 s** | Kompletní fáze velkých rozkmitů (vynuceného / silně nelineárního kmitání). | • Fázový portrét $(y(t), \dot{y}(t))$<br>

<br>• Metoda plochy vratné síly (*Restoring Force Surface*) pro odhad nelinearity tuhosti $k(y)$ | • Celkový silový impuls $I_{120}$ a efektivní hodnota budicí síly $F_{\text{RMS}, 120}$<br>

<br>• Koeficient nelinearity tuhosti (rozdíl tuhosti při malých a velkých výchylkách)<br>

<br>• Poměr energie disipované tlumením vs. akumulované |
| **0 – 900 s** | Celý děj včetně přechodu do lineárního doznívání ($t > 470\text{ s}$, $T \approx 19{,}5\text{ s}$) a ustálení ($y_\infty$). | • **Nelineární regrese** (`scipy.optimize.curve_fit`) tlumené sinusovky na ustálenou část ($t - t_0 > 330\text{ s}$)<br>

<br>• Analýza skokové odezvy | • **Statická složka budicí síly:** $F_{\text{stat}} = k \cdot (y_\infty - y_0)$<br>

<br>• Vlastní frekvence soustavy $f_0$ a vlastní perioda $T_0$<br>

<br>• Poměrný útlum $\zeta$ a koeficient tlumení $c$<br>

<br>• Doba ustálení $t_s$ (vstup do pásma $\pm 2\text{ \%}$ kolem $y_\infty$) a relativní překmit $\sigma_{\max}$ |

---

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
* Oddělení výpočetní logiky (funkce/třída `SignalAnalyzer`) od vykreslování grafů.