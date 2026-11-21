# Analýza časování vstupního signálu (Ophir StarLab / Vega PD300-IR)

## 1. Diagnóza: Chyba měření času vs. chyba ukládání dat

Jedná se jednoznačně o **chybu měření / kvantování časové značky** (nikoli o chybu ukládání, zpoždění při zápisu nebo výpadek vzorků). Zároveň skutečná vzorkovací frekvence není 16 Hz, ale **15,00 Hz**.

### Ověřená zjištění (nad všemi 49 soubory v `data/*.txt`)

1. **Skutečná průměrná vzorkovací frekvence je ~15,003 Hz ($\Delta t \approx 66{,}65\text{ ms}$):**
   - Všech 49 souborů (`0 C.txt` až `48 C.txt`) má průměrnou vzorkovací frekvenci v rozmezí **15,0028 – 15,0039 Hz** (průměr 15,0033 Hz; 13 503–13 504 vzorků za ~900 s), což odpovídá nativnímu vzorkování měřiče Ophir Vega (15 Hz).
   - V žádném ze 49 souborů se nevyskytuje jediný výpadek vzorku — všechny časové rozdíly $\Delta t$ leží buď kolem **64 ms** (63–65 ms), nebo **96 ms** (95–97 ms).

2. **Příčina skoků „+32 ms každých ~12 vzorků“ (kvantování časové osy):**
   - Zaznamenané časové značky jsou kvantovány na krok **32 ms** ($\Delta t \approx 64\text{ ms} = 2 \times 32\text{ ms}$ a $\Delta t \approx 96\text{ ms} = 3 \times 32\text{ ms}$, s drobným $\pm 1\text{ ms}$ zaokrouhlením).
   - Protože reálný interval vzorkování je $66{,}65\text{ ms}$, každý 64ms krok „ztratí“ ~2,65 ms. Po 12 vzorcích se skluz nasčítá na ~32 ms a časová značka skočí o 96 ms:
     $$11 \times 64\text{ ms} + 1 \times 96\text{ ms} = 800\text{ ms} \implies \frac{800\text{ ms}}{12} = 66{,}67\text{ ms} = \frac{1}{15}\text{ s}$$
   - Rezidua po odečtení lineárního trendu $r_i = t_i - (a \cdot i + b)$ tvoří pravidelnou **pilu o amplitudě $\pm 16\text{ ms}$** (rozpětí špička–špička ~33 ms; např. pro `12 C.txt`: min $-16{,}72\text{ ms}$, max $+17{,}00\text{ ms}$).

3. **Důkaz ze samotného signálu (`Channel A`):**
   - Kdyby senzor nebo ukládání skutečně vyčkalo 96 ms místo 64 ms, byla by změna signálu $\Delta y$ na hladkém sklonu během 96ms kroku **1,5× větší** než v sousedních 64ms krocích.
   - Analýza 574 hladkých úseků napříč všemi 49 soubory ukázala poměr $\Delta y(96\text{ ms}) / \Delta y(64\text{ ms})$:
     - **medián = 1,0000** (průměr 0,96, nikoli 1,50).
   - Fyzikální vzorky `Channel A` jsou tedy nasnímány **rovnoměrně po ~66,65 ms** — zkreslená je pouze časová značka přiřazená při čtení/časování v PC/StarLabu.

---

## 2. Návrh opravy dat

> **Zásada:** Originální soubory v `data/*.txt` **nepřepisovat** (i hlavička StarLab před jejich úpravou varuje). Opravu provádět při načítání dat v paměti (např. v modulu načítání / `SignalAnalyzer`), případně ukládat opravené kopie do samostatného adresáře.

### Doporučený postup: Lineární rekonstrukce časové osy (bez interpolace hodnot `y`)

1. **Proložení časové osy přímkou** přes index vzorku $i = 0, 1, \dots, N-1$ metodou nejmenších čtverců (`np.polyfit(np.arange(N), t, 1)`):
   $$t_{\text{clean}}[i] = a \cdot i + b \quad (\text{případně s posunem od nuly } t_{\text{clean}}[i] = a \cdot i)$$
   kde sklon $a \approx 0{,}066653\text{ s}$ je skutečná perioda vzorkování ($\Delta t_{\text{eff}}$) a efektivní vzorkovací frekvence je $f_s = 1 / a \approx 15{,}003\text{ Hz}$.

2. **Hodnoty signálu `Channel A` ponechat beze změny:**
   - **Neprovádět** interpolaci/převzorkování přímo ze surové (kvantované, pilovité) časové osy. Protože byly hodnoty `y` fyzicky nasnímány rovnoměrně, interpolace přes chybný čas by do signálu uměle zanesla fázovou modulaci s periodou ~0,8 s (frekvence ~1,25 Hz).
   - Použitím `t_clean` se zároveň odstraní umělé 33% propady v numerické derivaci $\Delta y / \Delta t$ na každém 96ms skoku a zajistí se přesně rovnoměrná mřížka nutná pro `savgol_filter`, FFT a CWT spektrogram.

3. **Bezpečnostní validace při načítání (pojistka pro budoucí data):**
   - Spočítat rezidua $r_i = t_i - (a \cdot i + b)$.
   - Pokud $\max |r_i| < 25\text{ ms}$ (platí pro všech 49 současných souborů, kde je max ~17 ms), bezpečně nahradit `t` rovnoměrnou osou `t_clean`.
   - Pokud $\max |r_i| \ge 25\text{ ms}$ (skutečný výpadek vzorků nebo drift), vypsat varování a ošetřit daný soubor / segment zvlášť.

4. **Volitelný krok (pouze při požadavku na jednotnou frekvenci napříč soubory):**
   - Teprve **z opravené osy** `t_clean` převzorkovat signál (např. pomocí `scipy.interpolate.CubicSpline`) na přesně `15.0 Hz` (nebo jinou společnou frekvenci). Pro většinu výpočtů dle `README_PROGRAM.MD` (`savgol_filter`, derivace $\dot{y}, \ddot{y}$, CWT/FFT) ale plně vyhovuje nativní rovnoměrná osa `t_clean` s $\Delta t = a$.
