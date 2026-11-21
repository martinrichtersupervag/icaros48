# `CorrectTimeStampOphirVega.py` – Oprava časové osy Ophir StarLab / Vega PD300-IR

## Účel skriptu

Skript `CorrectTimeStampOphirVega.py` opravuje kvantovanou časovou osu (`Timestamp`) v log souborech z měřiče Ophir StarLab / Vega PD300-IR, přičemž všechny naměřené hodnoty (`Channel A`) i strukturu souborů ponechává beze změny.

## Parametry a chování

- **Vstupní adresář:** `data/` (výchozí; originální soubory zůstávají nedotčené)
- **Výstupní adresář:** `data_correctedTimeStamp/` (výchozí)
- **Vzorkovací frekvence:** `15.0 Hz` ($\Delta t = 1/15\text{ s} = 0{,}066667\text{ s}$)
- **Výpočet časové značky:** `t[i] = i / 15.0`, formátováno na **6 desetinných míst** se zachováním šířky sloupce (14 znaků) a konců řádků (`CRLF`)
- **Sloupec `Channel A`:** zachován **bajt po bajtu beze změny** (žádná interpolace ani spline)
- **Hlavička souboru:** zachována **bajt po bajtu beze změny**, pouze do řádku `;Notes:` je doplněna poznámka:
  ```text
  ;Notes: 21112026a Corrected timestamps, sample frequency 15 Hz
  ```

### Spuštění

```bash
python CorrectTimeStampOphirVega.py
```

Volitelné přepínače:

```bash
python CorrectTimeStampOphirVega.py --input-dir data --output-dir data_correctedTimeStamp --fs 15.0
```

## Výsledek zpracování a ověření (21. 11. 2026)

- **Zpracováno:** 49 souborů (`0 C.txt` – `48 C.txt`), celkem **661 676 vzorků**
- **Kontrola integrity:** automaticky ověřeno nad všemi 49 soubory:
  - shodný počet souborů i řádků s `data/`,
  - všechny řádky hlavičky kromě `;Notes:` jsou bajt po bajtu identické,
  - řádek `;Notes:` obsahuje přesně `21112026a Corrected timestamps, sample frequency 15 Hz`,
  - celý sloupec `Channel A` (včetně tabulátorů a konců řádků) je bajt po bajtu identický na všech 661 676 řádcích,
  - všechny hodnoty `Timestamp` přesně odpovídají `i / 15.0` na 6 desetinných míst (`0.000000`, `0.066667`, `0.133333`, …).
- **Poznámka ke kontrole reziduí:** vestavěná pojistka upozornila u 7 souborů (`0 C.txt`, `4 C.txt`, `16 C.txt`, `19 C.txt`, `20 C.txt`, `30 C.txt`, `32 C.txt`) na fázový posun původní pilovité časové osy 40–64 ms; následná kontrola potvrdila, že v žádném souboru nechybí žádný vzorek (všechna původní $\Delta t$ leží výhradně v 63–65 ms nebo 95–97 ms).

### Ukázka výstupu (`data_correctedTimeStamp/1 C.txt`)

```text
;Notes: 21112026a Corrected timestamps, sample frequency 15 Hz
...
;First Pulse Arrived : 13/12/2024 at 10:20:12.915000
  Timestamp   	  Channel A  	
      0.000000	  2.400e-006 	
      0.066667	  2.400e-006 	
      0.133333	  2.400e-006 	
...
    900.000000	  1.440e-005 	
    900.066667	  1.450e-005 	
    900.133333	  1.450e-005 	
```
