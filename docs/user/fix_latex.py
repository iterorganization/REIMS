import re

with open('_build/latex/reims.tex', 'r', encoding='utf-8') as f:
    content = f.read()

# Fix longtable {lll}/{llll}: colonnes 'l' non-wrappantes -> debordement a droite
# \X{1}{N} = largeur fixe = sphinxcolwidth{1}{N} -> coherence avec varwidth
count_lll = content.count(r'\begin{longtable}{lll}')
content = content.replace(
    r'\begin{longtable}{lll}',
    r'\begin{longtable}{\X{1}{3}\X{1}{3}\X{1}{3}}'
)

count_llll = content.count(r'\begin{longtable}{llll}')
content = content.replace(
    r'\begin{longtable}{llll}',
    r'\begin{longtable}{\X{1}{4}\X{1}{4}\X{1}{4}\X{1}{4}}'
)

# Fix tabulary TTT/TTTT: multirow dans col1 -> largeur ~0 -> noms de prop. tronques
count_ttt = content.count(r'\begin{tabulary}{\linewidth}[t]{TTT}')
content = content.replace(
    r'\begin{tabulary}{\linewidth}[t]{TTT}',
    r'\begin{tabulary}{\linewidth}[t]{p{4.5cm}LL}'
)

count_tttt = content.count(r'\begin{tabulary}{\linewidth}[t]{TTTT}')
content = content.replace(
    r'\begin{tabulary}{\linewidth}[t]{TTTT}',
    r'\begin{tabulary}{\linewidth}[t]{p{4.0cm}LLL}'
)

# Fix tabular X/X/X and X/X/X/X schema tables: equal-width columns squeeze
# property names and long references in PDF. Use asymmetric wrapped columns.
count_tabular3 = content.count(r'\begin{tabular}[t]{*{3}{\X{1}{3}}}')
content = content.replace(
    r'\begin{tabular}[t]{*{3}{\X{1}{3}}}',
    r'\begin{tabular}[t]{p{4.2cm}p{2.6cm}p{6.7cm}}'
)

count_tabular4 = content.count(r'\begin{tabular}[t]{*{4}{\X{1}{4}}}')
content = content.replace(
    r'\begin{tabular}[t]{*{4}{\X{1}{4}}}',
    r'\begin{tabular}[t]{p{3.9cm}p{2.0cm}p{2.0cm}p{5.1cm}}'
)

with open('_build/latex/reims.tex', 'w', encoding='utf-8') as f:
    f.write(content)

print(
    'Fixed '
    f'{count_lll} lll, {count_llll} llll longtables; '
    f'{count_ttt} TTT, {count_tttt} TTTT tabulary; '
    f'{count_tabular3} 3-col and {count_tabular4} 4-col tabular tables'
)
