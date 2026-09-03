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

with open('_build/latex/reims.tex', 'w', encoding='utf-8') as f:
    f.write(content)

print(f'Fixed {count_lll} lll, {count_llll} llll longtables; {count_ttt} TTT, {count_tttt} TTTT tabulary')
