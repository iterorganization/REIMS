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
    r'\begin{tabulary}{\linewidth}[t]{p{3.2cm}p{2cm}p{9cm}}'
)

count_tttt = content.count(r'\begin{tabulary}{\linewidth}[t]{TTTT}')
content = content.replace(
    r'\begin{tabulary}{\linewidth}[t]{TTTT}',
    r'\setlength{\tymax}{5cm} \begin{tabulary}{\linewidth}[t]{p{3.2cm}p{2cm}p{3cm}L}'
)

count_ttttt = content.count(r'\begin{tabulary}{\linewidth}[t]{TTTTT}')
content = content.replace(
    r'\begin{tabulary}{\linewidth}[t]{TTTTT}',
    r'\begin{tabulary}{\linewidth}[t]{p{3.2cm}p{2cm}LLL}'
)

# Reduce the horizontal indentation of itemize lists, including those in tables.
list_preamble = r'''\usepackage{enumitem}
\setlist[itemize]{leftmargin=1.1em,labelsep=0.25em}
'''
count_document = content.count(r'\begin{document}')
content = content.replace(
    r'\begin{document}',
    list_preamble + r'\begin{document}'
)

# Fix spacing issue in 'time_between_2D_writes'
content = content.replace(
    r'time\_between\_2D\_writes',
    r'time\_between\_ 2D\_writes'
)

# inject new chapters for better organization
content = content.replace(
    r'\subsubsection{channel}',
    r'\subsection{State components}\subsubsection{channel}'
)
content = content.replace(
    r'\subsubsection{junction}',
    r'\subsection{Link components}\subsubsection{junction}'
)
content = content.replace(
    r'\subsubsection{0D\_signal}',
    r'\subsection{Common definitions}\subsubsection{0D\_signal}'
)


with open('_build/latex/reims.tex', 'w', encoding='utf-8') as f:
    f.write(content)

print(
    'Fixed '
    f'{count_lll} lll, {count_llll} llll longtables; '
    f'updated itemize indentation in {count_document} document preamble; '
    f'{count_ttt} TTT, {count_tttt} TTTT tabulary; '
    f'{count_ttttt} TTTTT tabulary'
)
