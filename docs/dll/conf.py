from pathlib import Path

project = 'REIMS DLL Developer Guide'
copyright = '2026, ITER Organization'
author = 'Authors: D. Furfaro, J. Kosek'
release = 'V2.1'

extensions = ['sphinx.ext.autodoc']

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
html_theme = 'alabaster'
_static_dir = Path(__file__).parent / '_static'
html_static_path = ['_static'] if _static_dir.is_dir() else []
html_css_files = ['custom.css'] if (_static_dir / 'custom.css').is_file() else []

# -- Options for LaTeX/PDF output --------------------------------------------

latex_elements = {
    'papersize': 'a4paper',
    'sphinxsetup': 'hmargin={2.5cm,2.5cm}, vmargin={2.5cm,2.5cm}',
    'preamble': r'\def\_{\textunderscore\penalty0}',
    'extraclassoptions': 'openany',
}
