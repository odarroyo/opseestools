import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = 'opseestools'
author = 'opseestools contributors'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
]

autodoc_mock_imports = ['openseespy', 'matplotlib', 'numpy', 'scipy', 'pandas']

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
