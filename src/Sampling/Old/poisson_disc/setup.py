# -*- coding: utf-8 -*-
from setuptools import setup

modules = \
['poisson_disc']
install_requires = \
['numpy>=1.20.0,<2.0.0', 'scipy>=1.6.0,<2.0.0']

setup_kwargs = {
    'name': 'poisson-disc',
    'version': '0.2.1',
    'description': 'Poisson disc sampling in arbitrary dimensions using Bridson\'s algorithm, implemented in python using numpy and scipy. Generates so-called "blue noise" that prevents clustering by ensuring each two points are at least "radius" apart.',
    'long_description': '# Poisson disc sampling\nPoisson disc sampling in arbitrary dimensions using Bridson\'s algorithm, implemented in python using ``numpy`` and ``scipy``.\n\nGenerates so-called "blue noise" that prevents clustering by ensuring each two points are at least ``radius`` apart.\n\nhttps://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf\n\nImplementation is located in ``poisson_disc.py``, while ``poisson_disc_sampling.ipynb`` contains some examples. \n\nAvailable through PyPI as ``poisson_disc``, https://pypi.org/project/poisson-disc/',
    'author': 'Pavel Zun',
    'author_email': 'pavel.zun@gmail.com',
    'maintainer': None,
    'maintainer_email': None,
    'url': 'https://github.com/diregoblin/poisson_disc_sampling',
    'py_modules': modules,
    'install_requires': install_requires,
    'python_requires': '>=3.7,<4.0',
}


setup(**setup_kwargs)
