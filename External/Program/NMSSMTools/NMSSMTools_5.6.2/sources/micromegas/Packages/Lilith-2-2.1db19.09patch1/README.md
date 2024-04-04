# Lilith-2

Lilith (Light Likelihood Fit for the Higgs) is a light and easy-to-use Python tool for constraining new physics from signal strength measurements of the 125 GeV Higgs boson. Lilith is provided with the latest experimental measurements from the ATLAS and CMS collaborations at the LHC. The Higgs likelihood is based on experimental results stored in an easily extensible XML database; it is evaluated from the user input, given in XML format in terms of reduced couplings or signal strengths. 

New from version 2.0 onwards (more in changelog.txt):

- Use of variable Gaussian and generalized Poisson likelihoods for a __better treatment of assymetric uncertainties__. 
The generalized Poisson likelihood can be used for experimental results in 1 or 2 dimensions (the latter with correlation), while the variable Gaussian approximation is available for results of any dimension. 

- Use of __N-dim correlation matrices__ for ordinary and variable Gaussian likelihoods.

- Database 19.09 contains the __published ATLAS and CMS Run 2 results for 36 fb-1__ as of Sep 2019.


### Download 

The latest official version can be found on the [releases page](https://github.com/sabinekraml/Lilith-2/releases). Note: 

- Version 2.0 is in Python 2 (2.7.4 or higher)
- Version 2.1 is in Python 3 (3.6 or higher)

Other prerequisites are SciPy and NumPy; 
the example codes doing a likelihood profile analysis require iminuit.


### Usage

- The usage is explained in the original Lilith manual [arXiv:1502.04138](https://arxiv.org/abs/1502.04138) and the presentation of Lilith-2 in [arXiv:1908.03952](https://arxiv.org/abs/1908.03952). For usage in micrOMEGAs, see [arXiv:1606.03834](https://arxiv.org/abs/1606.03834).

- A __tutorial__ is available from the [Tools2020 workshop](https://indico.cern.ch/event/955391/contributions/4086275/).

**IMPORTANT: please cite the above papers** when you use Lilith for your work. Thank you.


### Troubleshooting

- In case of problems running the code, check whether the `__init.py__` file exists in lilith/internal/ and is executable. If not, create it (as an empty file) and declare it as executable. If the code still does not work, check that all the Lilith Python (`.py`) files are executable.  

- If the example codes don't run, check whether they are executable; if they are not, do `chmod u+x` ...

- If you get an error `ImportError: No module named lilith`, your path to Lilith is probably not correct. (e.g, when this happens with the Python example codes, check where `lilith_dir` points to)


### Ongoing developments

The current developments towards inclusion and validation of the ATLAS and CMS results for full Run 2 luminosity (~140/fb) are done on the __py3-fullRun2__ branch.

