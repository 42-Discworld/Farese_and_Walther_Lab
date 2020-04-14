## **Lipidomics Project**

### Introduction
This analysis pipeline for mass spectrometry-based lipidomics (DDA) studies is being continuously developed and improved in the laboratory of Drs. Robert Farese, Jr. and Tobias Walther at the Harvard T.H. Chan School of Public Health.

### Installation
Requires R version 3.6.0 or newer. Any additional packages/libraries will be installed and loaded automatically. At present, pipeline supports mass-spectrometry derived data generated from LipidSearch (Thermo Fischer Scientific/Mitsui Knowledge Industry) or similar software for post-processing such quality control functions, statistical quantification and pathway visualization on experimental data.

Suggesting using homebrew for installing XQaurtz and gfortran. Open your terminal and copy the code below.

For Homebrew:   
`bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"`

For XQaurtz: 

- `ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" < /dev/null 2> /dev/null ; brew install caskroom/cask/brew-cask 2> /dev/null`

- `brew cask install xquartz`

For gfortran:
`brew install gcc`

### Caveats
For Mac users, please make sure XQaurtz and gfortran are installed for prevent crashing.

### Feedback
wlyu@hsph.harvard.edu

lyuwting@gmail.com

### Contributors
Wenting Lyu, Manuele Piccolis, Niklas Mejhert, Kun Wang, Zon Weng Lai, Laura Bond & Sebastian Boland

Farese & Walther Lab

Department of Genetics and Complex Diseases, Harvard T.H. Chan School of Public Health

Department of Cell Biology, Harvard Medical School

### Acknowledgements

Carson, Nicole, Whitney
