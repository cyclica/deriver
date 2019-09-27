# Deriver

Documentation in progress. See `contrib.examples.basic_usage_example.py` for a basic usage example.

# Requirements
Requires `selfies` and `rdkit` to be installed.

# Installation
rdkit can be installed on Ubuntu 18.04 using the provided script: 

`mkdir ./rdkit`

Then source your virtualenv. Then:

`./install-rdkit ./rdkit/ -p 3.6` 

Otherwise, follow the instructions here: https://github.com/rdkit/rdkit

Then you may install deriver using pip: `pip install deriver`

# To-do:
* implement a random SELIFES generator method
* create a molecule object that can be passed around