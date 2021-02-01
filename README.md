# Deriver

Documentation in progress. See `contrib.examples.basic_usage_example.py` for a basic usage example.

## Paper

[Reeves, S., DiFrancesco, B., Shahani, V., MacKinnon, S., Windemuth, A., & Brereton, A. E. (2020). Assessing methods and obstacles in chemical space exploration. Applied AI Letters.](https://onlinelibrary.wiley.com/doi/full/10.1002/ail2.17)

## Requirements
Requires `selfies` and `rdkit` to be installed.

## Installation
rdkit can be installed on Ubuntu 18.04 or 20.04 using the provided script: 

https://github.com/cyclica/rdkit-installer

Source that environment once created, then you may install deriver using pip: `pip install deriver`, or by cloning this 
repo and from the cloned directory running: `pip install -e .`

## To-do:
* refactor the filtering module to be polymorphic and support operators
