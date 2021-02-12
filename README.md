# Deriver

Documentation in progress. See `contrib.examples.basic_usage_example.py` for a basic usage example.

# Requirements
Requires `rdkit` to be installed, along with everything in `requirements.txt`.

To use the `derive_local_space()` method, you will need to make or use a prexisting fragment database. You can find more information here: https://github.com/DrrDom/crem


# Installation
rdkit can be installed on Ubuntu 18.04 using the provided script: 

https://github.com/cyclica/rdkit-installer

Source that environment once created, then you may install deriver using pip: `pip install deriver`, or by cloning this 
repo and from the cloned directory running: `pip install -e .`

# To-do:
* unit testings and code coverage
* refactor filter module to be better and more polymorphic
* create a means to combine multiple `derive_*` methods into a single function, polymorphically
