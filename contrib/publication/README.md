dme for replicating the results from "Assessing Methods and Obstacles in Chemical Space Exploration"!

## Getting Started
To benchmark Deriver using the Guacamol framework, you will need an environment with the `guacamol` python package:
```bash
pip install guacamol
```
`guacamol` requires the [RDKit library](http://rdkit.org/) (version `2018.09.1.0` or newer).
The publication uses `deriver` version 2.3.4:
```bash
pip install deriver==2.3.4
```
If you wish to use BRICS generators, you will need to create a BRICS fragment database from the subset of Chembl provided by Guacamol. First download the data:
```bash
wget https://ndownloader.figshare.com/files/13612745 -O guacamol_v1_all.smi
```
Then use Deriver's built-in fragment database generator:
```bash
python deriver/src/deriver/fragment.py -i guacamol_v1_all.smi -o chembl.db
```
The authors note that this process is memory and time-intensive. It is advisable to split the list into many parts, run in parallel on virtual machines, and then recombine the final set of fragments. This was done for the publication.
The authors also advise using the following folder structure when running benchmarks to avoid pathing errors:
.
├ results.py
├ plotting.py
├ data
│   ├ alert_collection.csv
│   ├ rules.json (for using rd_filters to check for passing molecules)
│   └ guacamol_v1_all.smi
└ deriver_goal
    ├ chembl.db
    ├ chembl.bin (created after first run from chembl.db)
    └ goal_directed_generation.py
## Benchmarking Deriver
After this setup, any benchmark described in the paper can be replicated by using the correct arguments described at the end of goal_directed_generation.py:
```bash
nohup python deriver_goal/goal_directed_generation.py
```
`nohup` must be used! `nohup` saves the console output of the script to `nohup.out` which can be later used by plotting software. Many of the results are also saved by CSVs for user analysis.
Note that running the script again in the same location will overwrite results.

## Reproducing Figures
The non-post-processed version of the publication figures require plotly:
```bash
pip install plotly
```
Static svg generation requires `orca`. Please follow the configuration instructions here: https://github.com/plotly/orca
After running `goal_directed_generation.py` there is some post-processing to be done by `results.py`:
```bash
python results.py --folder deriver_goal
```
The authors assert the following folder structure when plotting using plotting.py. Each folder is to contain the contents of deriver_goal for a single benchmark, with nohup.out copied from the parent folder and chembl.db removed:
all_results
├ data
│   ├ alert_collection.csv
│   └ rules.json
├ director
│   ├ BRICS + naive SELFIES + Scanner Greedy
│   ├ BRICS + naive SELFIES + Scanner Linear
│   └ BRICS + naive SELFIES + Scanner Metropolis
├ external
│   ├ BRICS + naive SELFIES + Scanner Greedy
│   ├ CReM
│   ├ Deriver Optimized
│   ├ MSO reported
│   └ graph GA reported
├ generator
│   ├ BRICS + naive SELFIES
│   ├ BRICS + naive SELFIES + Scanner
│   ├ BRICS only
│   ├ Scanner only
│   └ naive SELFIES only
├ graph_GA_discriminator
│   ├ graph GA
│   ├ graph GA counterscreen
│   ├ graph GA delayed filter
│   └ graph GA persistent filter
└ mixed_discriminator
    ├ BRICS + naive SELFIES + Scanner delayed
    ├ BRICS + naive SELFIES + Scanner persistent
    └ BRICS + naive SELFIES + Scanner unfiltered
This will allow a user to reproduce all graphical figures in an unrefined state when calling this from the parent folder:
```bash
python plotting.py all_results
```
The `main()` function of plotting.py demonstrates how individual figures may be generated from single benchmarking experiments or from groups as seen in comparative plots and indicated in the grouped folder structure. The parameters used in each folder shown in the tree correspond to one experiment which can be identified using the plublication.






