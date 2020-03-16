"""
Script to dynamically get the BRICS mating rules, in case they change/update
"""

# this is the reaction definiations
from rdkit.Chem.BRICS import reactionDefs
from collections import defaultdict

# this is the format I want to store them in, a dict of lists
brics_mating_rules = defaultdict(list)

for reaction_set in reactionDefs:
    for reaction in reaction_set:

        # get the two pseudoatom types
        x, y = reaction[0], reaction[1]

        """
        deal with (read: ignore) their special definition of pseudoatom 7
        this causes some problems later, but ultimately is still simpler to deal with
        we need the value to be an integer for the database
        """
        if x == "7a" or x == "7b":
            x = 7
        if y == "7a" or y == "7b":
            y = 7

        # I really must insist that these be integers
        x = int(x)
        y = int(y)

        # store symmetrically
        brics_mating_rules[x].append(y)
        brics_mating_rules[y].append(x)

# remove redundant entries in each list
for key in brics_mating_rules:
    tmp_set = set(brics_mating_rules[key])
    brics_mating_rules[key] = list(tmp_set)
