#Author Shawn Reeves
#Contact sareeves@uwaterloo.ca
from __future__ import print_function

import argparse
import json
import logging
from typing import List, Optional
from rdkit import rdBase, Chem
from joblib import delayed
import joblib
import gc
import pandas as pd
import random

from guacamol.utils.chemistry import canonicalize
from guacamol.assess_goal_directed_generation import assess_goal_directed_generation
from guacamol.goal_directed_generator import GoalDirectedGenerator
from guacamol.scoring_function import ScoringFunction

from deriver.api import Deriver
from deriver.child_filter import apply_filter
from deriver.fragment_index import frag_index_cache  #pylint: disable=unused-import
import numpy as np
import os
from guacamol.utils.chemistry import canonicalize_list

# !!!
# publication uses deriver v 2.3.4
# !!!

def make_mating_pool(scored_population, selection_size, method, best, temperature):
    """Select molecules from the population to derive from"""

    print('Sampling from population')
    scored_population = sorted(scored_population, key=lambda x: x[1], reverse=True)

    # this method has a linear relationship between score and probability of selection. It is used in graph_GA
    if method in ['original', 'probs', 'linear']:
        sum_scores = sum([s[1] for s in scored_population])
        population_probs = np.divide([s[1] for s in scored_population], sum_scores)
        indices = np.random.choice(range(len(scored_population)), p=population_probs, size=selection_size, replace=True)
        scored_seeds = [scored_population[i] for i in indices]

    # include all molecules which beat the previous score, and probabilistically sample the remainder to return
    elif method == 'metropolis':
        assert best is not None and temperature is not None
        # get unnormalized selection probabilities
        probs = [min([1, p]) for p in np.exp(-1 * (np.subtract(best, [s[1] for s in scored_population])) / temperature)]

        # check how many children are going to be selected for having the same or better score than the previous best
        n_better_children = sum(np.equal(probs, 1))
        print('Greater or equal children: ',  n_better_children)
        if n_better_children >= selection_size:
            # greed
            scored_seeds = scored_population[:selection_size]
        else:
            # get those scoring 1
            scored_seeds = scored_population[:n_better_children]

            # get those not scoring 1
            lesser_children = scored_population[n_better_children:]
            # and their probs
            lesser_probs = probs[n_better_children:]

            # normalize
            try:
                normalizer = 1/sum(lesser_probs)
                probs_normalized = np.multiply(lesser_probs, normalizer)
            except ZeroDivisionError:
                print('No children with score above 0. Sampling greedily')
                scored_seeds = scored_population[:selection_size]
                scored_seeds = sorted(scored_seeds, key=lambda x: x[1], reverse=True)
                print(f'Sampled {len(scored_seeds)} mols. {len(set(scored_seeds))} unique')
                return scored_seeds

            # sample without replacement accoring to normalized probs, indexing lesser children
            try:
                sampled_children = np.random.choice(range(len(lesser_children)),
                                                    p=probs_normalized,
                                                    size=selection_size - n_better_children,
                                                    replace=False)
                print(f'Sampled {len(sampled_children)} lesser children probabilistically')
            # happens when fewer p > 0 than size, for fewer possible children than size
            except ValueError as e:
                print(e)
                # the lesser of the remaining requested candidates and the remaining available candidates
                sampled_children = range(min([selection_size-n_better_children, len(lesser_children)]))
                print(f'Sampled {len(sampled_children)} lesser children greedily')

            # get their actual (smiles, score) tuples
            chosen_children = [lesser_children[c] for c in sampled_children]
            # extend better children
            scored_seeds = scored_seeds + chosen_children
            # sanity check
            assert len(scored_seeds) == min([selection_size, len(scored_population)])

    # this is the 'greedy' case
    else:
        scored_seeds = scored_population[:selection_size]

    scored_seeds = sorted(scored_seeds, key=lambda x: x[1], reverse=True)
    print(f'Sampled {len(scored_seeds)} mols. {len(set(scored_seeds))} unique')
    return scored_seeds


def summarize_results(generation, population_scores):
    p_max = np.max(population_scores)
    p_avg = np.mean(population_scores)
    p_min = np.min(population_scores)
    p_std = np.std(population_scores)
    p_sum = np.sum(population_scores)

    print(f'{generation} | '
          f'max: {p_max:.3f} | '
          f'avg: {p_avg:.3f} | '
          f'min: {p_min:.3f} | '
          f'std: {p_std:.3f} | '
          f'sum: {p_sum:.3f}')
    return p_max, p_avg, p_min, p_std, p_sum


def save_scores_by_gen(scores, kind, location='deriver_goal'):
    with open(os.path.join(location, f'{kind}_scores_by_generation.txt'), 'a+') as f:
        f.write(', '.join([str(s) for s in scores]) + '\n')
        f.close()


def save_top_mols(children, location, number=100):
    with open(location, 'a+') as f:
        for i in range(number):
            try:
                f.write(children[i] + f' mol{i}\n')
            except IndexError:
                print(f'Less than {number} children were saved')
                break
        f.close()


def save_and_exit(scored_population, number_molecules, mean_scores_by_gen, best_scores_by_gen, worst_scores_by_gen):

    best_children = [s[0] for s in scored_population]
    scores = [s[1] for s in scored_population]
    mean_scores_by_gen.append(np.mean(scores))
    best_scores_by_gen.append(max(scores))
    print('entire final population')
    print(best_children)
    # for doing 'quality' the way guacamol does
    save_top_mols(best_children, location='deriver_goal/top_100_output_smiles.smi', number=100)
    # for capturing the entire output
    save_top_mols(best_children, location='deriver_goal/all_output_smiles.smi', number=len(best_children))
    # for capturing the molecules which will be returned
    save_top_mols(best_children, location='deriver_goal/selected_output_smiles.smi', number=number_molecules)
    save_scores_by_gen(best_scores_by_gen, 'best')
    save_scores_by_gen(mean_scores_by_gen, 'mean')
    save_scores_by_gen(worst_scores_by_gen, 'worst')
    relevant_scores = scores[:number_molecules]
    relevant_children = best_children[:number_molecules]
    print(f"Return len: {len(relevant_scores)}")
    print('children returned to guacamol')
    print(relevant_children)

    return relevant_children


def derive(deriver, seeds, mut_rate, n_brics, n_selfies, n_smiles_gb, n_selfies_gb, scanner):
    print('Deriving new mols...')
    all_mols = set()
    deriver.set_seeds(seeds)

    if scanner:
        good_scanner, _ = deriver.scan_selfies()
        print(f'Generated {len(good_scanner)} scanner mols.')
        all_mols.update(good_scanner)

    if n_brics > 0:
        try:
            good_brics, _ = deriver.derive_brics(n_children=int(n_brics))
            print(f'Generated {len(good_brics)} brics mols.')
            all_mols.update(good_brics)
        except ZeroDivisionError:
            print(f'No valid seed fragments could be generated for BRICs.')

    if n_selfies > 0:
            good_selfies, _ = deriver.derive_selfies(n_children=int(n_selfies))
            print(f'Generated {len(good_selfies)} selfies mols.')
            all_mols.update(good_selfies)

    if n_selfies_gb > 0:
        good_selfies_gb, _ = deriver.derive_gb(n_children=int(n_selfies_gb),
                                               mut_rate=mut_rate, kind='selfies')
        print(f'Generated {len(good_selfies_gb)} selfies_gb mols.')
        all_mols.update(good_selfies_gb)

    if n_smiles_gb > 0:
        good_smiles_gb, _ = deriver.derive_gb(n_children=int(n_smiles_gb),
                                              mut_rate=mut_rate, kind='smiles')
        print(f'Generated {len(good_smiles_gb)} smiles_gb mols.')
        all_mols.update(good_smiles_gb)

    return canonicalize_list(list(all_mols), False)


rdBase.DisableLog('rdApp.error')


class DeriverGenerator(GoalDirectedGenerator):
    """Class which inherits a benchmarking class. Sets up parameters to run generate_optimized_molecules, a class req"""

    def __init__(self,
                 smi_file,
                 generations,
                 population=100,
                 selection_size=200,
                 selection_method='linear',
                 derive_size=100,
                 brics_fragment_db=None,
                 selfies_proportion=0,
                 brics_proportion=0,
                 smiles_gb_proportion=0,
                 selfies_gb_proportion=0,
                 mutation_rate=0.5,
                 enable_scanner=False,
                 enable_filter=False,
                 delayed_filtering=0,
                 patience=5,
                 temperature=1,
                 temp_decay=0.95,
                 start_task=0,
                 random_start=False,
                 derive_population=False,
                 counterscreen=False):

        self.smi_file = smi_file
        # won't be parallel on small vms, but will still work
        self.pool = joblib.Parallel(n_jobs=-1)
        self.generations = generations
        self.selection_size = selection_size
        self.selection_method = selection_method
        self.derive_size = derive_size
        self.population = population
        self.patience = patience
        self.brics_fragment_db = brics_fragment_db
        self.brics_proportion = brics_proportion
        self.selfies_proportion = selfies_proportion
        self.selfies_gb_proportion = selfies_gb_proportion
        self.smiles_gb_proportion = smiles_gb_proportion
        self.mut_rate = mutation_rate
        self.enable_scanner = enable_scanner
        self.enable_filter = enable_filter
        self.delayed_filtering = delayed_filtering
        # task is used to skip a step where a list is scored, opting to check saved scores instead
        self.task = 0
        # start at this benchmarking task, skipping previous ones
        self.start_task = start_task
        # use the population size to determine how many molecules to deriver, as is done in graph_GA in Guacamol
        self.derive_population = derive_population
        # return only members of the population which pass filters
        self.counterscreen = counterscreen

        # unused params
        self.initial_temperature = temperature
        self.temperature_decay = temp_decay
        self.random_start = random_start

    def load_smiles_from_file(self, smi_file):
        with open(smi_file) as f:
            # canonicalize smiles from file and read into list
            return self.pool(delayed(canonicalize)(s.strip()) for s in f)

    def get_precomputed_scores(self, k):
        """If we have already ranked and scored the molecules for a single task, load them to save computation"""
        with open(f'data/ranked_list_task_{self.task}.csv', 'r') as f:
            score_list = f.readlines()
            scored_smiles = [(s[:-1].split(', ')[0], float(s[:-1].split(', ')[1])) for s in score_list[:k]]
            return scored_smiles[:k]

    def filter_mols(self, scored_population, number_mols, add_bad=True):
        print('Filtering final mols for quality...')
        best_children = [Chem.MolFromSmiles(s[0]) for s in scored_population]
        all_filtered_children = apply_filter(self.deriver.data.filter_params, best_children,
                                             must_have_patterns=self.deriver.data.must_have_patterns,
                                             must_not_have_patterns=self.deriver.data.must_not_have_patterns)
        good_children = []
        for child in all_filtered_children:
            if all_filtered_children[child]["is_good"]:
                good_children.append(child)

        return_mols = good_children
        print(f'Filtered mols: {len(best_children)}')
        print(f'Good mols: {len(return_mols)}')
        print(f'Bad mols: {len(best_children)-len(return_mols)}')

        if add_bad:
            i = 0
            while len(return_mols) < number_mols and i < len(scored_population):
                if scored_population[i][0] not in good_children:
                    print('Adding a filter-failing molecule to the output, since there were not enough filtered molecules')
                    return_mols.append(scored_population[i][0])
                i += 1

            scored_return_mols = [smile for (smile, score) in scored_population if smile in return_mols]

        else:
            scored_return_mols = [smile for (smile, score) in scored_population if smile in return_mols][:number_mols]

        print(f'Scored return mols: {len(scored_return_mols)}')
        print(f'Fraction passing filter: {min(len(good_children)/number_mols, 1)}')
        return scored_return_mols

    def rank_and_score(self, smiles, scoring_function):
        print('Scoring population...')
        joblist = (delayed(scoring_function.score)(s) for s in smiles)
        scores = self.pool(joblist)
        scored_smiles = list(zip(smiles, scores))
        scored_smiles = sorted(scored_smiles, key=lambda x: x[1], reverse=True)
        return scored_smiles

    def init_deriver(self):
        # Get a fresh deriver object with nothing saved
        self.deriver = Deriver()
        # Get a basic drug-likeness filter
        if self.enable_filter:
            self.deriver.enable_and_expand_filter()
            alerts = pd.read_csv('data/alert_collection.csv')
            sure_chembl = set(alerts.loc[alerts['rule_set_name'] == 'SureChEMBL', 'smarts'])
            bai = set(alerts.loc[alerts['rule_set_name'] == 'BAI', 'smarts'])
            self.deriver.set_must_not_have_patterns(list(sure_chembl.union(bai)))

        if self.brics_fragment_db:
            self.deriver.set_fragment_source_db(self.brics_fragment_db)

    def clean_up(self):
        if self.brics_proportion > 0:
            global frag_index_cache
            frag_index_cache['deriver_goal/chembl'].clear_cache()
        gc.collect()

    def generate_optimized_molecules(self, scoring_function: ScoringFunction, number_molecules: int,
                                     starting_population: Optional[List[str]] = None) -> List[str]:
        """The function called by the benchmarking software: All backend has to be controlled here"""
        # starting population is provided for some benchmarks.
        if number_molecules > self.population:
            self.population = number_molecules
            print(f'Benchmark requested more molecules than expected: new population is {number_molecules}')

        self.task += 1
        if self.task < self.start_task:
            return ['CCC']
        self.init_deriver()

        scored_population = []
        if starting_population is None:
            print('selecting initial population...')
            if self.random_start:
                all_smiles = self.load_smiles_from_file(self.smi_file)
                selected_smiles = np.random.choice(all_smiles, self.population)
                scored_population = [(s, scoring_function.score(s)) for s in selected_smiles]
            else:
                # we are just going to get the top scoring mols, and we've checked before so we'll just load a file
                scored_population = self.get_precomputed_scores(self.population)
            self.deriver.set_seeds([s[0] for s in scored_population])
        elif len(starting_population) == 1:
            self.deriver.set_seeds(starting_population)
            scored_population = [(s, scoring_function.score(s)) for s in starting_population]

        # allow self-mating in deriver for new methods
        #if len(scored_population) == 1:
        #    scored_population = [scored_population[0], scored_population[0]]

        best = max([s[1] for s in scored_population])
        p_max = best
        no_progess_counter = 0
        old_avg = 0
        mean_scores_by_gen = []
        best_scores_by_gen = []
        worst_scores_by_gen = []
        temperature = self.initial_temperature
        early_stop_annealing = False
        anneal_counter = 0
        filter_enabled = False
        current_population = self.population
        if self.derive_population:
            self.derive_size = self.population

        for generation in range(self.generations):
            # filter annealing #######################################################################
            if ((generation >= (self.generations * self.delayed_filtering)) and self.delayed_filtering) or early_stop_annealing:
                if (not anneal_counter) and (not filter_enabled):
                    filter_enabled = True
                    print(f'Enabling filter at generation {generation}')
                    self.deriver.enable_and_expand_filter()
                    alerts = pd.read_csv('data/alert_collection.csv')
                    sure_chembl = set(alerts.loc[alerts['rule_set_name'] == 'SureChEMBL', 'smarts'])
                    bai = set(alerts.loc[alerts['rule_set_name'] == 'BAI', 'smarts'])
                    self.deriver.set_must_not_have_patterns(list(sure_chembl.union(bai)))
                    print('Expanding population to introduce filtered candidates')
                    current_population *= 2
                print('Generating filtered candidates...')
                anneal_counter += 1
            ###########################################################################################

            scored_seeds = make_mating_pool(scored_population=scored_population,
                                            selection_size=self.selection_size,
                                            method=self.selection_method,
                                            best=best,
                                            temperature=temperature)

            # we want the best score from the previous generation, else there can only be ties
            best = p_max

            good_children = derive(deriver=self.deriver,
                                   seeds=[s[0] for s in scored_seeds],
                                   mut_rate=self.mut_rate,
                                   n_brics=self.derive_size*self.brics_proportion,
                                   n_selfies=self.derive_size*self.selfies_proportion,
                                   n_smiles_gb=self.derive_size*self.smiles_gb_proportion,
                                   n_selfies_gb=self.derive_size*self.selfies_gb_proportion,
                                   scanner=self.enable_scanner)

            scored_children = self.rank_and_score(good_children, scoring_function)
            scored_population = list((set(scored_children)).union(set(scored_population)))
            scored_population = sorted(scored_population, key=lambda x: x[1], reverse=True)[:current_population]
            relevant_scores = [s[1] for s in scored_population][:max([100, number_molecules])]

            # summarization
            p_max, p_avg, p_min, p_std, p_sum = summarize_results(generation, relevant_scores)
            mean_scores_by_gen.append(p_avg)
            best_scores_by_gen.append(p_max)
            worst_scores_by_gen.append(p_min)

            if early_stop_annealing:
                p_avg = np.mean([s[1] for s in scored_population])
                print(f'Population mean: {p_avg}')
            else:
                print(f'Population mean: {np.mean([s[1] for s in scored_population])}')

            # early stopping
            if p_avg == old_avg:
                no_progess_counter += 1
            else:
                no_progess_counter = 0
                anneal_counter = 0
            if self.task < 4 and max(relevant_scores) == 1:
                if self.delayed_filtering:
                    early_stop_annealing = True
                else:
                    print('Finished early on a rediscovery benchmark!')
                    break
            if (no_progess_counter >= self.patience) or (p_avg == 1 and len(scored_population) > 1):
                if self.delayed_filtering:
                    early_stop_annealing = True
                else:
                    print("Finished early!")
                    break
            if (anneal_counter == self.patience) and self.delayed_filtering:
                print("Converged after filtering!")
                break

            old_avg = p_avg
            temperature *= self.temperature_decay
            self.clean_up()

        if self.delayed_filtering:
            scored_population = self.filter_mols(scored_population, number_molecules)
        if self.counterscreen:
            self.deriver.enable_and_expand_filter()
            alerts = pd.read_csv('data/alert_collection.csv')
            sure_chembl = set(alerts.loc[alerts['rule_set_name'] == 'SureChEMBL', 'smarts'])
            bai = set(alerts.loc[alerts['rule_set_name'] == 'BAI', 'smarts'])
            self.deriver.set_must_not_have_patterns(list(sure_chembl.union(bai)))
            scored_population = self.filter_mols(scored_population, number_molecules, add_bad=False)

        return save_and_exit(scored_population, number_molecules, mean_scores_by_gen, best_scores_by_gen, worst_scores_by_gen)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smiles_file', default='data/guacamol_v1_all.smi')
    parser.add_argument('--seed', type=int, default=24)
    parser.add_argument('--output_dir', type=str, default=None)
    parser.add_argument('--suite', default='v2')

    parser.add_argument('--generations', type=int, default=1000)
    parser.add_argument('--population', type=int, default=100)
    parser.add_argument('--selection_size', type=int, default=200)
    parser.add_argument('--selection_method', type=str, default='linear')
    parser.add_argument('--derive_size', type=int, default=100)

    parser.add_argument('--brics_fragment_db', type=str, default='deriver_goal/chembl.db')
    parser.add_argument('--selfies_proportion', type=float, default=0)
    parser.add_argument('--brics_proportion', type=float, default=0)
    parser.add_argument('--selfies_gb_proportion', type=float, default=0)
    parser.add_argument('--smiles_gb_proportion', type=float, default=0)
    parser.add_argument('--mutation_rate', type=float, default=0.5)

    parser.add_argument('--enable_scanner', action='store_true', default=False)
    parser.add_argument('--enable_filter', action='store_true', default=False)
    parser.add_argument('--delayed_filtering', type=float, default=0)

    parser.add_argument('--patience', type=int, default=5)

    parser.add_argument('--temperature', type=float, default=1)
    parser.add_argument('--temp_decay', type=float, default=0.95)
    parser.add_argument('--start_task', type=int, default=0)
    parser.add_argument('--random_start', action='store_true', default=False)

    parser.add_argument('--derive_population', action='store_true', default=False)
    parser.add_argument('--counterscreen', action='store_true', default=False)

    args = parser.parse_args()

    #assert os.path.exists('data/guacamol_v1_all.smiles')

    try:
        os.remove('deriver_goal/all_output_smiles.smi')
        os.remove('deriver_goal/mean_scores_by_generation.txt')
        os.remove('deriver_goal/best_scores_by_generation.txt')
        os.remove('deriver_goal/worst_scores_by_generation.txt')
        os.remove('deriver_goal/selected_output_smiles.smi')
        os.remove('deriver_goal/top_100_output_smiles.smi')
    except FileNotFoundError:
        pass

    np.random.seed(args.seed)
    random.seed = args.seed

    # turn off the millions of lines of useless text
    logging.basicConfig(format='%(levelname)s : %(message)s', level=logging.CRITICAL)

    if args.output_dir is None:
        args.output_dir = os.path.dirname(os.path.realpath(__file__))

    # save command line args
    with open(os.path.join(args.output_dir, 'goal_directed_params.json'), 'w') as jf:
        json.dump(vars(args), jf, sort_keys=True, indent=4)

    optimiser = DeriverGenerator(smi_file=args.smiles_file,
                                 generations=args.generations,
                                 population=args.population,
                                 selection_size=args.selection_size,
                                 selection_method=args.selection_method,
                                 derive_size=args.derive_size,
                                 brics_fragment_db=args.brics_fragment_db,
                                 selfies_proportion=args.selfies_proportion,
                                 brics_proportion=args.brics_proportion,
                                 selfies_gb_proportion=args.selfies_gb_proportion,
                                 smiles_gb_proportion=args.smiles_gb_proportion,
                                 mutation_rate=args.mutation_rate,
                                 enable_scanner=args.enable_scanner,
                                 enable_filter=args.enable_filter,
                                 delayed_filtering=args.delayed_filtering,
                                 patience=args.patience,
                                 temperature=args.temperature,
                                 temp_decay=args.temp_decay,
                                 start_task=args.start_task,
                                 random_start=args.random_start,
                                 derive_population=args.derive_population,
                                 counterscreen=args.counterscreen
                                 )

    json_file_path = os.path.join(args.output_dir, 'goal_directed_results.json')
    assess_goal_directed_generation(optimiser, json_output_file=json_file_path, benchmark_version=args.suite)


if __name__ == "__main__":
    main()
