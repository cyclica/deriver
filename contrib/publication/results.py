import os
import pandas as pd
from subprocess import call
import re


benchmarks = ['Celecoxib rediscovery', 'Troglitazone rediscovery', 'Thiothixene rediscovery',
              'Aripiprazole similarity', 'Albuterol similarity', 'Mestranol similarity', 'C11H24',
              'C9H10N2O2PF2Cl', 'Median molecules 1', 'Median molecules 2', 'Osimertinib MPO',
              'Fexofenadine MPO', 'Ranolazine MPO', 'Perindopril MPO', 'Amlodipine MPO', 'Sitagliptin MPO',
              'Zaleplon MPO', 'Valsartan SMARTS', 'Scaffold Hop', 'Deco Hop']

selected_lengths = [1, 1, 1, 100, 100, 100, 159, 250, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
top_100_lengths = [100] * 20


def scrub_population_scores(nohup):
    """Use to remove population scores from nohup when they don't need to appear in graph"""
    lines = open(nohup, 'r').readlines()
    lines_out = []
    for line in lines:
        if 'Population mean' not in line:
            lines_out.append(line)
    open(nohup, 'w').writelines(lines_out)


def break_into_many(f_name, rule='selected_output_smiles'):
    """Separate output smiles into many files to use rd_filters"""
    new_path = re.sub(os.path.basename(f_name), rule, f_name)
    os.makedirs(new_path, exist_ok=True)
    if rule == 'selected_output_smiles':
        rules = selected_lengths
    elif rule == 'top_100_output_smiles':
        rules = top_100_lengths
    else:
        raise UnboundLocalError

    with open(f_name, 'r') as f:
        lines = f.readlines()
        line_no = 0
        for name, l in zip(benchmarks, rules):
            try:
                write_lines = ''.join(lines[line_no:line_no + l])
            except IndexError:
                print(f_name + ' failed at ' + name)
                break
            line_no += l
            new_name = re.sub(' ', '_', name) + '_smiles.smi'
            new_name = os.path.join(new_path, new_name)
            with open(new_name, 'w') as w:
                w.write(write_lines)
                w.close()
        f.close()


def describe_failures(f_name):
    """Use alert collection csv to figure out what rule set tripped a failure."""
    alerts = pd.read_csv('data/alert_collection.csv')
    df = pd.read_csv(f_name)
    df.loc[df['FILTER'] != 'OK', 'FILTER'] = df.loc[df['FILTER'] != 'OK', 'FILTER'].apply(
        lambda x: x.split('>')[0][:-1])
    mapper = alerts[['description', 'rule_set_name']].set_index('description').to_dict()['rule_set_name']
    mapper.update({'OK': None})
    df['rule_set_name'] = df['FILTER'].map(mapper)
    return df


def create_experiment_filter_dataframe(folder, kind='selected'):
    """Get a dataframe from the rd_filters filtered molecules describing individual failures"""
    df = None
    for benchmark in benchmarks:
        benchmark = re.sub(' ', '_', benchmark)
        new_path = os.path.join(folder, kind + '_output_smiles', benchmark + '_smiles.smi_f.csv')
        new_df = describe_failures(new_path)
        new_df['benchmark'] = benchmark
        if df is None:
            df = new_df
        else:
            df = pd.concat([df, new_df], axis=0)
    df.to_csv(os.path.join(folder, f'all_{kind}_filtered_smiles.csv'))


def combine_scores_and_filter_passing(folder='deriver_goal'):
    """Use rd_filters to assess filter passing molecules for each benchmark and combine these results with scores"""
    try:
        os.remove('nohup.out')
    except FileNotFoundError:
        pass
    filter_passing = pd.DataFrame(index=benchmarks,
                                  columns=['selected_output_smiles', 'top_100_output_smiles'])

    # Get filter failure percentages
    for suffix in ['selected_output_smiles.smi', 'top_100_output_smiles.smi']:
        entry = suffix.split('.')[0]

        new_name = os.path.join(folder, suffix)
        break_into_many(new_name, entry)

        for benchmark in benchmarks:
            b = re.sub(' ', '_', benchmark)
            benchmark_file_name = os.path.join(folder, entry, b + '_smiles.smi')
            prefix = benchmark_file_name + '_f'
            print(benchmark_file_name)

            c = ['nohup', 'python', '../rd_filters/rd_filters/rd_filters.py', 'filter',
                 '--in', benchmark_file_name,
                 '--prefix', prefix,
                 '--rules', 'data/rules.json',
                 '--alerts', 'data/alert_collection.csv']
            call(c)

            with open('nohup.out', 'r') as f:
                s = f.read()
                filter_passing.at[benchmark, entry] = float(re.search('.*passed filters (.*)%', s)[1])
                f.close()

            os.remove('nohup.out')

    df = pd.read_csv(os.path.join(folder, 'benchmark_scores.csv'), index_col=0)
    df[['frac selected', 'frac top 100']] = filter_passing[['selected_output_smiles', 'top_100_output_smiles']]

    df = df.reset_index()
    create_experiment_filter_dataframe(folder, kind='selected')
    create_experiment_filter_dataframe(folder, kind='top_100')

    print(os.path.join(folder, 'results.csv'))
    df.to_csv(os.path.join(folder, 'results.csv'))


def parse_scores_from_nohup(f_name):
    """Although the CSVs produced here are the same as those this script deletes, they are regenerated as part of the pipeline."""
    print(f_name)
    try:
        os.remove(re.sub('nohup.out', 'mean_scores_by_generation.txt', f_name))
        os.remove(re.sub('nohup.out', 'best_scores_by_generation.txt', f_name))
        os.remove(re.sub('nohup.out', 'worst_scores_by_generation.txt', f_name))
        os.remove(re.sub('nohup.out', 'std_scores_by_generation.txt', f_name))
        #os.remove(re.sub('nohup.out', 'population_mean_scores_by_generation.txt', f_name))
        os.remove(re.sub('nohup.out', 'benchmark_scores.csv', f_name))
    except FileNotFoundError as e:
        print(e)
        pass

    def save_scores_by_gen(scores, kind, location='deriver_goal'):
        with open(os.path.join(location, f'{kind}_scores_by_generation.txt'), 'a+') as f:
            f.write(', '.join([str(s) for s in scores]) + '\n')
            f.close()

    with open(f_name, 'r') as f:
        lines = f.readlines()
        f.close()

    i = 0
    df = pd.DataFrame(index=benchmarks, columns=['score', 'frac selected', 'frac top 100', 'generations',
                                                 'top n mean', 'top n min', 'top 1', 'top n std'])
    benchmark_no = 1
    # skip lines until we get to a new benchmark
    while True:
        mean_scores = []
        best_scores = []
        worst_scores = []
        std_scores = []
        population_scores = []
        while f'Running benchmark {benchmark_no}/20' not in lines[i]:
            if '| max: ' in lines[i]:
                try:
                    current_mean = re.match('.*\| avg: ([0-1]\.\d+) \|.*', lines[i])[1]
                    current_best = re.match('.*\| max: ([0-1]\.\d+) \|.*', lines[i])[1]
                    current_worst = re.match('.*\| min: ([0-1]\.\d+) \|.*', lines[i])[1]
                    current_std = re.match('.*\| std: ([0-1]\.\d+) \|.*', lines[i])[1]
                    generation = re.match('^(\d+) \|.*', lines[i])[1]
                    mean_scores.append(current_mean)
                    best_scores.append(current_best)
                    worst_scores.append(current_worst)
                    std_scores.append(current_std)
                except Exception as e:
                    raise e
            elif 'Score:' in lines[i]:
                df.at[benchmarks[benchmark_no - 2], 'score'] = re.match('.* ([0-1]*[.0-9]+)', lines[i])[1]
                df.at[benchmarks[benchmark_no - 2], 'top 1'] = current_best
                df.at[benchmarks[benchmark_no - 2], 'top n mean'] = current_mean
                df.at[benchmarks[benchmark_no - 2], 'top n min'] = current_worst
                df.at[benchmarks[benchmark_no - 2], 'top n std'] = current_std
                df.at[benchmarks[benchmark_no - 2], 'generations'] = generation
            elif 'Population mean' in lines[i]:
                population_scores.append(re.match('Population mean: (\d+.\d+)$', lines[i])[1])
            i += 1
            try:
                lines[i]
            except IndexError:
                break
        if benchmark_no > 1:
            save_scores_by_gen(mean_scores, 'mean', re.sub('nohup.out', '', f_name))
            save_scores_by_gen(best_scores, 'best', re.sub('nohup.out', '', f_name))
            save_scores_by_gen(worst_scores, 'worst', re.sub('nohup.out', '', f_name))
            save_scores_by_gen(std_scores, 'std', re.sub('nohup.out', '', f_name))
            #save_scores_by_gen(population_scores, 'population_mean', re.sub('nohup.out', '', f_name))
        benchmark_no += 1
        i += 1
        try:
            lines[i]
        except IndexError:
            break
    df.to_csv(re.sub('nohup.out', 'benchmark_scores.csv', f_name))


def tabulate_total_scores_and_generations(folder='/home/shawnreeves/final_results_final_final'):
    df = pd.DataFrame(columns=['folder', 'experiment', 'total_score', 'mean_score'])
    for subfolder in ['external', 'generator', 'graph_GA_discriminator', 'mixed_discriminator', 'director']:
        new_path = os.path.join(folder, subfolder)
        for experiment in os.listdir(new_path):
            if experiment == 'plots' or not os.path.isdir(os.path.join(new_path, experiment)):
                continue
            scores = pd.read_csv(os.path.join(new_path, experiment, 'results.csv'))
            total = sum(scores['score'])
            mean = total/20
            try:
                gens = sum(scores['generations'])
            except:
                gens = None
            df = df.append(pd.Series({'folder': subfolder, 'experiment': experiment,
                                      'total_score': total, 'mean_score': mean, 'total_generations': gens}),
                           ignore_index=True)
    df.to_csv(os.path.join(folder, 'all_scores.csv'))


def tabulate_filter_and_scores(folder='/home/shawnreeves/final_results_final_final/graph_GA_discriminator'):
    inner = ['score', 'frac selected', 'frac top 100', 'top n mean', 'top n min', 'top 1']
    outer = ['graph GA', 'graph GA persistent filter', 'graph GA delayed filter', 'graph GA counterscreen']
    df = pd.DataFrame(index=benchmarks, columns=pd.MultiIndex.from_product([outer, inner])).T
    for o in outer:
        for i in inner:
            vals = pd.read_csv(os.path.join(folder, o, 'results.csv'))
            df.loc[(o, i), :] = list(vals[i].round(3))
    df.index = pd.MultiIndex.from_product([['unfiltered', 'filtered', 'delayed', 'counterscreen'], inner])
    df.T.to_csv('/home/shawnreeves/final_results_final_final/filter_scores.csv')


def calculate_mean_filter_passing(folder):
    filter_passing = pd.DataFrame()
    for sub in os.listdir(folder):
        if sub != 'plots':
            df = pd.read_csv(os.path.join(folder, sub, 'results.csv'))
            filter_passing = filter_passing.append(pd.Series(name=sub, data=list(df[['frac selected', 'frac top 100']].mean())))
    filter_passing.columns = ['selected', 'top_100']
    return filter_passing


def get_results():
    parser = argparse.ArgmentParser()
    parser.add_argument('--folder', default='deriver_goal')
    args = parser.parse_args()
    parse_scores_from_nohup(os.path.join(args.path, 'nohup.out'))
    combine_scores_and_filter_passing(args.path)

if __name__ == "__main__":
    get_results()
