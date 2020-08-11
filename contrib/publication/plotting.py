import pandas as pd
from plotly import graph_objects as go
from plotly.subplots import make_subplots
import os
import plotly
import re


benchmarks = ['Celecoxib rediscovery', 'Troglitazone rediscovery', 'Thiothixene rediscovery',
              'Aripiprazole similarity', 'Albuterol similarity', 'Mestranol similarity', 'C11H24',
              'C9H10N2O2PF2Cl', 'Median molecules 1', 'Median molecules 2', 'Osimertinib MPO',
              'Fexofenadine MPO', 'Ranolazine MPO', 'Perindopril MPO', 'Amlodipine MPO', 'Sitagliptin MPO',
              'Zaleplon MPO', 'Valsartan SMARTS', 'Scaffold Hop', 'Deco Hop']

benchmarks_cap = ['Celecoxib Rediscovery', 'Troglitazone Rediscovery', 'Thiothixene Rediscovery',
              'Aripiprazole Similarity', 'Albuterol Similarity', 'Mestranol Similarity', 'C11H24',
              'C9H10N2O2PF2Cl', 'Median Molecules 1', 'Median Molecules 2', 'Osimertinib MPO',
              'Fexofenadine MPO', 'Ranolazine MPO', 'Perindopril MPO', 'Amlodipine MPO', 'Sitagliptin MPO',
              'Zaleplon MPO', 'Valsartan SMARTS', 'Scaffold Hop', 'Deco Hop']


def load_and_parse_benchmark_scores(location, generations=False):
    df = pd.read_csv(location)
    if not generations:
        df = df[['index', 'score']]
    else:
        df = df[['index', 'generations']]
    df = df.set_index('index')
    return df


def gather_all_benchmark_results(location='all_results', generations=False):
    df = None
    for folder in os.listdir(location):
        if folder=='plots' or 'counterscreen' in folder or 'reported' in folder:
            continue
        current_df = load_and_parse_benchmark_scores(os.path.join(location, folder, 'results.csv'),
                                                     generations=generations)
        current_df.columns = [folder]
        if isinstance(df, pd.DataFrame):
            df = df.join(current_df)
        else:
            df = current_df
    return df


def gather_all_benchmark_generations(location='all_results'):
    df = None
    for folder in os.listdir(location):
        if folder=='plots':
            continue
        current_df = load_and_parse_benchmark_scores(os.path.join(location, folder, 'results.csv'), generations=True)
        current_df.columns = [folder]
        if isinstance(df, pd.DataFrame):
            df = df.join(current_df)
        else:
            df = current_df
    return df


def graph_compare_benchmarks(df, title, out_loc, space=0.15, height=900, width=1600):
    layout = go.Layout(
        plot_bgcolor='rgba(0,0,0,0)'
    )
    markers = ['diamond', 'x', 'circle', 'triangle-up', 'cross', 'square']
    fig = go.Figure(layout=layout)
    colour_set = ['#8aff00', '#005bed', '#ff0000', '#ff00ff', '#ff7f00']
    x_labels = benchmarks_cap
    deltas = [(-(len(df.columns)-1)*space + 2*space*k) for k in range(len(df.columns))]
    i = 0
    df.to_csv(out_loc + '.csv')
    for idx, column in enumerate(df.columns):

        fig.add_trace(
            go.Scatter(
                name=column,
                mode='markers',
                legendgroup=column,
                #x0=deltas[idx]+1,
                x0=1,
                dx=1,
                y=df[column],
                marker=dict(
                    symbol='line-ew',
                    #color=plotly.colors.qualitative.Set1[idx],
                    size=40,
                    opacity=0.5,
                    line=dict(
                        color=colour_set[idx],
                        width=5
                    ),
                ),
                showlegend=True,
            )
        )

        fig.add_trace(
            go.Scatter(
                name=column,
                mode='markers',
                legendgroup=column,
                x0=deltas[idx]+1,
                #x0=1,
                dx=1,
                y=df[column],
                marker=dict(
                    symbol=markers[idx],
                    color='#000000',
                    size=7,
                    opacity=1,
                    ),
                ),
            )

        i += 1

    # grey rectangular sections for each benchmark
    for j in range(1, len(x_labels)+1):
        fig.add_shape(
            # filled Rectangle
            type="rect",
            x0=j-0.45,
            y0=-0.05,
            x1=j+0.45,
            y1=1.05,
            line=dict(
                width=0,
            ),
            fillcolor='rgba(0.1,.1,.10,.10)',
        )

    # black outline for each category
    for k, l in [(1, 3), (4, 6), (7, 8), (9, 10), (11, 17), (18, 20)]:
        fig.add_shape(
            type="rect",
            x0=k-0.5,
            y0=-0.05,
            x1=l+0.5,
            y1=1.05,
            line=dict(
                width=5,
                color='Black'
            ),
        )

    fig.update_layout(
        yaxis_title='Score',
        xaxis_title='Benchmark',
        title=title,
        height=height,
        width=width
    )
    fig.update_xaxes(showgrid=False, ticks='outside', ticktext=x_labels, tickmode='array', tickvals=list(range(1, 21)),
                     tickangle=45, range=[0.4, 20.6])
    fig.update_yaxes(ticks='outside', range=[-0.05, 1.05])
    fig.write_image(out_loc+'.png')
    fig.write_image(out_loc+'.svg')


def graph_compare_benchmark_generations(df, title, out_loc, space=0.15, height=900, width=1600):
    layout = go.Layout(
        plot_bgcolor='rgba(0,0,0,0)'
    )
    y_max = max(df.to_numpy().flatten())
    colour_set = ['#8aff00', '#005bed', '#ff0000', '#ff00ff', '#ff7f00']
    #markers = ['diamond', 'x', 'circle', 'triangle-up', 'cross', 'square']
    fig = go.Figure(layout=layout)
    x_labels = benchmarks_cap
    #deltas = [(-(len(df.columns)-1)*space + 2*space*k) for k in range(len(df.columns))]
    i = 0
    for idx, column in enumerate(df.columns):
        fig.add_trace(
            go.Scatter(
                name=column,
                mode='markers',
                legendgroup=column,
                #x0=deltas[idx]+1,
                x0=1,
                dx=1,
                y=df[column],
                marker=dict(
                    symbol='line-ew',
                    size=40,
                    opacity=0.5,
                    line=dict(
                        color=colour_set[idx],
                        width=5
                    ),
                ),
                showlegend=True,
            )
        )
        i += 1

    # grey rectangular sections for each benchmark
    for j in range(1, len(x_labels)+1):
        fig.add_shape(
            # filled Rectangle
            type="rect",
            x0=j-0.45,
            y0=0,
            x1=j+0.45,
            y1=y_max+(y_max*0.05),
            line=dict(
                width=0,
            ),
            fillcolor='rgba(0.1,.1,.10,.10)',
        )

    # black outline for each category
    for k, l in [(1, 3), (4, 6), (7, 8), (9, 10), (11, 17), (18, 20)]:
        fig.add_shape(
            type="rect",
            x0=k-0.5,
            y0=0,
            x1=l+0.5,
            y1=y_max+(y_max*0.05),
            line=dict(
                width=5,
                color='Black'
            ),
        )

    fig.update_layout(
        yaxis_title='Score',
        xaxis_title='Benchmark',
        title=title,
        height=height,
        width=width
    )
    fig.update_xaxes(showgrid=False, ticks='outside', ticktext=x_labels, tickmode='array', tickvals=list(range(1, 21)),
                     tickangle=45, range=[0.4, 20.6], showticklabels=True)
    fig.update_yaxes(ticks='outside', range=[0, y_max+(y_max*0.05)], showticklabels=True)
    fig.write_image(out_loc+'.png')
    fig.write_image(out_loc+'.svg')


def load_and_parse_generational_scores(location):
    index = pd.MultiIndex.from_product(
        [pd.Series(benchmarks, name='benchmark'),
         pd.Series(['best', 'worst', 'mean', 'std', 'population_mean'], name='kind'),
         pd.Series([os.path.split(location)[-1]], name='name')])

    df = pd.DataFrame(index=index, columns=range(1, 1001))
    df = df.fillna(0)

    for kind in ['best', 'worst', 'mean', 'std', 'population_mean']:
        try:
            with open(os.path.join(location, f'{kind}_scores_by_generation.txt'), 'r') as f:
                lines = f.readlines()
                if not len(lines) == 20:
                    print(os.path.join(location, f'{kind}_scores_by_generation.txt'), 'does not have 20 lines')
                    continue
                for (bm, line) in zip(benchmarks, lines):
                    vals = line.strip().split(', ')
                    df.loc[(bm, kind, os.path.split(location)[-1]), 1:len(vals)] = vals if len(vals) > 1 else vals[0]
        except FileNotFoundError:
            print(os.path.join(location, f'{kind}_scores_by_generation.txt'), 'does not exist')

    df = df.reindex(pd.MultiIndex.from_tuples(df.index)).reset_index().set_index('level_0')
    return df.rename({'level_1': 'kind', 'level_2': 'experiment'}, axis=1)


def gather_all_generational_results(location='all_results/graph_GA_filtering'):
    df = pd.DataFrame(columns=range(1, 1001))
    for folder in os.listdir(location):
        if 'best_scores_by_generation.txt' in os.listdir(os.path.join(location, folder)):
            print(folder)
            current_df = load_and_parse_generational_scores(os.path.join(location, folder))
            df = pd.concat([df, current_df], axis=0)
    return df


def subplot_compare_generational_results(df, benchmark, experiments, out_loc, n_row, n_col, config='std',
                                         include_population=True, add_zoom=False):
    df = df.loc[benchmark, :]
    df = df.loc[df['experiment'].isin(experiments)]
    df.reset_index().to_csv(out_loc + '.csv')
    title = df.index[0]
    legend = [e[0] for e in df.groupby('experiment')]
    if add_zoom:
       new_legend = []
       for e in legend:
           new_legend.append(e)
           new_legend.append(e)
       legend = new_legend
    #legend = ['<b>'+l+'<\b>' for l in legend]

    fig = make_subplots(n_row, n_col, subplot_titles=legend,
                        shared_xaxes=True, shared_yaxes=True, horizontal_spacing=0.05,
                        vertical_spacing=0.2 if not add_zoom else 0.05)
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)')
    selected_lengths = [100, 100, 100, 100, 100, 100, 159, 250, 100, 100,
                        100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
    n = selected_lengths[benchmarks.index(benchmark)]
    x = list(range(1, 1000))

    row = 1
    col = 1
    max_gen = 0
    max_score = 0
    min_score = 1

    kind_dict = {'best': 'highest score in population',
                 'mean': f'mean score of top {n}',
                 'worst': f'lowest score in top {n}',
                 'std': f'standard deviation in top {n}',
                 'population_mean': 'full population mean'}

    for experiment, data_by_experiment in df.groupby('experiment'):
        j = -1
        grouped = list(data_by_experiment.groupby('kind', sort=True))
        for zoom_flag in [0, 1] if add_zoom else [0]:
            if add_zoom:
                col = 1 + zoom_flag
            if zoom_flag:
                j = -1
            for kind in ['best', 'mean', 'worst', 'population_mean']:
                j += 1
                if (kind in ['best', 'worst', 'mean', 'population_mean']) and config == 'std':
                    continue
                if kind == 'std' and config != 'std':
                    continue
                if kind == 'population_mean' and not include_population:
                    continue
                data = [g[1] for g in grouped if g[0] == kind][0]
                data = list(data[x].values[0].astype(float))
                data = [d for d in data if d > 0]
                if len(data) == 0:
                    print('Skipped', experiment, kind)
                    continue
                max_gen = max(max_gen, len(data))
                max_score = max(max_score, max(data))
                min_score = min(min_score, min(data))
                fig.add_trace(
                        go.Scatter(
                        name=kind_dict[kind],
                        mode='lines',
                        legendgroup=kind_dict[kind],
                        x=x,
                        y=data[:add_zoom] if add_zoom and zoom_flag else data,
                        marker=dict(
                            size=0,
                            color=plotly.colors.qualitative.Set1[j],
                            line=dict(
                                color=plotly.colors.qualitative.Set1[j],
                                width=1.5
                            ),
                        ),
                        showlegend=(row, col) == (1, 1)
                    ), row=row, col=col
                )

            col += 1
            if col > n_col:
                row += 1
                if row > n_row:
                    break
                col = 1

    for row in range(1, n_row+1):
        for col in range(1, n_col+1):
            fig.update_xaxes(showline=True, linecolor='Black', ticks='outside',
                             range=[1, max_gen + 1], row=row, col=col, title='Generation',
                             tickwidth=2, linewidth=2, title_font=dict(color='Black'))
            if add_zoom and col == 2:
                fig.update_xaxes(showline=True, linecolor='Black', ticks='outside',
                                     range=[1, add_zoom +1], row=row, col=2, title='Generation')
            fig.update_yaxes(showline=True, linecolor='Black', ticks='outside', linewidth=2, tickwidth=2,
                             range=[min_score - 0.05, max_score + 0.05], row=row, col=col, title='Score',
                             title_font=dict(color='Black'))

    fig.update_layout(title_text=title, height=900, width=1600, title_font=dict(color='Black'))
    if add_zoom:
        fig.update_layout(title_text=title, height=3000, width=1050)
    if config == 'std':
        out_loc += '_std'
    fig.write_image(out_loc + '.png')
    fig.write_image(out_loc + '.svg')


def main():
    plotly.io.orca.config.use_xvfb = True
    #plotly.io.orca.config.executable = 'orca-1.3.1.AppImage'
    path = 'all_results'
    exclude = ['data', 'get_results.py', 'all_scores.csv', 'filter_scores.csv']#, 'generator', 'discriminator', 'LD_filtering']

    n_rows = [2, 1, 1, 1, 2]
    n_cols = [3, 3, 3, 3, 3]

    i = 0
    folders = ['generator', 'graph_GA_discriminator', 'mixed_discriminator', 'director', 'external']
    for folder in folders:
        if folder in exclude:
            continue
        print(folder)
        new_path = os.path.join(path, folder)
        experiments = os.listdir(new_path)
        try:
            experiments.remove('plots')
        except:
            pass
        os.makedirs(os.path.join(new_path, 'plots'), exist_ok=True)
        benchmark_df = gather_all_benchmark_results(location=new_path)
        graph_compare_benchmarks(benchmark_df, folder, os.path.join(new_path, 'plots', 'compare_benchmarks'),
                                 space=0.5/len(experiments))
        try:
            total_n_gens_df = gather_all_benchmark_results(location=new_path, generations=True)

            graph_compare_benchmark_generations(total_n_gens_df, folder,
                                                os.path.join(new_path, 'plots', 'compare_benchmark_generations'),
                                                space=0.5/len(experiments))
        except:
            print(folder)
            pass

        if folder == 'external':
            continue
        gen_df = gather_all_generational_results(location=new_path)
        for benchmark in benchmarks:
            subplot_compare_generational_results(gen_df, benchmark, experiments=experiments,
                                                 out_loc=os.path.join(new_path, 'plots', benchmark),
                                                 n_row=n_rows[i], n_col=n_cols[i], config='summ',
                                                 include_population=folder == 'graph_GA_filtering')
            if benchmark == 'Ranolazine MPO' and folder == 'generator':
                subplot_compare_generational_results(gen_df, benchmark, experiments=experiments,
                                                 out_loc=os.path.join(new_path, 'plots', benchmark),
                                                 n_row=5, n_col=2, config='summ',
                                                 include_population=folder == 'graph_GA_discriminator',
                                                 add_zoom=10)
            #subplot_compare_generational_results(gen_df, benchmark, experiments=experiments,
            #                                     out_loc=os.path.join(new_path, 'plots', benchmark),
            #                                     n_row=n_rows[i], n_col=n_cols[i], config='std')
        i += 1


if __name__ == '__main__':
    main()
