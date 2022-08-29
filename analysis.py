import os
from operator import itemgetter

import pandas as pd

import util

HEADER = ['Nuc', 'Acid', 'Base', 'D1', 'D2', 'fitness']
HEADER_OCCURRENCES = ['Nuc', 'Acid', 'Base', 'D1', 'D2', 'fitness', 'occurences']


def store_triad_count(directory, output_directory, similar=True):
    count = []
    output_file = open(os.path.join(output_directory, "triad_count.csv"), 'w+')

    HEADER = ["name", "triad_count"]
    if similar:
        HEADER.append("similar_count")

    for filename in os.listdir(directory):
        if not filename.startswith("similar_"):
            line = []

            file = open(os.path.join(directory, filename), 'r')

            line.append(filename.replace(".csv", ""))
            line.append(str(sum(1 for _ in file)))

            if similar:
                similar_file = open(os.path.join(directory, "similar_" + filename), 'r')
                line.append(str(sum(1 for _ in similar_file)))

            count.append(line)

    count = sorted(count, key=itemgetter(0))
    count_text = [','.join(v) for v in count]

    output_file.write(f"{','.join(HEADER)}\n")
    output_file.write('\n'.join(count_text))


def store_best_individual_occurrences(ga_directory, population_directory, output_directory):
    for filename in os.listdir(ga_directory):
        best_df = pd.read_csv(os.path.join(ga_directory, filename), header=0)
        populations_df = pd.DataFrame(columns=HEADER)
        result_df = pd.DataFrame(columns=HEADER_OCCURRENCES)

        for population_file in os.listdir(population_directory):
            if filename.strip(".csv") not in population_file:
                continue
            else:
                iteration_csv = pd.read_csv(os.path.join(population_directory, population_file), header=0)
                iteration_csv = iteration_csv.drop_duplicates(keep='first')

                populations_df = populations_df.append(iteration_csv)

        for index, row in best_df.iterrows():
            row_df = pd.DataFrame(row).T
            row_df.columns = HEADER

            iteration_result = row_df.iloc[0].values.tolist()

            occurrences = row_df.merge(populations_df,
                                       on=HEADER,
                                       how='inner').shape[0]

            iteration_result.append(occurrences)

            iteration_result_df = pd.DataFrame(iteration_result).T
            iteration_result_df.columns = HEADER_OCCURRENCES
            result_df = result_df.append(iteration_result_df)
        result_df.to_csv(os.path.join(output_directory, filename), header=HEADER_OCCURRENCES)


def store_similarity(lists, output_file):
    similarity_lists = []

    for name in lists.keys():
        triad_list = lists[name]
        similarity = util.list_similarities(triad_list)

        similarity_lists.append([name, similarity])

    similarity_df = pd.DataFrame(similarity_lists, columns=['name', 'similarity'])
    similarity_df = similarity_df.sort_values(by=['name'])
    similarity_df.to_csv(output_file)


def store_similarity_best(directory, output_directory):
    lists = {}

    for filename in os.listdir(directory):
        df = pd.read_csv(os.path.join(directory, filename), header=0)
        df = df.drop(columns=['fitness'])
        triad_list = df.values.tolist()

        lists[filename.strip(".csv")] = triad_list

    store_similarity(lists, os.path.join(output_directory, "best.csv"))


def store_similarity_population(directory, output_directory):
    lists = {}

    for filename in os.listdir(directory):
        df = pd.read_csv(os.path.join(directory, filename), header=0)
        df = df.head(10)
        df = df.drop(columns=['fitness'])
        triad_list = df.values.tolist()

        lists[filename.strip(".csv")] = triad_list

    store_similarity(lists, os.path.join(output_directory, "population.csv"))
