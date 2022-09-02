import os
from operator import itemgetter

import pandas as pd

import util

HEADER = ['Nuc', 'Acid', 'Base', 'D1', 'D2', 'fitness']
HEADER_OCCURRENCES = ['Nuc', 'Acid', 'Base', 'D1', 'D2', 'fitness', 'occurrences']


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


# TODO: KeyError: 'occurrences'
def store_best_individual_occurrences(ga_directory, header, output_directory):
    for filename in os.listdir(ga_directory):
        best_df = pd.read_csv(os.path.join(ga_directory, filename), header=0)
        count_df = best_df.groupby(header).size().reset_index(name='Count')
        count_df.columns = HEADER_OCCURRENCES

        sorted_df = count_df.sort_values(['occurrences', 'fitness'], ascending=[False, False])

        count_df.to_csv(os.path.join(output_directory, filename), header=HEADER_OCCURRENCES)


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
