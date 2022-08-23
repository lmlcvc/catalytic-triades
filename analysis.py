import os
from operator import itemgetter


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
