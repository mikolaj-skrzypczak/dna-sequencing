import os
import re
import time

import numpy as np

from dna_sequencing.graph.graph import Graph
from dna_sequencing.greedy.greedy import GreedySolver
from dna_sequencing.simulated_annealing.simulated_annealing import SimulatedAnnealingSolver
from dna_sequencing.utils.file_handing import read_file

TEST_INSTANCES_FOLDER = 'test_instances'


def __save_str_to_csv(filename: str, csv_str: str) -> None:
    with open(f"{filename}.csv", "w") as f:
        f.write(csv_str)


def __get_original_sequence_length_and_original_nucleotides_count(
        filename: str) -> tuple[int, int]:
    match_obj = re.match("(.*)\.(.*)([+\-])(.*)\.txt", filename)

    before_errors = int(match_obj.group(2))
    error_type = match_obj.group(3)
    errors_count = int(match_obj.group(4))

    original_sequence_length = before_errors + 9
    should_use = before_errors - errors_count if error_type == "-" else before_errors
    return original_sequence_length, should_use


def __get_error_category_to_files_paths_mapping() -> dict[str, list]:
    return {
        error_category: [
            f"{TEST_INSTANCES_FOLDER}/{error_category}/{filename}"
            for filename in os.listdir(f"{TEST_INSTANCES_FOLDER}/{error_category}")
        ] for error_category in os.listdir(TEST_INSTANCES_FOLDER)}


def main() -> None:
    os.chdir("..") if str(os.getcwd()).split("\\")[-1] != "dna-sequencing" else None

    out_csv = "test_instance,oligonucleotides used (best),accuracy (best),sequence length (best)," \
              "average overlap(best), accuracy(average - 10 runs),sequence (best)\n"

    error_categories_to_files = __get_error_category_to_files_paths_mapping()

    for error_category, all_paths in error_categories_to_files.items():
        out_csv += f"{error_category}\n"
        for test_instance_path in all_paths:
            print(test_instance_path)
            best_accuracy = 0
            best_solution = None
            all_accuracies = []

            filename = test_instance_path.split("/")[-1]
            optimal_sequence_length, optimal_nucleotides_count = \
                __get_original_sequence_length_and_original_nucleotides_count(filename)

            test_instance = read_file(test_instance_path)

            graph = Graph(test_instance)
            greedy_solution = GreedySolver(graph).solve()
            _time = 0
            for _ in range(10):
                start = time.time()
                sa_solution = SimulatedAnnealingSolver(
                    graph=graph,
                    initial_solution=greedy_solution,
                    iterations=1000,
                    initial_temperature=50,
                    original_sequence_length=optimal_sequence_length
                ).solve()
                end = time.time()
                _time += end - start
                accuracy = sa_solution.get_accuracy(optimal_nucleotides_count)
                all_accuracies.append(accuracy)
                if accuracy > best_accuracy:
                    best_solution = sa_solution
                    best_accuracy = accuracy

            out_csv += f"{best_solution.to_csv_str(filename, optimal_nucleotides_count)}," \
                       f"{np.mean(all_accuracies)},{_time/10},{best_solution.get_sequence()}\n"

    __save_str_to_csv("results", out_csv)


if __name__ == '__main__':
    main()
