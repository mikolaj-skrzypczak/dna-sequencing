import os

from dna_sequencing.graph.graph import Graph
from dna_sequencing.greedy.greedy import GreedySolver
from dna_sequencing.utils.file_handing import read_file

FILEPATH = "test_instances/neg_random/9.200-40.txt"


def main():
    os.chdir("..") if str(os.getcwd()).split("\\")[-1] != "dna-sequencing" else None

    test_instance = read_file(FILEPATH)

    if test_instance:
        graph = Graph(test_instance)
        greedy_solution = GreedySolver(graph).solve()
        print(greedy_solution)

    else:
        print("Test instance not found!")


if __name__ == '__main__':
    main()
