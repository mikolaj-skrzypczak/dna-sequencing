import os

from dna_sequencing.graph.graph import Graph
from typing import List


def read_lines(filepath: str) -> List[str]:
    try:
        with open(f"{os.getcwd()}/{filepath}", "r") as f:
            return f.read().splitlines()
    except FileNotFoundError:
        print("Given file does not exist")
        return []


def main(append=True):
    filepath = "test_instances/neg_random/9.200-40.txt"

    if append:
        filepath = "dna_sequencing/" + filepath

    test_instance = read_lines(filepath)

    if test_instance:
        graph = Graph(test_instance)


if __name__ == '__main__':
    main(append=False)
