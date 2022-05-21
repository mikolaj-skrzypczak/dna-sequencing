import os

from dna_sequencing.graph.graph import Graph
from dna_sequencing.greedy.greedy import GreedySolver
from dna_sequencing.simulated_annealing.simulated_annealing import SimulatedAnnealingSolver
from dna_sequencing.utils.file_handing import read_file

FILEPATH = "test_instances/pos_random/9.200+80.txt"


def main():
    os.chdir("..") if str(os.getcwd()).split("\\")[-1] != "dna-sequencing" else None

    test_instance = read_file(FILEPATH)

    if test_instance:
        graph = Graph(test_instance)
        greedy_solution = GreedySolver(graph).solve()
        print("Greedy:")
        print(greedy_solution)
        print("evaluation", greedy_solution.evaluate_solution())
        sa_solution = SimulatedAnnealingSolver(
            graph=graph,
            solution=greedy_solution,
            iterations=1000,
            initial_temperature=50,
            optimal_sequence_length=209
        ).solve()
        print("SA:")
        print(sa_solution)
        print("evaluation: ", sa_solution.evaluate_solution())
    else:
        print("Test instance not found!")


if __name__ == '__main__':
    main()
