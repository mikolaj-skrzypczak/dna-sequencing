import numpy as np

from dna_sequencing.graph.graph import Graph
from dna_sequencing.solution.solution_container import SolutionContainer


def find_better_solution_with_smaller_n(graph: Graph, solution: SolutionContainer, n: int) -> SolutionContainer or None:
    oligonucleotide_length = solution.get_oligonucleotide_length()
    overlaps = np.asarray(solution.get_overlaps(), dtype=np.intc)
    mean_overlap = int(np.round(np.mean(overlaps)))

    window_width = int(1 + (n - oligonucleotide_length) / (oligonucleotide_length - mean_overlap))

    iterator = 0
    # TODO
    return None
