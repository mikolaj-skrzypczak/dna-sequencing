import numpy as np

from dna_sequencing.graph.graph import Graph
from dna_sequencing.solution.solution_container import SolutionContainer


def __find_starting_index_of_best_subspace(subspace_width: int, overlaps: np.array) -> int:
    best_eval = 0
    best_ind = 0
    for i in range(len(overlaps) + 1 - subspace_width):
        subspace_overlaps = overlaps[i: i + subspace_width]
        if (new_eval := np.sum(subspace_overlaps)) > best_eval:
            best_ind = i
            best_eval = new_eval
    return best_ind


def __create_new_solution(old_solution: SolutionContainer, oligonucleotide_length: int,
                          graph: Graph, starting_index: int, optimal_sequence_length: int) -> SolutionContainer:
    new_solution = SolutionContainer(oligonucleotide_length)
    new_solution.init_ids(old_solution.get_vertex_id_by_index(starting_index))
    new_solution.init_sequence(graph.get_vertex_by_ind(starting_index).oligonucleotide)

    current = starting_index
    while new_solution.get_sequence_length() <= optimal_sequence_length:
        left_vertex = graph.get_vertex_by_ind(old_solution.get_vertex_id_by_index(current))
        right_vertex = graph.get_vertex_by_ind(old_solution.get_vertex_id_by_index(current + 1))
        overlap = graph.compute_overlap(left_vertex, right_vertex)

        if (oligonucleotide_length - overlap + new_solution.get_sequence_length()) > optimal_sequence_length:
            break

        new_solution.add_vertex_to_solution(
            _id=right_vertex.id,
            overlap=overlap,
            oligonucleotide=right_vertex.oligonucleotide
        )
        current += 1

    return new_solution


def shrink_solution_to_fit_optimal_sequence_length(
        graph: Graph, old_solution: SolutionContainer, optimal_sequence_length: int) -> SolutionContainer:
    oligonucleotide_length = old_solution.get_oligonucleotide_length()
    overlaps = np.asarray(old_solution.get_overlaps(), dtype=np.intc)
    mean_overlap = int(np.round(np.mean(overlaps)))

    subspace_width = \
        int(1 + (optimal_sequence_length - oligonucleotide_length) / (oligonucleotide_length - mean_overlap))
    starting_index = __find_starting_index_of_best_subspace(subspace_width, overlaps)

    new_solution = \
        __create_new_solution(old_solution, oligonucleotide_length, graph, starting_index, optimal_sequence_length)
    return new_solution
