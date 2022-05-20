import copy
import random
from math import exp

import numpy as np

from dna_sequencing.graph.graph import Graph
from dna_sequencing.solution.solution_container import SolutionContainer


class SimulatedAnnealingSolver:
    __graph: Graph
    __solution: SolutionContainer
    __iterations: int
    __initial_temperature: int
    __n: int

    def __init__(self, graph: Graph, solution: SolutionContainer,
                 iterations: int, initial_temperature: int, n: int) -> None:
        self.__graph = graph
        self.__solution = solution
        self.__iterations = iterations
        self.__initial_temperature = initial_temperature
        self.__n = n

    def solve(self):
        for iteration in range(self.__iterations):
            _T = self.__update_temperature(iteration)
            neighbor_solution = self.__get_new_solution()
            current_evaluation = self.__solution.evaluate_solution()
            neighbor_evaluation = neighbor_solution.evaluate_solution()
            if self.__acceptance_probability(current_evaluation, neighbor_evaluation, _T) >= random.random():
                self.__solution = neighbor_solution
        # todo handle positive mistakes
        return self.__solution

    @staticmethod
    def __acceptance_probability(current: np.ndarray, neighbor: np.ndarray, t: float) -> float:
        diff = neighbor - current
        if diff < 0:
            return 1
        return exp(-diff / t)

    def __update_temperature(self, iteration: int):
        return self.__initial_temperature * pow(0.85, iteration)

    def __get_new_solution(self) -> SolutionContainer:
        first_ind, second_ind = self.__choose_pair_to_swap()
        new_solution = copy.deepcopy(self.__solution)

        new_solution.swap_in_ids_list(first_ind, second_ind)

        first_id = new_solution.get_vertex_id_by_index(first_ind)
        second_id = new_solution.get_vertex_id_by_index(second_ind)

        if first_ind != 0:
            left_id = new_solution.get_vertex_id_by_index(first_ind - 1)
            self.__recalculate_and_set_new_overlap(new_solution, left_id, first_id, first_ind - 1)

        self.__recalculate_and_set_new_overlap(new_solution, first_id, second_id, first_ind)

        if second_ind + 1 != new_solution.get_vertices_count():
            right_id = new_solution.get_vertex_id_by_index(second_ind + 1)
            self.__recalculate_and_set_new_overlap(new_solution, second_id, right_id, first_ind + 1)

        return new_solution

    def __recalculate_and_set_new_overlap(
            self, solution: SolutionContainer, left_id: int, right_id: int, set_on_ind: int) -> None:
        left_vertex = self.__graph.get_vertex_by_id(left_id)
        right_vertex = self.__graph.get_vertex_by_id(right_id)
        new_overlap = self.__graph.compute_overlap(left_vertex, right_vertex)
        solution.set_overlap(ind=set_on_ind - 1, overlap=new_overlap)

    def __choose_pair_to_swap(self) -> tuple[int, int]:
        overlaps = self.__solution.get_overlaps()
        overlaps_keys_set = set(overlaps) - {self.__solution.get_oligonucleotide_length() - 1}
        # divide count by overlap length to prioritize swapping weekly aligned oligonucleotides:
        overlaps_frequencies_modified = \
            {overlap: overlaps.count(overlap) / (overlap + 1) for overlap in overlaps_keys_set}

        keys = np.asarray(list(overlaps_frequencies_modified.keys()))
        values = np.asarray(list(overlaps_frequencies_modified.values()))
        normalized_values = values / np.linalg.norm(values, ord=1)

        chosen_overlap = np.random.choice(keys, p=normalized_values)
        considered_indexes = np.asarray([i for i in range(len(overlaps)) if overlaps[i] == chosen_overlap])
        chosen_index = np.random.choice(considered_indexes)

        return chosen_index, chosen_index + 1
