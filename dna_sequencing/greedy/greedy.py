import copy

import numpy as np

from dna_sequencing.graph.graph import Graph
from dna_sequencing.solution.solution_container import SolutionContainer


class GreedySolver:
    __graph: np.array
    __visited: set[int]

    def __init__(self, graph: Graph) -> None:
        self.__graph = graph
        self.__visited = set()

    def solve(self) -> SolutionContainer:
        start = np.random.randint(self.__graph.get_vertices_no())
        return copy.deepcopy(self.__solve_greedy(start))

    def __get_sorted_candidates(self, i: int) -> np.array:
        # reverse since best vertices are the ones with highest overlap
        return self.__graph.get_overlaps_of_adjacent_vertices(i).argsort()[::-1]

    def __are_all_candidates_visited(self, candidates: np.array) -> bool:
        return any([_id not in self.__visited for _id in candidates])

    def __get_best_candidate(self, candidates: np.array) -> int or None:
        for i in range(len(candidates)):
            if candidates[i] not in self.__visited:
                return candidates[i]

    def __solve_greedy(self, _id: int) -> SolutionContainer:
        solution = SolutionContainer(self.__graph.get_oligonucleotide_length())
        solution.init_ids(_id)
        solution.init_sequence(self.__graph.get_vertex_by_id(_id).oligonucleotide)

        self.__visited.add(_id)
        candidates = self.__get_sorted_candidates(_id)
        are_there_candidates_left = self.__are_all_candidates_visited(candidates)

        current = _id
        while are_there_candidates_left:
            _next = self.__get_best_candidate(candidates)
            self.__visited.add(_next)

            overlap = self.__graph.get_overlap_between_vertices(current, _next)

            solution.add_vertex_to_solution(
                _id=_next,
                overlap=overlap,
                oligonucleotide=self.__graph.get_vertex_by_id(_next).oligonucleotide
            )

            current = _next
            candidates = self.__get_sorted_candidates(current)
            are_there_candidates_left = self.__are_all_candidates_visited(candidates)

        return solution
