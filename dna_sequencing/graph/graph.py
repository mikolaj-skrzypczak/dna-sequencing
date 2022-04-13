from typing import List
from difflib import SequenceMatcher

from dna_sequencing.graph.vertex import Vertex

import numpy as np


class Graph:
    def __init__(self, test_instance: List[str]):
        self.__oligonucleotides_list = Graph.__create_vertices_list(test_instance)
        self.__vertices_no = len(self.__oligonucleotides_list)
        self.__oligonucleotide_length = len(self.__oligonucleotides_list[0].oligonucleotide)
        self.__lysov_graph_matrix = self.__create_lysov_graph_matrix(self.__oligonucleotides_list)
        print(self.__lysov_graph_matrix)

    @staticmethod
    def __create_vertices_list(_test_instance: List[str]) -> List[Vertex]:
        return [Vertex(i, oligonucleotide) for i, oligonucleotide in enumerate(_test_instance, start=1)]

    @staticmethod
    def __get_edge_value(vertex1: Vertex, vertex2: Vertex) -> int:
        match = SequenceMatcher(None, vertex1.oligonucleotide, vertex2.oligonucleotide).find_longest_match()
        m_len = len(vertex1.oligonucleotide[match.a:match.a + match.size])
        return m_len + 1

    def __create_lysov_graph_matrix(self, oligonucleotides_list) -> np.array:
        return np.array(
            [[self.__get_edge_value(oligonucleotides_list[row], oligonucleotides_list[col])
              for col in range(self.__vertices_no)
              ] for row in range(self.__vertices_no)])
