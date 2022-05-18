import numpy as np

from dna_sequencing.graph.vertex import Vertex


class Graph:
    __vertices_list: list[Vertex]
    __vertices_no: int
    __oligonucleotide_length: int
    __lysov_graph_matrix: np.array

    def __init__(self, test_instance: list[str]) -> None:
        self.__vertices_list = Graph.__create_vertices_list(test_instance)
        self.__vertices_no = len(self.__vertices_list)
        self.__oligonucleotide_length = len(self.__vertices_list[0].oligonucleotide)
        self.__lysov_graph_matrix = self.__create_lysov_graph_matrix(self.__vertices_list)

    def get_overlap_between_vertices(self, i: int, j: int) -> int:
        return self.__lysov_graph_matrix[i][j]

    def get_overlaps_of_adjacent_vertices(self, i: int) -> np.array:
        return self.__lysov_graph_matrix[i, :]

    @staticmethod
    def __create_vertices_list(_test_instance: list[str]) -> list[Vertex]:
        return [Vertex(i, oligonucleotide) for i, oligonucleotide in enumerate(_test_instance, start=1)]

    def __compute_edge_value(self, vertex1: Vertex, vertex2: Vertex) -> int:
        i = 0
        while i < self.__oligonucleotide_length:
            if vertex2.oligonucleotide.startswith(vertex1.oligonucleotide[i:]):
                break
            i += 1
        return self.__oligonucleotide_length - i

    def __create_lysov_graph_matrix(self, oligonucleotides_list: list[Vertex]) -> np.array:
        return np.array(
            [[self.__compute_edge_value(oligonucleotides_list[row], oligonucleotides_list[col]) if row != col else 0
              for col in range(self.__vertices_no)
              ] for row in range(self.__vertices_no)])

    def get_oligonucleotide_length(self) -> int:
        return self.__oligonucleotide_length

    def get_vertex_by_id(self, _id: int) -> Vertex:
        return self.__vertices_list[_id]

    def get_vertices_no(self) -> int:
        return self.__vertices_no

    def __str__(self) -> str:
        out_str = ""
        for row in range(self.__vertices_no):
            for col in range(self.__vertices_no):
                out_str += f"{self.__lysov_graph_matrix[row][col]}\t"
            out_str += "\n"
        return out_str
