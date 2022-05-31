import copy

import numpy as np


class SolutionContainer:
    __oligonucleotide_length: int
    __vertices_ids: list[int]
    __overlaps: list[int]
    __sequence: str

    def __init__(self, oligonucleotide_length: int) -> None:
        self.__oligonucleotide_length = oligonucleotide_length
        self.__vertices_ids = []
        self.__overlaps = []
        self.__sequence = ""

    def evaluate_solution(self) -> np.ndarray:  # actually np.intc
        actual_overlaps = np.asarray(self.__overlaps, dtype=np.intc)
        ideal_overlaps = np.full(
            shape=actual_overlaps.shape,
            fill_value=self.__oligonucleotide_length - 1,
            dtype=actual_overlaps.dtype
        )
        return np.sum(ideal_overlaps - actual_overlaps)

    def init_ids(self, _id: int) -> None:
        self.__vertices_ids = [_id]

    def init_sequence(self, oligonucleotide: str):
        self.__sequence = oligonucleotide

    def swap_in_ids_list(self, first_ind: int, second_ind: int) -> None:
        (self.__vertices_ids[first_ind], self.__vertices_ids[second_ind]) =\
            (self.__vertices_ids[second_ind], self.__vertices_ids[first_ind])

    def set_overlap(self, ind: int, overlap: int) -> None:
        self.__overlaps[ind] = overlap

    def add_vertex_to_solution(
            self, _id: int, overlap: int, oligonucleotide: str
    ):
        self.__vertices_ids.append(_id)
        self.__overlaps.append(overlap)
        self.__sequence += oligonucleotide[overlap:]

    def get_vertex_id_by_index(self, ind: int) -> int:
        return self.__vertices_ids[ind]

    def get_overlaps(self) -> list[int]:
        return copy.deepcopy(self.__overlaps)

    def get_vertices_count(self) -> int:
        return len(self.__vertices_ids)

    def get_oligonucleotide_length(self) -> int:
        return copy.deepcopy(self.__oligonucleotide_length)

    def get_sequence_length(self) -> int:
        return len(self.__sequence)

    def get_sequence(self) -> str:
        return self.__sequence

    def get_accuracy(self, n_optimal_used_oligonucleotides: int) -> float:
        n_used_oligonucleotides = len(self.__vertices_ids)
        accuracy = n_used_oligonucleotides / n_optimal_used_oligonucleotides
        return accuracy

    def __str__(self) -> str:
        return f"Used oligonucleotides: {len(self.__vertices_ids)}\nGenerated sequence:\n{self.__sequence}"

    def to_csv_str(self, record_name: str, n_optimal_used_oligonucleotides: int) -> str:
        n_used_oligonucleotides = len(self.__vertices_ids)
        accuracy = n_used_oligonucleotides / n_optimal_used_oligonucleotides
        sequence_length = len(self.__sequence)
        mean_overlap = np.mean(self.__overlaps)
        return f"{record_name},{n_used_oligonucleotides},{accuracy},{sequence_length},{mean_overlap}"
