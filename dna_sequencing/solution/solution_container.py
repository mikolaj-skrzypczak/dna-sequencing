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

    def get_minimized_function_value(self) -> np.ndarray:  # actually np.intc
        actual_overlaps = np.asarray(self.__overlaps, dtype=np.intc)
        ideal_overlaps = np.full(
            shape=actual_overlaps.shape,
            fill_value=self.__oligonucleotide_length,
            dtype=actual_overlaps.dtype
        )
        return np.sum(ideal_overlaps - actual_overlaps)

    def init_ids(self, _id: int) -> None:
        self.__vertices_ids = [_id]

    def init_sequence(self, oligonucleotide: str):
        self.__sequence = oligonucleotide

    def add_vertex_to_solution(
            self, _id: int, overlap: int, oligonucleotide: str
    ):
        self.__vertices_ids.append(_id)
        self.__overlaps.append(overlap)
        self.__sequence += oligonucleotide[overlap:]

    def get_overlaps(self) -> list[int]:
        return copy.deepcopy(self.__overlaps)

    def get_oligonucleotide_length(self) -> int:
        return copy.deepcopy(self.__oligonucleotide_length)

    def delete_last(self) -> None:
        self.__sequence = self.__sequence[:-(self.__oligonucleotide_length - self.__overlaps[-1])]
        self.__vertices_ids.pop()
        self.__overlaps.pop()

    def __str__(self) -> str:
        return f"""Solution:\nUsed oligonucleotides: {len(self.__vertices_ids)}\nGenerated sequence:{self.__sequence}"""
