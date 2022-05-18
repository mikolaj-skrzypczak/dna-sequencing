from dna_sequencing.graph.graph import Graph
from dna_sequencing.solution.solution_container import SolutionContainer


class SimulatedAnnealingSolver:
    __graph: Graph
    __solution: SolutionContainer
    __steps: int
    __initial_temperature: int
    __n: int

    def __init__(self, graph: Graph, solution: SolutionContainer, steps: int, initial_temperature: int, n: int) -> None:
        self.__graph = graph
        self.__solution = solution
        self.__steps = steps
        self.__initial_temperature = initial_temperature
        self.__n = n

    # TODO
