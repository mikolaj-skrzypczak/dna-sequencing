class Vertex:
    def __init__(self, _id: int, oligonucleotide: str):
        self.id = _id
        self.oligonucleotide = oligonucleotide

    def __str__(self) -> str:
        return f"[{self.id}] {self.oligonucleotide}"
