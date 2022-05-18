import os


def read_file(filepath: str) -> list[str] or None:
    try:
        with open(f"{os.getcwd()}/{filepath}", "r") as f:
            return f.read().splitlines()
    except FileNotFoundError:
        print("Given file does not exist")
        return None
