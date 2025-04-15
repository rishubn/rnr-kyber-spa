from __future__ import annotations
from typing import Optional
from dataclasses import dataclass
import numpy.typing as npt
import pickle
from intt import INTTImpl


@dataclass
class Artifact:
    intt_values: npt.NDArray
    sigma: float
    num_expr: int
    iterations: int
    non_zero_sparse: int
    word_size: int
    priors: npt.NDArray
    attack_results: npt.NDArray
    intt_impl: INTTImpl

    @classmethod
    def save_artifact(cls, file_path: str, artifact: Artifact) -> None:
        with open(file_path, "wb") as f:
            pickle.dump(artifact, f)

    @classmethod
    def load_artifact(cls, file_path: Optional[str]) -> Optional[Artifact]:
        if file_path is not None:
            with open(file_path, "rb") as f:
                return pickle.load(f)
