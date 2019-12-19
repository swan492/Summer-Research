"""
The :mod:`skmultiflow.drift_detection` module includes methods for Concept Drift Detection.
"""

from .adwin import ADWIN
from .ddm import DDM
from .eddm import EDDM
from .page_hinkley import PageHinkley
from .seed import SEED
from .angle import ANGLE

__all__ = ["ADWIN", "DDM", "EDDM", "PageHinkley", "SEED", "ANGLE"]
