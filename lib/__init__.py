from ga import GPolySolver 
from raxmlib import RAxMLModel, LklModel
from TreeLib import TreeClass

VERSION = "0.0.1rc"
DESC = "gaperm find the best tree (leaf labeling) from a list of unlabeled trees obtained using polytomysolver. The algorithm proceed by using a genetic algorithm to find the best candidates."

__all__ = ["TreeClass", "GPolySolver", "RAxMLModel", "LklModel", "VERSION", "DESC"]