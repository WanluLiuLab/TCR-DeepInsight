import os
from ._dataloader import *

PATH = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(PATH, 'datasets'), exist_ok=True)
os.makedirs(os.path.join(PATH, 'pretrained_weights'), exist_ok=True)