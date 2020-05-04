import os
import sys
sys.path.insert(0, os.path.realpath('__file__'))

import pyidi

filename_2 = 'data/data_synthetic.cih'
video_2 = pyidi.pyIDI(filename_2)

Points = pyidi.selection.Select(video_2)
print(Points.points)
print(Points.polygon)