
import pymol
from pymol import cmd, util
from glob import glob
import re

#cmd.load('nucleophilic.pse')

seconds = 4
fps = 30
degrees = 30
cmd.mset(1, )
for loop in range(10):
    cmd.movie.add_nutate(seconds,degrees)
nframes = seconds * fps
cmd.mview('store', 1, object='nucleophilic_displacement_enumeration_for_FEP-sorted-x10789', state=1)
cmd.mview('store', 120*10, object='nucleophilic_displacement_enumeration_for_FEP-sorted-x10789', state=120*10)


#movie.produce('nucleophilic_displacement.mov', mode='ray')
