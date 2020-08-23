
import pymol
from pymol import cmd, util
from glob import glob
import re

#cmd.load('nucleophilic.pse')

nloop = 5
seconds = 4
fps = 30
degrees = 30
#cmd.mset(1, )
for loop in range(nloop):
    cmd.movie.add_nutate(seconds,degrees)
nframes = seconds * fps
cmd.mview('store', 1, object='2020-08-20-benzotriazoles-dockscores-x10876', state=1)
cmd.mview('store', nloop*seconds*fps, object='2020-08-20-benzotriazoles-dockscores-x10876', state=nloop*seconds*fps)

cmd.viewport(720, 720)
#movie.produce('nucleophilic_displacement.mov', mode='ray')
