
import pymol
from pymol import cmd, util
from glob import glob
import re

#cmd.load('nucleophilic.pse')

cmd.mset()
cmd.rewind()
nloop = 4
seconds = 4
fps = 30
degrees = 30
#cmd.mset(1, )
nframes = seconds * fps
for loop in range(nloop):
    cmd.movie.add_nutate(seconds,degrees)

objname = '2020-08-20-benzotriazoles-dockscores-x10876'
cmd.load(objname + '.sdf')
cmd.hide('sticks', 'hydrogen')
cmd.mview('store', 1, object=objname, state=1)
cmd.mview('store', nloop*seconds*fps, object=objname, state=nloop*seconds*fps)

cmd.viewport(720, 720)
#movie.produce('nucleophilic_displacement.mov', mode='ray')
