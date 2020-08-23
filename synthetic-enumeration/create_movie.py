
import pymol
from pymol import cmd, util
from glob import glob
import re

#cmd.load('nucleophilic.pse')

nloop = 4
seconds = 4
fps = 30
degrees = 30
#cmd.mset(1, )
nframes = seconds * fps
for loop in range(nloop):
    cmd.movie.add_nutate(seconds,degrees)

for frame in range(nframes):
    print(frame)
    cmd.mview('store', frame+1, object='2020-08-20-benzotriazoles-dockscores-x10876', state=frame+1)
#cmd.mview('store', nloop*seconds*fps, object='2020-08-20-benzotriazoles-dockscores-x10876', state=nloop*seconds*fps)

cmd.viewport(720, 720)
#movie.produce('nucleophilic_displacement.mov', mode='ray')
