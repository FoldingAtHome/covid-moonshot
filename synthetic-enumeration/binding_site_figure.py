import pymol
from pymol import cmd, util
from glob import glob
import re

cmd.reset()
cmd.delete('all')

path = "/Users/choderaj/github/foldingathome/covid-moonshot/receptors/monomer/"

# Aim 1 : 3-aminopyridines (cyan)
fragments = ['x2646', 'x10201', 'x10387', 'x10789', 'x10019', 'x10236', 'x10237', 'x10959']
for fragment in fragments:
  cmd.load(path + f'Mpro-{fragment}_0A_bound-ligand.mol2', f'aminopyridines-{fragment}-ligand')
  cmd.load(path + f'Mpro-{fragment}_0A_bound-protein.pdb', f'aminopyridines-{fragment}-protein')
util.cbac(f'aminopyridines-*')

# Aim 2 : quinolones (magenta)
fragments = ['x2910', 'x3080', 'x3303']
for fragment in fragments:
  cmd.load(path + f'Mpro-{fragment}_0A_bound-ligand.mol2', f'quinolones-{fragment}-ligand')
  cmd.load(path + f'Mpro-{fragment}_0A_bound-protein.pdb', f'quinolones-{fragment}-protein')
util.cbam(f'quinolones-*')

# Aim 3 : benzotriazoles (green)
fragments = ['x10820', 'x10871', 'x10876']
for fragment in fragments:
  cmd.load(path + f'Mpro-{fragment}_0A_bound-ligand.mol2', f'benzotriazoles-{fragment}-ligand')
  cmd.load(path + f'Mpro-{fragment}_0A_bound-protein.pdb', f'benzotriazoles-{fragment}-protein')
util.cbag(f'benzotriazoles-*')

# N3 ligand (white)
cmd.load('6lu7.pdb')
cmd.align('(6lu7 and chain A)', 'aminopyridines-x10237-protein')
cmd.select('N3-ligand', '6lu7 and chain C')
cmd.select('N3-protein', '6lu7 and (chain A or chain B)')
util.cbaw('N3-ligand')

# All fragments
show_fragments = False
if show_fragments:
    files = glob(f'{path}/Mpro-x*.mol2')
    for filename in files:
        print(filename)
        match = re.search('Mpro-(?P<fragment>x\d+)_', filename)
        fragment = match.group('fragment')
        cmd.load(filename, f'fragment-{fragment}')
    cmd.select('fragments', 'fragments-*')

# remove waters
cmd.remove('resn HOH')
cmd.deselect()

# Show molecular representation
cmd.hide('all')
cmd.dss('6lu7')
cmd.bg_color('white')
util.cbaw('*-protein')

cmd.show('sticks', f'*-protein and not hydrogen')
cmd.show('surface', f'*-protein')
cmd.disable('*-protein')
cmd.enable('N3-protein')

#cmd.show('sticks', f'({molecule}-protein and not hydrogen) within 7 of N3-ligand')
cmd.set('surface_color', 'white')

molecule = 'N3'
cmd.show('sticks', f'{molecule}-ligand and not hydrogen')
#cmd.show('cartoon', f'{molecule}-protein')
if show_fragments:
    cmd.show('lines', 'fragments and not hydrogen')
cmd.set('surface_mode', 3)
cmd.set('transparency', 0.2)
cmd.set('transparency', 0.8, 'resi 145 and resn CYS')
cmd.set('transparency', 0.8, 'resi 41 and resn HIS')
# Show CYS145 and HIS41
#cmd.show('sticks', '(resi 145 and resn CYS) and not hydrogen')
#cmd.show('sticks', '(resi 41 and resn HIS) and not hydrogen')

# DEBUG
show_conformers = False
if show_conformers:
    cmd.load('/Users/choderaj/github/foldingathome/covid-moonshot/perses-figure/test.mol2', 'conformers')
    cmd.load('/Users/choderaj/github/foldingathome/covid-moonshot/perses-figure/boronic_ester_enumeration_for_chodera_lab_FEP-permuted-dockscores-x10789.sdf', 'conformer')
    cmd.set('all_states', 1)

# DEBUG
#cmd.load('nucleophilic_displacement_enumeration_for_FEP-sorted-x10789.mol2')
#cmd.load('sprint1-winners.mol2')
cmd.viewport(720,720)

cmd.orient('N3-ligand')
cmd.turn('z', -90) # TODO: Shift camera instead
cmd.turn('y', +20) # TODO: Shift camera instead
cmd.zoom('N3-ligand')
#cmd.clip('far', -10)
cmd.move('y', +1)

cmd.show('sticks', 'not hydrogen and *-ligand')

 # Hide all fragments
cmd.deselect()

cmd.disable('aminopyridines-*')
cmd.disable('quinolones-*')
cmd.disable('benzotriazoles-*')
cmd.disable('fragment-*')

if show_conformers:
    cmd.hide('sticks', 'hydrogen')
    cmd.hide('sticks', 'conformers')
    cmd.show('lines', 'conformers and not hydrogen')
    util.cbag('conformer*')
    cmd.disable('6lu7')
    cmd.enable('aminopyridines-x10789-*')
    util.cbaw('aminopyridines-x10789-*')


cmd.hide('sticks', 'hydrogen')

#P1'
#show lines, (byres fragment-* within 2.5 of resi 25) and not hydrogen
# P1
#show lines, (byres fragment-* within 3 of resi 172) and not hydrogen
# P4
#show lines, (byres fragment-* within 2.5 of resi 192) and not hydrogen

pymol.finish_launching()
