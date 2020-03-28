#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Script to Split initial SDF into individual ligands, and create protonated SDF's and PDB's for further parameterization for Coronavirus FAH stuff
# Written by Dylan Novack 03/12/20 
# Contact:tuj66686@temple.edu for Q's

# Imports
import os,glob
import openbabel
from openbabel import pybel
import argparse
import warnings
warnings.filterwarnings("ignore")


# In[ ]:


# Arguments-Raw SDF File for input
parser = argparse.ArgumentParser(description="COVID-19 Ligand SDF Parser: Split initial combo ligand SDF into individual ligands and create protonated SDF's and PDB's for further parameterization for Coronavirus FAH stuff")
parser.add_argument('-s','--sdf', action='store', help=' Raw SDF file to be converted into protonated SDF and PDB')
arg=parser.parse_args()
print(f'{arg.sdf} will be Converted into simulation ready SDF and PDB' )


# In[ ]:


# Define Path and previous parameterizations
path = os.getcwd()
print(f'Current Directory is {path}')
LIG=[]
for directory in glob.glob(f'{path}/LIG*'):
    for filename in os.listdir(directory):
                if filename.endswith('h.sdf'):
                   LIG.append(os.path.join(directory,filename))
prior= len(LIG)
print(f'Prior Ligand Count is {prior}')


# In[ ]:


### Create New Directories and Unprotonated SDF's for new ligands
list=[]
for num, mol in enumerate(pybel.readfile("sdf", f'{path}/{arg.sdf}'),prior+1):
    print(f'Making Directory for LIG{num}')
    os.makedirs(f'{path}/LIG{num}', exist_ok=True)
    dir = f'{path}/LIG{num}'
    list.append(dir)
    print(f'Writing Unprotonated SDF file to {dir}')
    mol.write("sdf", f'{dir}/LIG{num}.sdf', overwrite=True)


# In[ ]:


# Make List of SDF Files
LIG=[]
for dir in list:
    for filename in os.listdir(dir):
        if filename.endswith('.sdf'):
           LIG.append(os.path.join(dir,filename))


# In[ ]:


# Create Protonated SDF Files
#for sdf in LIG:
#    obConversion = openbabel.OBConversion()
#    obConversion.SetInAndOutFormats("sdf", "sdf")
#    mol = openbabel.OBMol()
#    obConversion.ReadFile(mol, sdf)
#    base=os.path.basename(sdf)
#    name=os.path.splitext(base)[0]
#    dirname=os.path.dirname(sdf)
#    mol.AddHydrogens()
#    outMDL = obConversion.WriteFile(mol,f'{dirname}/{name}_h.sdf')
#    print(f'{name} Protonated and SDF Written to {dirname}/{name}_h.sdf ')
#
#
## In[ ]:
#
#
## Make List of Protonated SDF Files
#LIG=[]
#for dir in list:
#    for filename in os.listdir(dir):
#        if filename.endswith('h.sdf'):
#           LIG.append(os.path.join(dir,filename))
#
#
## In[ ]:
#
#
# Make Protonated PDB Files Ready to be Fed into Parameterization Pipeline
for sdf in LIG:
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("sdf", "pdb")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, sdf)
    base=os.path.basename(sdf)
    name=os.path.splitext(base)[0]
    dirname=os.path.dirname(sdf)
    outMDL = obConversion.WriteFile(mol,f'{dirname}/{name}.pdb')
    print(f'{name} PDB Written to {dirname}/{name}.pdb' )
print(f'Ready to Model.' )

