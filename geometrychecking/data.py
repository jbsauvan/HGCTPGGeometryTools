#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from enum import Enum
import uproot
import awkward as ak
import click

errortypes = Enum('Errors',
                  ['CellValidity',
                   'MissingCellInTC',
                   'InvalidCellInTC',
                   'MissingTCInModule',
                   'InvalidTCInModule',
                   'ConnectedModuleWithoutLpgbt',
                   'ModuleSplitInStage1',
                   'MissingModuleInStage1FPGA',
                   'InvalidModuleInStage1FPGA',
                   'MissingStage1InStage2FPGA',
                   'InvalidStage1InStage2FPGA'],
                  start=0
                  )

subdets = {'cells':[8,9,10],
           'triggercells':[1,2,10],
           'modules':[1,2,3]
           }

def check_error(error, errortype):
    return (error & (0x1 << errortype.value))!=0

def select_error(array, errortype):
    return array[(array.errorbits & (0x1<< errortype.value))>0]

def read_geometry(file_name, tree_names):
    file = uproot.open(file_name)
    trees = {}
    for name,tree_name in tree_names.items():
        trees[name] = file[tree_name]
    return trees

def to_awkward(trees, branches):
    arrays = {}
    for name,tree in trees.items():
        arrays[name] = tree.arrays(branches[name], library='ak')
    return arrays

def extract_layer(arrays, subdet, layer):
    layer_arrays = {}

    for name,array in arrays.items():
        if subdet=='cee':
            layer_arrays[name] = array[(array.layer==layer) & (array.subdet==subdets[name][0]) & (array.zside==1)]
        elif subdet=='ceh':
            layer_arrays[name] = array[(array.layer==layer) & ((array.subdet==subdets[name][1]) | (array.subdet==subdets[name][2])) & (array.zside==1)]
    return layer_arrays

def group_by(array, group):
    array_by = array[ak.argsort(array[group])]
    lengths = ak.run_lengths(array_by[group])
    array_by = ak.unflatten(
        array_by, 
        lengths
    )
    return array_by


def to_parquet(arrays, cells_by_module, cells_by_stc, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for name, array in arrays.items():
        ak.to_parquet(array, f"{output_dir}/{name}.parquet")
    ak.to_parquet(cells_by_module, f"{output_dir}/cells_by_module.parquet")
    ak.to_parquet(cells_by_stc, f"{output_dir}/cells_by_stc.parquet")



@click.command()
@click.option('--geometry', 
              'geometry_name',
              required=True,
              help='Trigger geometry name (name of the geometry ROOT file without the .root extension)'
              )
@click.option('--subdet',
              'subdet',
              required=True,
              type=click.Choice(['cee', 'ceh']),
              help='Subdetector: cee or ceh'
              )
def cli(geometry_name,subdet):
    click.echo('> Reading geometry ROOT tree')
    trees  = read_geometry(f"resources/{geometry_name}.root", 
                           {'cells':'L1THGCaltriggergeomtester/HGCalTriggerGeoTesterCells',
                            'triggercells':'L1THGCaltriggergeomtester/HGCalTriggerGeoTesterTriggerCells',
                            'modules':'L1THGCaltriggergeomtester/HGCalTriggerGeoTesterModules'
                            }
                           )
    branches = {
        'cells': ['id', 'stcid', 'modid', 'zside', 'subdet', 'layer', 'x', 'y', 'corners_x', 'corners_y', 'corners_z','errorbits', 'tcid'],
        'triggercells': ['id', 'modid', 'zside', 'subdet', 'layer', 'x', 'y','errorbits','cells'],
        'modules': ['id', 'cells', 'triggercells', 'zside', 'subdet', 'layer', 'errorbits', 'disconnected']
    }
    arrays = to_awkward(trees, branches)
    if subdet=='cee': maxlayer=26
    elif subdet=='ceh': maxlayer=21
    layers = range(1,maxlayer+1)
    click.echo(f'> Extracting {subdet} layers data and grouping cells by modules')
    with click.progressbar(layers) as bar:
        for layer in bar:
            layer_arrays = extract_layer(arrays, subdet=subdet, layer=layer)
            cells_by_module = group_by(layer_arrays['cells'], group='modid')
            cells_by_stc = group_by(layer_arrays['cells'], group='stcid')
            to_parquet(layer_arrays, cells_by_module, cells_by_stc, f'resources/{geometry_name}/{subdet}_layer{layer}/')

