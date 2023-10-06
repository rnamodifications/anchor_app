#!/usr/bin/python3

import sys
assert sys.version_info.major >= 3, 'Error: this package requires Python 3.'

import numpy as np
import os
from time import time
from statsmodels.nonparametric.smoothers_lowess import lowess
from multiprocessing import Manager
import multiprocessing as mp
from copy import copy
import re
import csv

from .classes import isnum, compoundset, ppm2dm, adductset, baseset
from treelib import Node, Tree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime

#linkedmass = []
# Determine processor/thread count to choose an optimal number of spawned processes
try:
    MAXTHREADS_ = mp.cpu_count()
except NotImplementedError:
    try:
        MAXTHREADS_ = os.cpu_count()
    except:
        MAXTHREADS_ = 4

# Global options
PARAMS_SET_ = False
PARAMS_ = {

    # ADDUCT CLUSTERING PARAMS
    'ADDUCTDB_CSV': None,
    'CLUSTER_PPM': None,
    'CLUSTER_RTWIDTH_FACTOR': None,
    'CLUSTER_MAX_DEPTH': None,

    # BASE MASS DB
    'BASEDB_CSV': None,

    # STATIC FILTERING PARAMS
    'FILTER_MIN_MASS': None,
    'FILTER_MIN_QUALITY': None,
    'FILTER_MIN_SCORE': None,
    'FILTER_MIN_INTENSITY': None,
    'FILTER_MIN_INTENSITYPCT': None,
    'FILTER_MIN_IONS': None,
    'FILTER_MIN_ZCOUNT': None,
    'FILTER_MIN_RT': None,
    'FILTER_MAX_RT': None,
    'FILTER_ORPHAN_WEIGHT': None,

    # LOCAL FITTING PARAMS
    'LOWESS_FRACTION': None,

    # RANDOM WALK PARAMS
    'WALK_DRT_MAX': None,
    'WALK_TRIALS': None,
    'WALK_STEPS_MIN_DRAFT': None,
    'WALK_STEPS_MIN_FINAL': None,
    'WALK_PPM': None,
    'WALK_INTERNAL_NOISE_VARIANCE': None,
    'WALK_STARTING_NOISE_VARIANCE': None,
    'WALK_TOPFRACTION': None,
    'WALK_BIDIR': None,

    # PLOTTING PARAMS
    'PLOT': None,
    'PLOT_SAVE': None,
    'PLOT_SAVE_FMT': None,
    'PLOT_LABELS': None,
    'PLOT_MASS_MIN': None,
    'PLOT_MASS_MAX': None,
    'PLOT_RT_MIN': None,
    'PLOT_RT_MAX': None,
    'PLOT_ALPHA_HIGH': None,
    'PLOT_ALPHA_LOW': None,
    'PLOT_ALL_WALKS': None,
    'PLOT_MARKERSIZE': None,
    'PLOT_USETEX': None
}

def read_params(cfgfilename):
    """
    Import parameter csv file and load parameters in gerbil scope.

    Takes only a filename that must be a csv containing the names of all of the __PARAMS__ global dict keys.
    """
    global PARAMS_SET_
    global PARAMS_

    try:
        with open(cfgfilename) as csvfile:
            csvfile.seek(0)
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                if not len(row):
                    continue
                if "#" in row[0]:
                    continue
                for pkey in PARAMS_:
                    if row[0] == pkey:
                        val = row[1].split('#')[0].strip()
                        if len(val) <= 0 or val is None or val == 'None' or val == 'none' or val == '':
                            PARAMS_[pkey] = None
                        elif val == 'True' or val == 'true' or val is True:
                            PARAMS_[pkey] = True
                        elif val == 'False' or val == 'false' or val is False:
                            PARAMS_[pkey] = False
                        elif isnum(val):
                            PARAMS_[pkey] = float(val)
                        else:
                            PARAMS_[pkey] = val
    except:
        PARAMS_SET_ = False
        print('Error importing parameters.  Fail at everything.')
        return False

    assert os.path.isfile(PARAMS_['ADDUCTDB_CSV']), 'Error: adducts database csv file not found.'
    assert os.path.isfile(PARAMS_['BASEDB_CSV']), 'Error: base database csv file not found.'

    PARAMS_SET_ = True
    return True


def define_startingpos(compounds):
    """
    Define the likely starting point as the compound with maximum integrated intensity.
    (This function usually finds the full length compound, assuming that it has maximal abundance in the dataset. This
     is useful primarily for trying to guess the read orientation after walk generation.)

    :compounds: list of compound dicts containing keys 'Mass', 'RT', 'Vol', 'Cpd', and 'Width'
    """

    # TODO: FOR THE TIME BEING, GUTTED ALL FEATURES THAT PERTAIN TO POOR INITIAL SAMPLE QUALITY

    maxintensity = np.max([c['Vol'] for c in compounds])
    pickedcompound = [c for c in compounds if c['Vol'] == maxintensity][0]

    # What if the sample is overdigested, and the starting material cannot be assumed to be one with highest intensity?
    #    Check to see if >20% of compound intensity is larger than the picked one.
    largercompounds = [c for c in compounds if c['Mass'] > pickedcompound['Mass']]
    largervol = np.sum([c['Vol'] for c in largercompounds])
    totalvol = np.sum([c['Vol'] for c in compounds])
    if largervol/totalvol > 0.2:
        # choose starting point as the one with maximum mass, to avoid truncating the dataset incorrectly on the mass
        # axis
        maxmass = np.max([c['Mass'] for c in compounds])
        pickedcompound = [c for c in compounds if c['Mass'] == maxmass][0]

    return pickedcompound


def findbymassrt(compounds, m, t, ppm, tstringency):
    """
    Find compounds near mass m and RT t within a mass ppm and RT stringency window

    :param compounds: list of compound dicts containing keys 'Mass', 'RT', and 'Width'
    :param m: float mass
    :param t: float retention time
    :return: list of compound dict matches
    """
    mstringency = ppm2dm(m, ppm)
    mcandidates = [c for c in compounds if ((c["Mass"] >= m - mstringency)
                                            & (c["Mass"] <= m + mstringency))]
    matches = [c for c in mcandidates if ((c["RT"] >= t - c["Width"] * tstringency)
                                          & (c["RT"] <= t + c["Width"] * tstringency))]
    return matches


def randstep(cpd, compounds, bases, dmmin, dmmax, dt, tdir, ppm, variance, manager, start=None, internal=None, end=None, ppm_multi=1):
    candidates = [c for c in compounds 
            if ((c['Mass'] >= cpd['Mass'] - dmmax)
               & (c['Mass'] <= cpd['Mass'] - dmmin)
               & (c['Cpd'] != cpd['Cpd']))]

    #print("===>{}".format(cpd['Cpd']))
    #print("--->{}".format(candidates))
    my_list = list()
    for i in range(len(candidates)):
        basecall = bases.findidbymass(cpd['Mass'] - candidates[i]['Mass'],
                      ppm_multi*ppm2dm(candidates[i]['Mass'], ppm),
                      #1.05,
                      start=start,
                      end=end,
                      internal=internal)
        if basecall:
            #my_list.append(((cpd['Mass'], cpd['RT']), (candidates[i]['Mass'], candidates[i]['RT'])))
            my_list.append((candidates[i]['Mass'], candidates[i]['RT']))

    if my_list:
        link = [(cpd['Mass'], cpd['RT']), my_list]
        if link not in manager:
            manager.append(link)

def biased_walk(cpd, compounds, dmmin, dmmax, dt, ppm, maxtrials, minsteps, maxsteps, startingpos, bases,
                startingvar, internalvar, manager, bidirectional=False, ppm_multi=1):
    trials = []
    pos = []
    pos.append(copy(cpd))
    pos[0]['WalkScore'] = pos[0]['Vol']  # assign initial score for this position as its intensity
    pos[0]['Call'] = 0
    tdir = 0
    nextpos = randstep(pos[-1], compounds, bases, dmmin, dmmax, dt, tdir, ppm,
                       variance=startingvar,
                       manager = manager,
                       start=True, ppm_multi=ppm_multi)
    trials.append(nextpos)
    return trials


def generate_trajectories(compounds, dt, ppm, maxtrials, minsteps, maxsteps, startingpos, bases,
                          startingvar, internalvar, manager, bidirectional=False, ppm_multi=1):
    dmmin = float(np.min(bases.mass)-1.0)
    dmmax = float(np.max(bases.mass)+1.0)
    print('Using %d cores.' % MAXTHREADS_)
    with mp.Pool(MAXTHREADS_) as p:
        trials = p.starmap(biased_walk, [(cpd, compounds, dmmin, dmmax, 
            dt, ppm, maxtrials, minsteps, maxsteps,
              startingpos, bases, startingvar, internalvar, manager,
              bidirectional,ppm_multi) for cpd in compounds])


    alltrials = []
    for trial in trials:
        alltrials.extend(trial)

    return alltrials


def nodes_with_tag(tree, leaf):
    nodes_list = list()
    for node in tree.all_nodes():
        if node.tag == leaf:
            nodes_list.append(node)

    return nodes_list

def preprocess_lines(lines):
    print("preprocess_lines {}".format(len(lines)))
    l_lines = []
    m_lines = []
    r_lines = []
    lines_x = [line[0] for line in lines]
    lines_y = [line[1] for line in lines]
    for line in lines:
        is_on_right = False
        is_on_left = False
        for line_y in lines_y:
            is_on_right = line[0] in line_y
            if is_on_right:
                break

        for line_y in line[1]:
            if line_y in lines_x:
                is_on_left = True
                break
        #if is_on_right:
            #r_lines.append(line)
        #'''
        #is_on_left = line[1] in lines_x
        if is_on_right and is_on_left:
            m_lines.append(line)
        elif is_on_right:
            r_lines.append(line)
        elif is_on_left:
            l_lines.append(line)
        #'''

    return l_lines, m_lines, r_lines

def process_lines(l_lines, m_lines, r_lines):
    l_dicts = dict()
    m_dicts = dict()
    r_dicts = dict()
    full_dicts = dict()
    for line in l_lines:
        l_dicts[line[0]] = line[1]
        full_dicts[line[0]] = line[1]
    for line in m_lines:
        m_dicts[line[0]] = line[1]
        full_dicts[line[0]] = line[1]
    for line in r_lines:
        r_dicts[line[0]] = line[1]
        full_dicts[line[0]] = line[1]
    return l_dicts, m_dicts, r_dicts, full_dicts

def process_tree(lines, ori_filename, cdb, bdb, params):
    print("process_tree {} total lines ".format(len(lines)))
    tree = Tree()
    tree.create_node(0, 0) # root
    for idx, line in enumerate(lines):
        nodes_list = nodes_with_tag(tree, line[0])
        if not nodes_list:
            head_node = tree.create_node(line[0], parent=0)
            nodes_list = list()
            nodes_list.append(head_node)
        for node in nodes_list:
            sub_nodes = tree.children(node.identifier)
            sub_nodes_tags = [sub_node.tag for sub_node in sub_nodes]
            for leaf in line[1]:
                if leaf not in sub_nodes_tags:
                    tree.create_node(leaf, parent=node)

    print("process_tree created tree done")
    print("max level of the tree {}".format(tree.depth()))
    max_level = tree.depth()
    ret_list = list()
    seq_collect = list()

    plt.figure(figsize=(9, 7))

    if params.scatter_all:
        points_x = [row['Mass'] for row in cdb.compounds if row['Mass'] < params.mass_max and row['RT'] < params.rt_max]
        points_y = [row['RT'] for row in cdb.compounds if row['Mass'] < params.mass_max and row['RT'] < params.rt_max]
        plt.scatter(points_x, points_y)

    idx = 1
    if len(seq_collect) == 0:
        for id_list in tree.paths_to_leaves():
            seq = [tree.get_node(nid).tag for nid in id_list]
            seq_collect.append(seq[1:])
        
    for p in seq_collect:
        if len(p) > max_level:
            max_level = len(p)

    seq_collect = [p for p in seq_collect if len(p) >= max_level-1]
    for idx, p in enumerate(seq_collect):
        plist = p
        if idx < 3:
            x = [node[0] for node in plist]
            y = [node[1] for node in plist]
            ret_list.append(x)
            lbl = 'Seq'+str(idx)
            plt.plot(x, y, label=lbl, marker='o')
            # add mark for each node
            x_len = len(x)
            mass_diff_list = []
            for mass_idx in range(x_len-1):
                m_diff = x[mass_idx] - x[mass_idx+1]
                mass_diff = "{0:.4f}".format(m_diff)
                mass_diff_list.append(mass_diff)
            mass_diff_list.append('')
            for p_idx, point in enumerate(zip(x, y)):
                if mass_diff_list[p_idx] != '':
                    mass_diff = float(mass_diff_list[p_idx])
                else:
                    mass_diff = 0.0
                stringency = params.ppm_multi*ppm2dm(mass_diff, PARAMS_['WALK_PPM'])
                mass_name = bdb.findnamebymass(mass_diff, stringency)
                plt.annotate(mass_name, point, size=20) 
            idx += 1

    tm = datetime.now()
    img_file = 'static/a%d_%d_%d_%d.png' % (tm.hour, tm.minute, tm.second, tm.microsecond)
    plt.title(ori_filename)
    plt.xlabel('Mass')
    plt.ylabel('RT')
    plt.legend()
    plt.grid(linestyle='--')
    plt.savefig(img_file)
    return (ret_list, img_file)


def process(compounds_csvfile, ori_filename, params):

    assert os.path.isfile(compounds_csvfile), 'Error: compounds dataset csv file not found.'

    assert PARAMS_SET_, 'Error: cannot process dataset without first loading parameters, buddy.  Run read_params().'

    print('Processing compounds file: ' + compounds_csvfile)
    start_time = time()

    print(params)
    ppm_multi = params.ppm_multi


    ## IMPORT SOURCE CSV FILES

    # define and import compoundset from csv file
    # compoundset.compounds is a list of compound dicts
    cdb = compoundset(compounds_csvfile, params)

    # define and import adducts from csv file
    adb = adductset(PARAMS_['ADDUCTDB_CSV'])

    # define and import base mass differences from csv file
    bdb = baseset(PARAMS_['BASEDB_CSV'], params)

    print('Found ' + str(len(cdb.compounds)) + ' compounds.')

    # Sort compoundset by abundance (i.e. intensity volume)
    cdb.sort('Vol')

    ## COMPOUND FILTERING

    # prefilter on quality and time-range before doing anything else (generally unused)
    if PARAMS_['FILTER_MIN_QUALITY'] is not None:
        cdb.filter(quality=(PARAMS_['FILTER_MIN_QUALITY'], 101))
    if PARAMS_['FILTER_MIN_RT'] is not None and PARAMS_['FILTER_MAX_RT']:
        cdb.filter(time=(PARAMS_['FILTER_MIN_RT'], PARAMS_['FILTER_MAX_RT']))

    # Agglomerative adduct clustering
    print('Cluster adducts, stringency: ppm = ' + str(PARAMS_['CLUSTER_PPM']) + ', RT peak width factor = ' +
          str(PARAMS_['CLUSTER_RTWIDTH_FACTOR']))

    print('Adducts clustering leaves ' + str(len(cdb.compounds)) + ' compounds.')

    # Find candidate starting point for full length material
    # (Note: This is a starting point to later attempt to append to a final set of walks.  Walks from biased_walk
    # start from every compound point, not just this one.)
    startingpos = define_startingpos(cdb.compounds)
    print('Chose ' + str(startingpos['Mass']) + ' as candidate full length mass.')

    # This feature down-weights any compound for which no corresponding compound can be found that would sum to the
    # correct starting point mass (not necessary for most datasets, and may incidentally exclude valuable points
    # in the high mass region due to lost fragments in the low mass region).
    if PARAMS_['FILTER_ORPHAN_WEIGHT'] is not None and PARAMS_['FILTER_ORPHAN_WEIGHT'] > 0:
        cdb.weight_orphans(startingpos, PARAMS_['FILTER_ORPHAN_WEIGHT'], PARAMS_['WALK_PPM']*10)

    # Try to predict number of fragments using max intensity fragment and the average native base mass
    fragcount = np.round(((startingpos['Mass'] / 317) * 2 - 1))

    # apply the selected filters
    if PARAMS_['FILTER_MIN_MASS'] is not None:
        cdb.filter(mass=(PARAMS_['FILTER_MIN_MASS'], startingpos['Mass']))

    if PARAMS_['FILTER_MIN_INTENSITY'] is not None and PARAMS_['FILTER_MIN_INTENSITY'] > 0:
        cdb.filter(intensity=(PARAMS_['FILTER_MIN_INTENSITY'], 1.0E10),)

    if PARAMS_['FILTER_MIN_INTENSITYPCT'] is not None and PARAMS_['FILTER_MIN_INTENSITYPCT'] > 0:
        cdb.filter(intensitypct=(PARAMS_['FILTER_MIN_INTENSITYPCT'], 100))

    if PARAMS_['FILTER_MIN_IONS'] is not None and PARAMS_['FILTER_MIN_IONS'] > 0:
        cdb.filter(ions=(PARAMS_['FILTER_MIN_IONS'], 999))

    if PARAMS_['FILTER_MIN_ZCOUNT'] is not None and PARAMS_['FILTER_MIN_ZCOUNT'] > 0:
        cdb.filter(charges=(PARAMS_['FILTER_MIN_ZCOUNT'], 999))

    if PARAMS_['FILTER_MIN_SCORE'] is not None and PARAMS_['FILTER_MIN_SCORE'] > 0:
        cdb.filter(score=(PARAMS_['FILTER_MIN_SCORE'], 999))

    #if PARAMS_['FILTER_EXPFRAGMENT_FACTOR'] is not None and PARAMS_['FILTER_EXPFRAGMENT_FACTOR'] > 0:
    #    loosefragcount = int(fragcount * (1 + PARAMS_['FILTER_EXPFRAGMENT_FACTOR']))
    #    cdb.filter(fragments=loosefragcount)

    # Resort compoundset by abundance (i.e. intensity volume)
    cdb.sort('Vol')


    ## DEFINE MIDLINE FOR CLUSTERING SEQUENCING WALKS

    print('Fitting midline for clustering...')

    m = np.array([d['Mass'] for d in cdb.compounds])
    t = np.array([d['RT'] for d in cdb.compounds])

    # find a midline of the compound data points by LOWESS fitting
    midline = lowess(t, m, frac=PARAMS_['LOWESS_FRACTION'])


    ## GENERATE RANDOM SEQUENCE TRAJECTORIES

    print('Generating walk trajectories...')

    # try to trace contiguous fragment ladder by random walking through the data from a random starting point
    mana = Manager()
    manager = mana.list()
    trials = generate_trajectories(cdb.compounds,
                                   dt=PARAMS_['WALK_DRT_MAX'],
                                   ppm=PARAMS_['WALK_PPM'],
                                   maxtrials=int(PARAMS_['WALK_TRIALS'] * len(cdb.compounds)),
                                   minsteps=int(PARAMS_['WALK_STEPS_MIN_DRAFT']),
                                   maxsteps=int(fragcount),
                                   startingpos=startingpos,
                                   bases=bdb,
                                   manager= manager,
                                   startingvar=PARAMS_['WALK_STARTING_NOISE_VARIANCE'],
                                   internalvar=PARAMS_['WALK_INTERNAL_NOISE_VARIANCE'],
                                   bidirectional=PARAMS_['WALK_BIDIR'], ppm_multi=params.ppm_multi)

    lines = list()
    for line in manager:
        lines.append(line)
    lines.sort(key=lambda x:x[0][0], reverse=True)
    l_lines, m_lines, r_lines = preprocess_lines(lines)
    l_dicts, m_dicts, r_dicts, full_dicts = process_lines(l_lines, m_lines, r_lines)

    print("list left {} mid {} right {}".format(len(l_lines), len(m_lines), len(r_lines)))
    return process_tree(m_lines, ori_filename, cdb, bdb, params)

