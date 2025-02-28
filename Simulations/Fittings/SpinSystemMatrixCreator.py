"""
Module Documentation here

Home for a handful of useful spin system matrix related functions.
Spin system matrices are awkward to handle and can be confusing for users.
These functions make working with spin system matrices easier.
"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
_copyright_ = ""
_credits_ = ""
_licence_ = ("")
_reference_ = ("")
#=========================================================================================
# Last code modification:
#=========================================================================================
_modifiedBy_ = "$modifiedBy: Morgan Hayward $"
_dateModified_ = "$dateModified: 2023-11-24 11:42:59 +0000 (Fri, November, 24, 2023) $"
_version_ = "$Revision: 3.2.0 $"
#=========================================================================================
# Created:
#=========================================================================================
_author_ = "$Author: Morgan Hayward $"
_date_ = "$Date: 2023-11-15 10:46:00 +0000 (Wed, November, 15, 2023) $"
#=========================================================================================
# Start of code
#=========================================================================================


from nmrsim import qm
import numpy as np
import pandas as pd
import os
import sys


def createSpinSystemMatrix(splitting_data, multiplet_data, frequency):
    shifts = [float(shift) * frequency for shift in multiplet_data.center]
    dimension = len(multiplet_data.index)
    ss_matrix = np.zeros((dimension, dimension))
    for cp_index, line in splitting_data.iterrows():
        j = float(line.j)
        x = int(line.multiplet_id1.split('.')[-1]) - 1
        y = int(line.multiplet_id2.split('.')[-1]) - 1
        # x = multiplet_data.index[multiplet_data['multiplet_id'] == line.multiplet_id1].tolist()[0]
        # y = multiplet_data.index[multiplet_data['multiplet_id'] == line.multiplet_id2].tolist()[0]
        ss_matrix[x, y] = j
    ss_matrix = ss_matrix + ss_matrix.T
    # if len(ss_matrix) > 11:
    #     print('spin system matrix is too large')
    #     return None
    for index, shift in enumerate(shifts):
        ss_matrix[index][index] = round(shift / frequency, 4)
    return ss_matrix


def peakListFromSpinSystemMatrix(spinSystemMatrix, frequency, width):
    breaks = [index + 1 for index in range(len(spinSystemMatrix) - 1) if not np.any(spinSystemMatrix[:index + 1, index + 1:] != 0)]
    breaks = [0] + breaks + [len(spinSystemMatrix)]
    matrices = [spinSystemMatrix[breakpoint:breaks[index + 1], breakpoint:breaks[index + 1]] for index, breakpoint in enumerate(breaks[:-1])]
    # Suppress the printing from qm
    output = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    peaklist = [item for sublist in [qm.qm_spinsystem([ssm[index][index] * frequency for index in range(len(ssm))], ssm) for ssm in matrices] for item in sublist]
    sys.stdout = output
    peaklist_out = [(peak[0] / frequency, peak[1], width, None) for peak in peaklist]
    df_out = pd.DataFrame(peaklist_out, columns=['chemical_shift', 'height', 'width', 'multiplet_id'])
    return df_out


def fullSpinSystemMatrixFromSimplified(couplingGrid, protonCounts, shifts):
    realSsm = []
    for i, row in enumerate(couplingGrid):
        for j in range(protonCounts[i]):
            newrow = []
            for k, value in enumerate(row):
                newrow += [value] * protonCounts[k]
            realSsm.append(newrow)
    ss_matrix = np.array(realSsm, dtype='float')
    ss_matrix = ss_matrix + ss_matrix.T
    RealShifts = []
    for i, shift in enumerate(shifts):
        RealShifts += [shift] * protonCounts[i]
    for i, shift in enumerate(RealShifts):
        ss_matrix[i][i] = shift
    return ss_matrix

