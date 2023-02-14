#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 18:21:40 2023

@author: psturm

Obin Sturm (psturm@usc.edu)
"""

import numpy as np
import os
import tensorflow as tf
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Input,Dense
from tensorflow.keras.callbacks import EarlyStopping


"""
Generate weighted, directed incidence (or biadjacency) matrix A
(also called the stoichiometric matrix)
As in Sturm and Wexler (2022): https://doi.org/10.5194/gmd-15-3417-2022
Note that networkx calls this a biadjacency matrix
"""
# Link to small stratospheric mechanism: https://kpp.readthedocs.io/en/latest/getting_started/02_running_kpp_sample_mech.html
# rows of sparse stoichiometric matrix from small_strato_StoichiomSP.f
rows = np.array([2,  2,  3,  2,  3,  2,  3,  1,  3,  1,  2,  
                 1,  3,  3,  4,  5,  2,  4,  5,  2,  4,  5 ])
# columns of sparse stoichiometric matrix from small_strato_StoichiomSP.f
cols = np.array([1,  2,  2,  3,  3,  4,  4,  5,  5,  6,  6,  
                 7,  7,  8,  8,  8,  9,  9,  9, 10, 10, 10 ])
# stoichiometric coefficients of sparse matrix from small_strato_StoichiomSP.f
vals =    np.array([ 2,   -1,    1,    1,   -1,   -1,   -1,    1,
                    -1,   -1,    1,   -1,   -1,   -1,   -1,    1,
                    -1,    1,   -1,    1,    1,   -1 ])

A = np.zeros((5,10))
A[rows-1,cols-1] = vals

n_spc = np.max(rows)
n_rxn = np.max(cols)

# Create edgelist, useful for graph analysis and visualization
elist = [(rows[i]-1, n_spc+cols[i]-1, -vals[i]) for i in range(0,len(rows))]
for i in range(0,len(rows)):
    if elist[i][2]<0:
        elist[i] = (elist[i][1],elist[i][0],-elist[i][2])


"""
Some tensorflow tricks: 
    seeding 
    setting early stopping
"""
seed_val = 1234
os.environ['PYTHONHASHSEED'] = str(seed_val)
np.random.seed(seed_val)
tf.random.set_seed(seed_val) 

# Stop early if a skill score does not improve for a set of validation data
es = EarlyStopping(monitor='val_loss', min_delta=1e-5 , verbose=1, patience=10)


"""
Neural networks:
    FluxNN: 
        trained to predict reaction fluxes S (given representative data for S)
        fluxes in this case are atom fluxes, aka extents of reaction
        tendencies can be found via ∆C = AS (Sturm and Wexler, 2020)
        if A is stoichiometrically balanced, atoms will be conserved
    GraphNN: 
        predicts tendencies ∆C
        the last layer weights are A
        stoichiometry is built into the neural network
        this is useful if S is unknown (Sturm and Wexler, 2022)
        open question: can the penultimate layer in a trained model resemble S?
    
"""

FluxNN = Sequential()
FluxNN.add(Input(shape=(n_spc,)))  # Note: might want to have more input than just species concentrations
FluxNN.add(Dense(20, activation='relu'))
FluxNN.add(Dense(12, activation ='relu'))
FluxNN.add(Dense(n_rxn, activation='relu')) # In this example, FluxNN predicts reaction fluxes of size n_rxn
FluxNN.compile(loss = 'mse', optimizer='adam', callbacks=[es],
                  metrics=['mean_absolute_error', tf.keras.metrics.RootMeanSquaredError()])
# if FluxNN predicts fluxes S, then tendencies can be found via ∆C = AS


GraphNN = Sequential()
GraphNN.add(Input(shape=(n_spc,)))  # Note: might want to have more input than just species concentrations
GraphNN.add(Dense(20, activation='relu'))
GraphNN.add(Dense(12, activation ='relu'))
GraphNN.add(Dense(n_spc,trainable = False))
GraphNN.compile(loss = 'mse', optimizer='adam', callbacks=[es],
                  metrics=['mean_absolute_error', tf.keras.metrics.RootMeanSquaredError()])
# These steps built the graph relational structure of the chemical mechanism into the neural network
weights = GraphNN.layers[2].get_weights()[0]
biases = GraphNN.layers[2].get_weights()[1]
GraphNN.layers[2].set_weights([(A).transpose(),biases])
GraphNN.summary()

