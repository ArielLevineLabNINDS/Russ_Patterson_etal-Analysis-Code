#!/usr/bin/env python3

import logging

# Configure logging
logging.basicConfig(
                    level=logging.INFO, 
                    filename="SVM.log", 
                    filemode="a+",
                    format='%(asctime)s - %(message)s'
                    )

import pickle as pkl

import numpy as np
import pandas as pd
import scipy.sparse as sparse
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, MaxAbsScaler
from sklearn.svm import LinearSVC


# Data not provided in repository
# This loads the sparse counts matrix
X = sparse.load_npz("all_cells_count_matrix.npz")
logging.info(
             "Data loaded\n" +
             f"\tdtype: {X.dtype}\n" +
             f"\tshape: {X.shape}\n" +
             f"\tnnz: {X.nnz}"
             )

# Data not provided in repository
# This loads the cluster labels
with open("cluster_labels.pkl", "rb") as file:
    y = pkl.load(file)
    y = y.id4.to_numpy()

logging.info(
             "Clusters loaded\n" +
             f"\thead: {y[:5]}\n" +
             f"\tdims: {y.shape}" +
             f"\tdtype: {y.dtype}"
             )

# Encode labels
le = LabelEncoder()
y = le.fit_transform(y)
y = y.astype(int)

with open('results/label_encoder.pkl', 'wb') as file:
        pkl.dump(le, file)

logging.info('Labels encoded and saved')
logging.info(f'{y[:5]}, {y.dtype}')

# Split data
skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
train_idx, test_idx = next(skf.split(X, y))

X_train, X_test = X[train_idx], X[test_idx]
y_train, y_test = y[train_idx], y[test_idx]
logging.info(f'Data split. Train: {len(train_idx)}, Test: {len(val_idx)}')

# Create model
scl = MaxAbsScaler(copy=False)
clf = LinearSVC(dual=False, random_state=42, max_iter=10000)
pl = Pipeline([
               ('Scale', scl),
               ('Classifier', clf)
               ])
logging.info(
             "Pipeline created\n" +
             f"\t{pl}"
             )

# Train and Evaluate
logging.info("Training begins")
pl.fit(X_train, y_train)
logging.info(f"Training completed")

# Save results
with open("predicted.pkl", "wb") as file:
        pkl.dump(pl.predict(X_test), file)
        logging.info("Predictions saved!")

with open("actual.pkl", "wb") as file:
        pkl.dump(y_test, file)
        logging.info("Actual saved!")

with open("svm_pipeline.pkl", "wb") as file:
        pkl.dump(pl, file)
        logging.info("Pipeline saved!")

logging.info("Training complete!")
