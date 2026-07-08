#!/usr/bin/env python3
"""
Quick test of inference notebook functionality
"""
import time
import numpy as np
import pandas as pd
import pickle
from keras.models import load_model
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("INFERENCE NOTEBOOK TEST")
print("="*60)

# Step 1: Load embeddings
print("\n[1/5] Loading embeddings...")
t0 = time.time()
with open('data/embeddings_node2vec_MechDist_dimensions_128_WalkLength_10_NumWalks_100.pickle', 'rb') as f:
    MechDist_embeddings = pickle.load(f)
config = MechDist_embeddings['config']
diseases = [k for k in MechDist_embeddings.keys() if k != 'config']
print(f"      ✅ Loaded {len(diseases)} diseases in {time.time()-t0:.2f}s")

# Step 2: Load model
print("\n[2/5] Loading model...")
t0 = time.time()
model = load_model('data/UPNA-DDAL_epochs_2__embeddings_node2vec_MechDist_dimensions_128_WalkLength_10_NumWalks_100.h5')
print(f"      ✅ Model loaded in {time.time()-t0:.2f}s")

# Step 3: Load disease names
print("\n[3/5] Loading disease names...")
t0 = time.time()
df = pd.read_csv('data/disease_id_names.csv', low_memory=False)
names = df.set_index('disease_id')['name'].to_dict()
print(f"      ✅ Loaded {len(names)} names in {time.time()-t0:.2f}s")

# Step 4: First prediction (TensorFlow init - slow)
print("\n[4/5] Running first prediction (TensorFlow init - may be slow)...")
t0 = time.time()
d1, d2 = diseases[0], diseases[1]
e1, e2 = MechDist_embeddings[d1], MechDist_embeddings[d2]
X1 = np.array([list(e1) + list(e2)])
pred1 = model.predict([X1[:,:128], X1[:,128:]], verbose=0)[0][0]
X2 = np.array([list(e2) + list(e1)])
pred2 = model.predict([X2[:,:128], X2[:,128:]], verbose=0)[0][0]
mean = (pred1 + pred2) / 2
std = np.std([pred1, pred2])
print(f"      ✅ First prediction in {time.time()-t0:.2f}s")
print(f"         {d1} ↔ {d2}")
print(f"         Mean: {mean:.6f}, Std: {std:.6f}")

# Step 5: Second prediction (should be fast)
print("\n[5/5] Running second prediction (should be fast)...")
t0 = time.time()
d1, d2 = diseases[2], diseases[3]
e1, e2 = MechDist_embeddings[d1], MechDist_embeddings[d2]
X1 = np.array([list(e1) + list(e2)])
pred1 = model.predict([X1[:,:128], X1[:,128:]], verbose=0)[0][0]
X2 = np.array([list(e2) + list(e1)])
pred2 = model.predict([X2[:,:128], X2[:,128:]], verbose=0)[0][0]
mean = (pred1 + pred2) / 2
std = np.std([pred1, pred2])
print(f"      ✅ Second prediction in {time.time()-t0:.2f}s")
print(f"         {d1} ↔ {d2}")
print(f"         Mean: {mean:.6f}, Std: {std:.6f}")

print("\n" + "="*60)
print("✅ ALL TESTS PASSED - NOTEBOOK IS FUNCTIONAL")
print("="*60)
