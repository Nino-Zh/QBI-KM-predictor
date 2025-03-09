import torch
import esm
import pandas as pd

# load csv data containing mutation sequences and km values
data = pd.read_csv("mutations.csv", sep="\t") 


mutations = data["Amino Acid Sequence"].tolist() 
km_values = data["KM (mM)"].tolist() 

# loading esm2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results


# formatting data as a list of tuples (label, sequence) for each mutation
data = [(f"protein{i+1}", mutation) for i, mutation in enumerate(mutations)]

# converting the sequences into model input format
batch_labels, batch_strs, batch_tokens = batch_converter(data)
batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

# extraacting per-residue representations
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=True)
token_representations = results["representations"][33]


sequence_representations = []
for i, tokens_len in enumerate(batch_lens):
    sequence_representations.append(token_representations[i, 1 : tokens_len - 1].mean(0))

# in theory we want to find the sequence with the lowest KM value as "best candidate"
best_candidate_index = km_values.index(min(km_values))  # finding index with lowest KM value
best_sequence = mutations[best_candidate_index]
best_km_value = km_values[best_candidate_index]

print(f"The best evolutionary candidate is protein{best_candidate_index + 1} with KM value {best_km_value}.")
print(f"Sequence: {best_sequence}")

# plot further here
