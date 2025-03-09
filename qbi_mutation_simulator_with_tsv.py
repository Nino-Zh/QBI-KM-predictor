import random
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

def mutate_sequence(original_seq, mutation_prob):
    """Introduce substitutions and track mutations"""
    mutated = []
    mutations = []
    for pos, aa in enumerate(original_seq, 1):  # 1-based position
        if aa not in AMINO_ACIDS:
            mutated.append(aa)
            continue
            
        if random.random() < mutation_prob:
            new_aa = random.choice([a for a in AMINO_ACIDS if a != aa])
            mutated.append(new_aa)
            mutations.append({
                'position': pos,
                'original': aa,
                'mutated': new_aa
            })
        else:
            mutated.append(aa)
    
    return ''.join(mutated), mutations

def process_fasta(input_file, mutation_prob, num_variants, fasta_output=None, tsv_output=None):
    """Process FASTA and generate mutation report"""
    variant_records = []
    mutation_report = []
    
    for record in SeqIO.parse(input_file, "fasta"):
        original_seq = str(record.seq)
        seq_id = record.id
        
        for variant_num in range(1, num_variants + 1):
            attempts = 0
            while True:
                mutated_seq, mutations = mutate_sequence(original_seq, mutation_prob)
                attempts += 1
                if len(mutations) > 0:
                    break

            variant_id = f"{seq_id}_mutated_variant{variant_num}"
            variant_records.append(SeqRecord(
                seq=mutated_seq,
                id=variant_id,
                description=f"Original: {seq_id} | MutationRate: {mutation_prob}"
            ))

            # Record each mutation in TSV
            for mutation in mutations:
                mutation_report.append({
                    'VariantID': variant_id,
                    'OriginalID': seq_id,
                    'Position': mutation['position'],
                    'OriginalAA': mutation['original'],
                    'MutatedAA': mutation['mutated'],
                    'MutationRate': mutation_prob,
                    'Attempts': attempts
                })
     # Write FASTA output
    if fasta_output is not None:
        SeqIO.write(variant_records, fasta_output, "fasta")
        
    # Write TSV output
    if tsv_output is not None:
        
        with open(tsv_output, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'VariantID', 'OriginalID', 'Position', 
                'OriginalAA', 'MutatedAA', 'MutationRate', 'Attempts'
            ], delimiter='\t')
            writer.writeheader()
            writer.writerows(mutation_report)

    return variant_records, mutation_report

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AA Mutation Simulator with Tracking")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-of", "--fasta_output", required=True, help="Output FASTA file")
    parser.add_argument("-ot", "--tsv_output", required=True, help="Mutation report TSV file")
    parser.add_argument("-m", "--mutation_rate", type=float, default=0.01,
                      help="Mutation probability per residue (default: 0.01)")
    parser.add_argument("-n", "--num_variants", type=int, default=1,
                      help="Number of variants per sequence (default: 1)")
    
    args = parser.parse_args()
    
    process_fasta(
        input_file=args.input,
        fasta_output=args.fasta_output,
        tsv_output=args.tsv_output,
        mutation_prob=args.mutation_rate,
        num_variants=args.num_variants
    )