pip install Bio
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Read the sequences using SeqIO
sequences = list(SeqIO.parse('/content/biology_dataset.txt', 'fasta'))

# Create a dictionary to store unique IDs
unique_ids = {}
processed_sequences = []

# Process sequences and ensure unique IDs
for i, seq_record in enumerate(sequences):
    # Create a unique ID if it's a duplicate
    original_id = seq_record.id
    if original_id in unique_ids:
        unique_ids[original_id] += 1
        seq_record.id = f"{original_id}_{unique_ids[original_id]}"
    else:
        unique_ids[original_id] = 0

    # Choose trimming or padding (uncomment the desired method)
    method = "trim"  # or "pad"

    # Process sequences based on the chosen method
    if method == "trim":
        # Trim sequences to the shortest length
        shortest_length = min(len(seq) for seq in sequences) # Calculate shortest length here
        trimmed_seq = seq_record.seq[:shortest_length]
        processed_sequences.append(SeqRecord(trimmed_seq, id=seq_record.id, description=seq_record.description))
    elif method == "pad":
        # Pad sequences to the longest length with gaps (-)
        longest_length = max(len(seq) for seq in sequences) # Calculate longest length here
        padded_seq = seq_record.seq.ljust(longest_length, '-')  # Use longest_length for padding
        processed_sequences.append(SeqRecord(padded_seq, id=seq_record.id, description=seq_record.description))

# Create a MultipleSeqAlignment object from the processed sequences
align = AlignIO.MultipleSeqAlignment(processed_sequences)

# Print the alignment
print(align)
