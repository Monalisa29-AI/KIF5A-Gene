from Bio import SeqIO
from collections import Counter

def read_fna_file(file_path):
    """Reads a .fna file and returns the nucleotide sequence."""
    for record in SeqIO.parse(file_path, "fasta"):
        return str(record.seq)  # Return the first sequence found

def count_repeated_sequences(sequence, length):
    """Counts repeated nucleotide sequences of a given length."""
    if length <= 0:
        raise ValueError("Length must be a positive integer.")
    
    # Use a Counter to count occurrences of each sequence
    sequence_counts = Counter()
    
    # Slide over the sequence to extract substrings of the specified length
    for i in range(len(sequence) - length + 1):
        subsequence = sequence[i:i + length]
        sequence_counts[subsequence] += 1
    
    # Filter out sequences that occur only once
    repeated_sequences = {seq: count for seq, count in sequence_counts.items() if count > 1}
    
    return repeated_sequences

def main():
    file_path = 'gene.fna'  # Path to your .fna file
    sequence = read_fna_file(file_path)

    if sequence:
        length = 20

         # Specify the length of the nucleotide sequence to count
        repeated_sequences = count_repeated_sequences(sequence, length)
        
        print(f"Total Length of Sequence: {len(sequence)}")
        print(f"Repeated Sequences of Length {length}:")
        for seq, count in repeated_sequences.items():
            print(f"{seq}: {count} times")
    else:
        print("No sequence found in the file.")

if __name__ == "__main__":
    main()