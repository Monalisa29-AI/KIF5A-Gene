from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import re

def read_fna_file(file_path):
    """Reads a .fna file and returns the nucleotide sequence."""
    for record in SeqIO.parse(file_path, "fasta"):
        return str(record.seq)  # Return the first sequence found

def calculate_gc_content(sequence):
    """Calculates the GC content of a nucleotide sequence."""
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    gc_content = (g_count + c_count) / len(sequence) * 100
    return gc_content  
def calculate_at_content(sequence):
    """Calculates the GC content of a nucleotide sequence."""
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    at_content = (a_count + t_count) / len(sequence) * 100
    return at_content

def count_nucleotides(sequence):
    """Counts the occurrences of each nucleotide in the sequence."""
    return Counter(sequence)

def main():
    file_path = 'gene.fna'  # Path to your .fna file
    sequence = read_fna_file(file_path)

    if sequence:
        print(f"Sequence Length: {len(sequence)}")
        print(f"GC Content: {calculate_gc_content(sequence):.2f}%")
        print(f"AT Content: {calculate_at_content(sequence):.2f}%")
        nucleotide_counts = count_nucleotides(sequence)
        print("Nucleotide Counts:")
        for nucleotide, count in nucleotide_counts.items():
            print(f"{nucleotide}: {count}")
    else:
        print("No sequence found in the file.")


def plot_nucleotide_frequencies(counts):
    labels = counts.keys()
    values = counts.values()

    plt.bar(labels, values, color=['blue', 'green', 'red', 'orange'])
    plt.title("Nucleotide Frequency")
    plt.xlabel("Nucleotide")
    plt.ylabel("Frequency")
    plt.show()

# Example usage

 
# Load the sequence from the gene.fna file
def load_sequence(file_path):
    with open(file_path, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            return record.seq  # Return the first sequence found

# Example usage
gene_sequence = load_sequence("gene.fna")
print(f"Loaded sequence: {gene_sequence[:100]}")  # Display the first 100 bases
 
#The line below had an extra space causing the IndentationError
def count_nucleotides(sequence):
    return Counter(sequence)

# Example usage
nucleotide_counts = count_nucleotides(gene_sequence)
print(f"Nucleotide counts: {nucleotide_counts}")
    
    
def plot_nucleotide_frequencies(counts):
    labels = counts.keys()
    values = counts.values()

    plt.bar(labels, values, color=['blue', 'green', 'red', 'orange'])
    plt.title("Nucleotide Frequency")
    plt.xlabel("Nucleotide")
    plt.ylabel("Frequency")
    plt.show()

# Example usage
plot_nucleotide_frequencies(nucleotide_counts)
    
    
def find_motifs(sequence, motif):
    pattern = re.compile(motif)
    matches = pattern.finditer(str(sequence))
    return [(match.start(), match.end()) for match in matches]


    
    
    
   
if __name__ == "__main__": 
    main()