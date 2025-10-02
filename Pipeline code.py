#Arin DAS Biopython CA1
import os
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.Align import PairwiseAligner
from Bio.Blast import NCBIWWW, NCBIXML
import matplotlib.pyplot as plt
from Bio.SeqFeature import SeqFeature, FeatureLocation

Entrez.email = "dasarin007@gmail.com" 

#  Utility Functions 

def fetch_record(acc_id, db="nucleotide", rettype="gb"):
    print(f"Fetching {acc_id} from {db}...")
    handle = Entrez.efetch(db=db, id=acc_id, rettype=rettype, retmode="text")
    record = SeqIO.read(handle, rettype)
    handle.close()
    return record

def count_bases(seq):
    total = sum(seq.count(c) for c in "ATGC")
    counts = {base: seq.count(base) for base in "ATGC"}
    return total, counts

def translate_cds(record):
    proteins = []
    cds_features = [f for f in record.features if f.type == "CDS"]
    if cds_features:
        for f in cds_features:
            cds_seq = f.extract(record.seq)
            try:
                prot = cds_seq.translate(to_stop=True)
                proteins.append((f.qualifiers.get("product", ["Unknown"])[0], prot))
            except Exception:
                proteins.append(("Error translating CDS", None))
    else:
        try:
            proteins.append(("Full sequence", record.seq.translate(to_stop=True)))
        except Exception as e:
            proteins.append((f"Translation error: {e}", None))
    return proteins

def run_alignment(seq1, seq2):
    aligner = PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    return alignments[0]

def run_blast(seq, record_id, db="nt", outdir="results", top_hits=5):
    os.makedirs(outdir, exist_ok=True)
    print(f"Running BLAST for {record_id} against {db}...")
    result_handle = NCBIWWW.qblast("blastn", db, str(seq))
    outpath = os.path.join(outdir, f"{record_id}_blast.xml")
    with open(outpath, "w") as f:
        f.write(result_handle.read())
    result_handle.close()
    print("BLAST results saved:", outpath)

    # Parse BLAST XML and print hits
    with open(outpath) as handle:
        blast_record = NCBIXML.read(handle)
        print(f"\nTop {top_hits} hits for {record_id}:")
        for alignment in blast_record.alignments[:top_hits]:
            hsp = alignment.hsps[0]
            print(f"Hit: {alignment.title}")
            print(f"Score: {hsp.score}, E-value: {hsp.expect}")
            print(f"Identities: {hsp.identities}/{hsp.align_length}\n")

def plot_gc_and_orfs(record, outdir="results", show_inline=True):
    import matplotlib.pyplot as plt
    from Bio.SeqUtils import gc_fraction
    import os

    os.makedirs(outdir, exist_ok=True)
    seq = record.seq
    seq_len = len(seq)

    # GC content (sliding window)
    window = 100
    gc_values = [gc_fraction(seq[i:i+window])*100 for i in range(0, seq_len-window, window)]

    # ---- GC plot ----
    plt.figure(figsize=(10,5))
    plt.plot(gc_values, label="GC % (sliding window)")
    plt.xlabel("Window (100bp)")
    plt.ylabel("GC %")
    plt.title(f"GC content for {record.id}")
    plt.legend()
    gc_path = os.path.join(outdir, f"{record.id}_gc.png")
    plt.savefig(gc_path)
    if show_inline:  
        plt.show()
    plt.close()

    # ORF map
    plt.figure(figsize=(10,2))
    for f in record.features:
        if f.type == "CDS":
            start, end = int(f.location.start), int(f.location.end)
            plt.plot([start, end], [1, 1], linewidth=6)
    plt.title(f"ORF map for {record.id}")
    plt.xlabel("Position")
    plt.yticks([])
    orf_path = os.path.join(outdir, f"{record.id}_orf.png")
    plt.savefig(orf_path)
    if show_inline:  
        plt.show()
    plt.close()

    print(f"Plots saved: {gc_path}, {orf_path}")


def annotate_genome(record, outdir="results"):
    os.makedirs(outdir, exist_ok=True)
    gb_path = os.path.join(outdir, f"{record.id}_annotated.gb")
    SeqIO.write(record, gb_path, "genbank")
    print("Annotated GenBank saved:", gb_path)

    # Print features summary
    print("\nFeatures:")
    for f in record.features[:10]:  # first 10 features
        print(f" - {f.type} at {f.location}, qualifiers: {f.qualifiers.get('product', '')}")

# Menu 

def pipeline():
    acc_ids = input("Enter accession IDs (comma separated): ").strip().split(",")
    db_choice = input("Database? (nucleotide/protein): ").strip().lower()
    rettype = "gb" if db_choice == "nucleotide" else "fasta"

    records = [fetch_record(acc_id.strip(), db_choice, rettype) for acc_id in acc_ids]

    while True:
        print("\nChoose an option:")
        print("1. Show sequence (FASTA)")
        print("2. Count ATGC bases")
        print("3. Translate DNA to protein (CDS-aware)")
        print("4. Align two sequences")
        print("5. Visualize GC & ORF maps")
        print("6. Run BLAST")
        print("7. Annotate genome (save GenBank)")
        print("0. Exit")

        choice = input("Enter choice: ").strip()

        if choice == "0":
            break
        elif choice == "1":
            for r in records:
                print(f"\n>{r.id}\n{str(r.seq)}")
        elif choice == "2":
            for r in records:
                total, counts = count_bases(str(r.seq))
                print(f"\n{r.id}:\nTotal ATGC = {total}\nCounts = {counts}")
        elif choice == "3":
            for r in records:
                if db_choice == "nucleotide":
                    prots = translate_cds(r)
                    for desc, prot in prots:
                        print(f"\n{r.id} {desc}:\n{prot}")
                else:
                    print(f"\n{r.id} is protein, translation not needed.")
        elif choice == "4":
            if len(records) >= 2:
                aln = run_alignment(records[0].seq, records[1].seq)
                print("\nAlignment:\n", aln)
            else:
                print("Need at least 2 sequences for alignment.")
        elif choice == "5":
            for r in records:
                if db_choice == "nucleotide":
                    plot_gc_and_orfs(r)
                else:
                    print(f"{r.id} is protein, skipping GC/ORF.")
        elif choice == "6":
            for r in records:
                run_blast(r.seq, r.id, db="nt" if db_choice=="nucleotide" else "nr")
        elif choice == "7":
            for r in records:
                if db_choice == "nucleotide":
                    annotate_genome(r)
                else:
                    print(f"{r.id} is protein, annotation limited.")
        else:
            print("Invalid choice!")

# Run 
if __name__ == "__main__":
    pipeline()
