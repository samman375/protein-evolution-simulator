import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import random

AAMutationMatrix = [["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"],
["A","9867","2","9","10","3","8","17","21","2","6","4","2","6","2","22","35","32","0","2","18"],
["R","1","9914","1","0","1","10","0","0","10","3","1","19","4","1","4","6","1","8","0","1"],
["N","4","1","9822","36","0","4","6","6","21","3","1","13","0","1","2","20","9","1","4","1"],
["D","6","0","42","9859","0","6","53","6","4","1","0","3","0","0","1","5","3","0","0","1"],
["C","1","1","0","0","9973","0","0","0","1","1","0","0","0","0","1","5","1","0","3","2"],
["Q","3","9","4","5","0","9876","27","1","23","1","3","6","4","0","6","2","2","0","0","1"],
["E","10","0","7","56","0","35","9865","4","2","3","1","4","1","0","3","4","2","0","1","2"],
["G","21","1","12","11","1","3","7","9935","1","0","1","2","1","1","3","21","3","0","0","5"],
["H","1","8","18","3","1","20","1","0","9913","0","1","1","0","2","3","1","1","1","4","1"],
["I","2","2","3","1","2","1","2","0","0","9871","9","2","12","7","0","1","7","0","1","33"],
["L","3","1","3","0","0","6","1","1","4","22","9947","2","45","13","3","1","3","4","2","15"],
["K","2","37","25","6","0","12","7","2","2","4","1","9924","20","0","3","8","11","0","1","1"],
["M","1","1","0","0","0","2","0","0","0","5","8","4","9875","1","0","1","2","0","0","4"],
["F","1","1","1","0","0","0","0","1","2","8","6","0","4","9944","0","2","1","3","28","0"],
["P","13","5","2","1","1","8","3","2","5","1","2","2","1","1","9924","12","4","0","0","2"],
["S","28","11","34","7","11","4","6","16","2","2","1","7","4","3","17","9840","38","5","2","2"],
["T","22","2","13","4","1","3","2","2","1","11","2","8","6","1","5","32","9869","0","2","9"],
["W","0","2","0","0","0","0","0","0","0","0","0","0","0","1","0","1","0","9976","1","0"],
["Y","1","0","3","0","3","0","1","0","4","1","1","0","0","21","0","1","1","2","9947","1"],
["V","13","2","1","1","3","2","2","3","3","57","11","1","17","1","3","2","10","0","2","9901"]]

columns = AAMutationMatrix[0]
rows = [row[1:] for row in AAMutationMatrix[1:]]
index = [row[0] for row in AAMutationMatrix[1:]]

df = pd.DataFrame(rows, columns=columns, index=index)

records = 1
for record in SeqIO.parse(sys.stdin, "fasta"):
    if records > 1:
        sys.exit(f"Error: More than one sequence provided. Exiting.")
    
    # Format and print original sequence
    iteration = 0
    record.id = ""
    record.description = "0"
    SeqIO.write(record, sys.stdout, "fasta")
    
    for i in range(500):
        iteration += 1
        newSeq = ""
        # For each amino acid, generate random num < 10000. 
        # Then see which AA mutation it corresponds to in matrix lookup.
        for aa in record.seq:
            if aa not in columns:
                sys.exit(f"Error: Unknown amino acid supplied: \"{aa}\". Exiting")

            psFrame = df[aa]
            rand = random.randrange(1,10000)
            acumulator = 0
            j = 0
            # Loop through probabilities to find corresponding AA
            for p in psFrame:
                acumulator += int(p)
                if rand <= acumulator:
                    # Set AA
                    newLetter = columns[j]
                    newSeq += newLetter
                    break
                j += 1
        
        # Format and print sequence
        record = SeqRecord(
            Seq(newSeq),
            id="",
            description=str(iteration)
        )
        SeqIO.write(record, sys.stdout, "fasta")
    records += 1
