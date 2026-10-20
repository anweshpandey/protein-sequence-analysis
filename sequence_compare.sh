#!/usr/bin/env bash
set -euo pipefail

SEQ1=""
SEQ2=""

while getopts ":1:2:" opt; do
  case $opt in
    1) SEQ1="$OPTARG" ;;
    2) SEQ2="$OPTARG" ;;
  esac
done

[[ -z "$SEQ1" || -z "$SEQ2" ]] && { echo "Usage: $0 -1 <sequence1> -2 <sequence2>"; exit 1; }

# Clean up sequences (remove whitespace and convert to uppercase)
PDBSEQ1=$(echo "$SEQ1" | tr -d ' \t\n\r' | tr '[:lower:]' '[:upper:]')
PDBSEQ2=$(echo "$SEQ2" | tr -d ' \t\n\r' | tr '[:lower:]' '[:upper:]')

python3 - <<'PY' "$PDBSEQ1" "$PDBSEQ2"
import sys
seq1=sys.argv[1].strip()
seq2=sys.argv[2].strip()

match,mismatch,gap=1,-1,-2
n,m=len(seq1),len(seq2)
score=[[0]*(m+1) for _ in range(n+1)]
trace=[[None]*(m+1) for _ in range(n+1)]

for i in range(1,n+1):
    score[i][0]=i*gap; trace[i][0]='U'
for j in range(1,m+1):
    score[0][j]=j*gap; trace[0][j]='L'

for i in range(1,n+1):
    for j in range(1,m+1):
        diag=score[i-1][j-1]+(match if seq1[i-1]==seq2[j-1] else mismatch)
        up=score[i-1][j]+gap
        left=score[i][j-1]+gap
        best=max(diag,up,left)
        score[i][j]=best
        trace[i][j]='D' if best==diag else ('U' if best==up else 'L')

aln_seq1=[]; aln_seq2=[]
i,j=n,m
while i>0 or j>0:
    t=trace[i][j]
    if t=='D': aln_seq1.append(seq1[i-1]); aln_seq2.append(seq2[j-1]); i-=1; j-=1
    elif t=='U': aln_seq1.append(seq1[i-1]); aln_seq2.append('-'); i-=1
    else: aln_seq1.append('-'); aln_seq2.append(seq2[j-1]); j-=1

aln_seq1=list(reversed(aln_seq1))
aln_seq2=list(reversed(aln_seq2))

mismatches=[]; gaps_in_seq1=[]; gaps_in_seq2=[]
seq1_pos=0
seq2_pos=0
for idx,(a,b) in enumerate(zip(aln_seq1, aln_seq2)):
    if a!='-': seq1_pos+=1
    if b!='-': seq2_pos+=1
    
    if a=='-' and b!='-': 
        gaps_in_seq1.append((seq2_pos, b))
    elif a!='-' and b=='-': 
        gaps_in_seq2.append((seq1_pos, a))
    elif a!='-' and b!='-' and a!=b:
        mismatches.append((seq1_pos, seq2_pos, a, b))
        # highlight in aligned sequences
        aln_seq1[idx]='['+a+']'
        aln_seq2[idx]='['+b+']'

print("")
print("=== Sequence Comparison Results ===")
print(f"Sequence 1 length: {len(seq1)}")
print(f"Sequence 2 length: {len(seq2)}")
print()

print("Aligned sequences (differences highlighted):")
print("Sequence 1:", ''.join(aln_seq1))
print("Sequence 2:", ''.join(aln_seq2))

print("\nDifferences:")
if not (mismatches or gaps_in_seq1 or gaps_in_seq2):
    print("Perfect match 🎉")
else:
    if mismatches:
        print("\nMismatches:")
        for pos1, pos2, a, b in mismatches:
            print(f"  - Sequence 1 pos {pos1}: {a} ↔ Sequence 2 pos {pos2}: {b}")
    if gaps_in_seq2:
        print("\nGaps in Sequence 2 (missing residues):")
        for pos, a in gaps_in_seq2:
            print(f"  - Sequence 1 pos {pos}: {a} is missing in Sequence 2")
    if gaps_in_seq1:
        print("\nGaps in Sequence 1 (extra residues in Sequence 2):")
        for pos, b in gaps_in_seq1:
            print(f"  - Sequence 2 pos {pos}: {b} is extra (not in Sequence 1)")

PY
