#!/usr/bin/env bash
set -euo pipefail

PDB1=""
PDB2=""
CHAIN1=""
CHAIN2=""

while getopts ":a:b:c:d:" opt; do
  case $opt in
    a) PDB1="$OPTARG" ;;
    b) PDB2="$OPTARG" ;;
    c) CHAIN1="$OPTARG" ;;
    d) CHAIN2="$OPTARG" ;;
  esac
done

[[ -z "$PDB1" || -z "$PDB2" ]] && { echo "Usage: $0 -a <pdb1file> -b <pdb2file> [-c chain1] [-d chain2]"; exit 1; }

# Function to extract SEQRES from PDB file
extract_seqres() {
    local pdb_file="$1"
    local want_chain="$2"
    
    awk -v want_chain="$want_chain" '
    BEGIN{
      aa["ALA"]="A"; aa["ARG"]="R"; aa["ASN"]="N"; aa["ASP"]="D"; aa["CYS"]="C";
      aa["GLN"]="Q"; aa["GLU"]="E"; aa["GLY"]="G"; aa["HIS"]="H"; aa["ILE"]="I";
      aa["LEU"]="L"; aa["LYS"]="K"; aa["MET"]="M"; aa["PHE"]="F"; aa["PRO"]="P";
      aa["SER"]="S"; aa["THR"]="T"; aa["TRP"]="W"; aa["TYR"]="Y"; aa["VAL"]="V";
      aa["MSE"]="M"; aa["SEC"]="U"; aa["PYL"]="O";
      aa["HSD"]="H"; aa["HSE"]="H"; aa["HSP"]="H";
    }
    function map1(x){ u=toupper(x); return (u in aa)?aa[u]:"X"; }
    $1=="SEQRES"{
      cid=(length($3)==1)?$3:substr($0,12,1);
      if(target_chain==""){
        if(want_chain!=""){ target_chain=want_chain } else { target_chain=cid }
      }
      if(cid!=target_chain) next;
      for(i=5;i<=NF;i++){ if($i!="") seq=seq map1($i) }
    }
    END{ if(seq=="") print "ERROR"; else print seq }
    ' "$pdb_file"
}

SEQ1=$(extract_seqres "$PDB1" "$CHAIN1")
[[ "$SEQ1" == "ERROR" ]] && { echo "No SEQRES found in $PDB1"; exit 1; }
PDBSEQ1=$(echo "$SEQ1" | tr -d ' \t\n\r')

SEQ2=$(extract_seqres "$PDB2" "$CHAIN2")
[[ "$SEQ2" == "ERROR" ]] && { echo "No SEQRES found in $PDB2"; exit 1; }
PDBSEQ2=$(echo "$SEQ2" | tr -d ' \t\n\r')

python3 - <<'PY' "$PDB1" "$PDB2" "$PDBSEQ1" "$PDBSEQ2"
import sys
pdb1_name=sys.argv[1]
pdb2_name=sys.argv[2]
seq1=sys.argv[3].strip().upper()
seq2=sys.argv[4].strip().upper()

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

print("=== PDB SEQRES Comparison Results ===")
print(f"{pdb1_name} length: {len(seq1)}")
print(f"{pdb2_name} length: {len(seq2)}")
print()

print("Aligned sequences (differences highlighted):")
print(f"{pdb1_name}:", ''.join(aln_seq1))
print(f"{pdb2_name}:", ''.join(aln_seq2))

print("\nDifferences:")
if not (mismatches or gaps_in_seq1 or gaps_in_seq2):
    print("Perfect match 🎉")
else:
    if mismatches:
        print("\nMismatches:")
        for pos1, pos2, a, b in mismatches:
            print(f"  - {pdb1_name} pos {pos1}: {a} ↔ {pdb2_name} pos {pos2}: {b}")
    if gaps_in_seq2:
        print("\nGaps in", pdb2_name, "(missing residues):")
        for pos, a in gaps_in_seq2:
            print(f"  - {pdb1_name} pos {pos}: {a} is missing in {pdb2_name}")
    if gaps_in_seq1:
        print("\nGaps in", pdb1_name, "(extra residues in", pdb2_name + "):")
        for pos, b in gaps_in_seq1:
            print(f"  - {pdb2_name} pos {pos}: {b} is extra (not in {pdb1_name})")

PY
