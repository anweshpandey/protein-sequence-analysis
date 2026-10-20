#!/usr/bin/env bash
set -euo pipefail

PDB=""
REFSEQ=""
CHAIN=""

while getopts ":p:s:c:" opt; do
  case $opt in
    p) PDB="$OPTARG" ;;
    s) REFSEQ="$(echo "$OPTARG" | tr -d ' \t\n\r' | tr '[:lower:]' '[:upper:]')" ;;
    c) CHAIN="$OPTARG" ;;
  esac
done

[[ -z "$PDB" || -z "$REFSEQ" ]] && { echo "Usage: $0 -p <pdbfile> -s <sequence> [-c chain]"; exit 1; }

SEQ=$(awk -v want_chain="$CHAIN" '
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
' "$PDB")

[[ "$SEQ" == "ERROR" ]] && { echo "No SEQRES found"; exit 1; }
PDBSEQ=$(echo "$SEQ" | tr -d ' \t\n\r')

python3 - <<'PY' "$REFSEQ" "$PDBSEQ"
import sys
ref=sys.argv[1].strip().upper()
pdb=sys.argv[2].strip().upper()

match,mismatch,gap=1,-1,-2
n,m=len(ref),len(pdb)
score=[[0]*(m+1) for _ in range(n+1)]
trace=[[None]*(m+1) for _ in range(n+1)]

for i in range(1,n+1):
    score[i][0]=i*gap; trace[i][0]='U'
for j in range(1,m+1):
    score[0][j]=j*gap; trace[0][j]='L'

for i in range(1,n+1):
    for j in range(1,m+1):
        diag=score[i-1][j-1]+(match if ref[i-1]==pdb[j-1] else mismatch)
        up=score[i-1][j]+gap
        left=score[i][j-1]+gap
        best=max(diag,up,left)
        score[i][j]=best
        trace[i][j]='D' if best==diag else ('U' if best==up else 'L')

aln_ref=[]; aln_pdb=[]
i,j=n,m
while i>0 or j>0:
    t=trace[i][j]
    if t=='D': aln_ref.append(ref[i-1]); aln_pdb.append(pdb[j-1]); i-=1; j-=1
    elif t=='U': aln_ref.append(ref[i-1]); aln_pdb.append('-'); i-=1
    else: aln_ref.append('-'); aln_pdb.append(pdb[j-1]); j-=1

aln_ref=list(reversed(aln_ref))
aln_pdb=list(reversed(aln_pdb))

missing=[]; substitutions=[]; insertions=[]
ref_pos=0
for idx,(a,b) in enumerate(zip(aln_ref,aln_pdb)):
    if a!='-': ref_pos+=1
    if a=='-' and b!='-': insertions.append((ref_pos,b))
    elif a!='-' and b=='-': missing.append((ref_pos,a))
    elif a!='-' and b!='-' and a!=b:
        substitutions.append((ref_pos,a,b))
        # highlight in aligned sequences
        aln_ref[idx]='['+a+']'
        aln_pdb[idx]='['+b+']'

print("=== Results ===")
print(f"Reference length: {len(ref)}")
print(f"PDB (SEQRES) length: {len(pdb)}\n")

print("Aligned (reference on top, differences highlighted):")
print("ref:", ''.join(aln_ref))
print("pdb:", ''.join(aln_pdb))

print("\nDifferences (relative to reference):")
if not (missing or substitutions or insertions):
    print("Perfect match 🎉")
else:
    if substitutions:
        print("\nSubstitutions:")
        for pos,r,p in substitutions:
            print(f"  - Pos {pos}: {r} → {p}")
    if missing:
        print("\nMissing in PDB:")
        for pos,r in missing:
            print(f"  - Pos {pos}: {r}")
    if insertions:
        print("\nExtra in PDB:")
        for pos,p in insertions:
            print(f"  - After ref pos {pos}: {p}")

PY

