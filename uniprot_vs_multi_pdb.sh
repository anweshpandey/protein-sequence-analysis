#!/usr/bin/env bash
set -euo pipefail

REFSEQ=""
PDB_FILES=()
CHAIN=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -s)
            REFSEQ="$(echo "$2" | tr -d ' \t\n\r' | tr '[:lower:]' '[:upper:]')"
            shift 2
            ;;
        -c)
            CHAIN="$2"
            shift 2
            ;;
        -*) 
            # Assume any argument starting with - is a PDB file flag followed by filename
            if [[ $# -ge 2 ]]; then
                PDB_FILES+=("$2")
                shift 2
            else
                echo "Error: Missing filename for flag $1"
                exit 1
            fi
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

[[ -z "$REFSEQ" || ${#PDB_FILES[@]} -eq 0 ]] && { 
    echo "Usage: $0 -s <sequence> -1 <pdb1> [-2 <pdb2> -3 <pdb3> ...] [-c chain]"
    echo "Example: $0 -s \"ACDEFGH\" -1 file1.pdb -2 file2.pdb -3 file3.pdb -c A"
    exit 1
}

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

echo ""
echo "===== Comparing reference sequence with ${#PDB_FILES[@]} PDB files ====="
echo "===== Reference sequence length: ${#REFSEQ}"
echo ""

for ((i=0; i<${#PDB_FILES[@]}; i++)); do
    PDB_FILE="${PDB_FILES[$i]}"
    PDB_NAME=$(basename "$PDB_FILE")
    
    echo ""
    echo "***** Processing: $PDB_NAME:"
    
    # Extract SEQRES from PDB
    SEQ=$(extract_seqres "$PDB_FILE" "$CHAIN")
    if [[ "$SEQ" == "ERROR" ]]; then
        echo "  No SEQRES found in $PDB_NAME"
        echo ""
        continue
    fi
    
    PDBSEQ=$(echo "$SEQ" | tr -d ' \t\n\r')
    
    # Compare using Python
    python3 - <<'PY' "$REFSEQ" "$PDBSEQ" "$PDB_NAME"
import sys
ref=sys.argv[1].strip().upper()
pdb=sys.argv[2].strip().upper()
pdb_name=sys.argv[3]

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
        aln_ref[idx]='['+a+']'
        aln_pdb[idx]='['+b+']'

print(f"  {pdb_name} length: {len(pdb)}")

if not (missing or substitutions or insertions):
    print("  Perfect match 🎉")
else:
    print("  Aligned sequences:")
    print("  ref:", ''.join(aln_ref))
    print("  pdb:", ''.join(aln_pdb))
    
    if substitutions:
        print("  Substitutions:")
        for pos,r,p in substitutions:
            print(f"    - Pos {pos}: {r} → {p}")
    if missing:
        print("  Missing in PDB:")
        for pos,r in missing:
            print(f"    - Pos {pos}: {r}")
    if insertions:
        print("  Extra in PDB:")
        for pos,p in insertions:
            print(f"    - After ref pos {pos}: {p}")

# Calculate identity
identical = sum(1 for a,b in zip(ref, pdb) if a == b)
total = min(len(ref), len(pdb))
if total > 0:
    identity_pct = (identical / total) * 100
    print(f"  Identity: {identity_pct:.1f}% ({identical}/{total} residues)")
else:
    print("  Identity: Cannot calculate (zero length sequence)")

PY
    
    echo ""
done
