#!/usr/bin/env python3
import sys
from pathlib import Path
from read import FileFormatError, is_bpseq_file, is_dotbracket_file, read_bpseq_file, read_dotbracket_file
from structure import build_structure_map, print_structure_types


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Uso: {sys.argv[0]} <archivo.bpseq|archivo.dbn>")
    path = Path(sys.argv[1])
    if not path.exists():
        sys.exit(f"Archivo no encontrado: {path}")
    try:
        if is_bpseq_file(path):
            bp, seq = read_bpseq_file(path)
            id_ = path.stem
        elif is_dotbracket_file(path):
            bp, seq = read_dotbracket_file(path)
            id_ = path.stem
        else:
            sys.exit("Formato desconocido. Use BPSEQ o dot-bracket.")
    except FileFormatError as e:
        sys.exit(f"Error de formato: {e}")
    dotb,s,k,types,page = build_structure_map(seq,bp)
    print_structure_types(id_,seq,dotb,s,k,types,page,"")

if __name__=="__main__":
    main()