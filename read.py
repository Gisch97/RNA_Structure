#!/usr/bin/env python3
"""
Módulo para detectar y leer archivos BPSEQ o dot-bracket (.db/.dbn).
"""
import re
from pathlib import Path
from typing import Dict, Tuple, List

class FileFormatError(Exception):
    pass


def is_bpseq_file(path: Path) -> bool:
    """Detecta formato BPSEQ: líneas con 3 columnas (i, base, j)."""
    with path.open() as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) != 3:
                return False
    return True


def is_dotbracket_file(path: Path) -> bool:
    """Detecta formato dot-bracket: 2 ó 3 líneas (opcional header), secuencia y bracket igual longitud."""
    lines = [ln.strip() for ln in path.open() if not ln.startswith('#') and ln.strip()]
    if len(lines) not in (2, 3):
        return False
    if len(lines) == 3 and not lines[0].startswith('>'):
        return False
    seq = lines[-2]
    dotb = lines[-1].split()[0]
    return len(seq) == len(dotb)


def read_bpseq_file(path: Path) -> Tuple[Dict[int,int], str]:
    """
    Lee .bpseq y retorna (bp_dict, seq_string).
    bp_dict: llave=posición 1-based, valor=posición pareja.
    """
    bp: Dict[int,int] = {}
    seen: Dict[int,int] = {}
    seq_chars: List[str] = []
    with path.open() as f:
        for num, line in enumerate(f, start=1):
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) != 3:
                raise FileFormatError(f"Línea {num}: se esperaban 3 columnas, hay {len(parts)}")
            i, b, j = int(parts[0]), parts[1], int(parts[2])
            seq_chars.append(b)
            if j:
                if i == j:
                    raise FileFormatError(f"Línea {num}: posición {i} emparejada consigo misma")
                if i in bp or j in seen:
                    raise FileFormatError(f"Línea {num}: emparejamiento inconsistente en {i},{j}")
                bp[i] = j
                seen[j] = i
    return bp, ''.join(seq_chars)


def read_dotbracket_file(path: Path) -> Tuple[Dict[int,int], str]:
    """
    Lee .db/.dbn y retorna (bp_dict, seq_string), usando pair mapping.
    """
    from collections import defaultdict
    bracket_pairs = {'(': ')', '[': ']', '{': '}', '<': '>'}
    open_stack: Dict[str,List[int]] = defaultdict(list)
    bp: Dict[int,int] = {}
    lines = [ln.strip() for ln in path.open() if not ln.startswith('#') and ln.strip()]
    seq = lines[-2]
    dotb = lines[-1].split()[0]
    for idx, ch in enumerate(dotb, start=1):
        if ch in bracket_pairs:
            open_stack[ch].append(idx)
        else:
            for oc, cc in bracket_pairs.items():
                if ch == cc:
                    if not open_stack[oc]:
                        raise FileFormatError(f"Bracket sin abrir en posición {idx}")
                    i = open_stack[oc].pop()
                    bp[i] = idx
                    bp[idx] = i
                    break
    return bp, seq