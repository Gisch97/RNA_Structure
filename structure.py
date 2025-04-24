import networkx as nx
from typing import Dict, List, Tuple

from segments import get_segments, separate_segments, filter_base_pairs
from graph import build_segment_graph

class RNAError(Exception):
    """Errores genéricos en mapeo de estructura."""
    pass


def compute_dotbracket(
    segments: List[List[Tuple[int, int]]],
    knots: List[List[Tuple[int, int]]],
    seq: str
) -> Tuple[str, int]:
    """
    Genera la notación dot-bracket con pseudonudos.
    Retorna (dotbracket, page_number).
    """
    n = len(seq)
    dotb = ['.'] * n
    # Asignar paréntesis a stems
    for stem in segments:
        for l, r in stem:
            dotb[l-1] = '('
            dotb[r-1] = ')'
    # Caracteres para pseudonudos
    lefts = '([{<' + ''.join(chr(i) for i in range(65, 91))  # A-Z
    rights = ')]}>' + ''.join(chr(i) for i in range(97, 123))  # a-z
    page_number = 1
    # Asignar brackets a cada pseudonudo
    for idx, knot in enumerate(knots):
        if idx >= len(lefts):
            raise RNAError(f"Demasiados pseudonudos (>{len(lefts)}) para representar.")
        lchar, rchar = lefts[idx], rights[idx]
        for l, r in knot:
            dotb[l-1] = lchar
            dotb[r-1] = rchar
        page_number = max(page_number, idx + 2)
    return ''.join(dotb), page_number


def fwd_finder(i: int, dotbracket: str, knot_brackets: set) -> Tuple[str, int]:
    """Encuentra siguiente símbolo emparejado o knot a la derecha."""
    n = len(dotbracket)
    while i < n:
        ch = dotbracket[i]
        if ch in '()' or ch in knot_brackets:
            return ch, i
        i += 1
    return '', i


def bwd_finder(i: int, dotbracket: str, knot_brackets: set) -> Tuple[str, int]:
    """Encuentra siguiente símbolo emparejado o knot a la izquierda."""
    while i >= 0:
        ch = dotbracket[i]
        if ch in '()' or ch in knot_brackets:
            return ch, i
        i -= 1
    return '', i


def between(fwd_ip: int, bwd_ip: int, dotbracket: str) -> str:
    """Determina si región entre índices contiene un paréntesis -> X o I."""
    fp = fwd_ip + 1
    while fp < bwd_ip:
        if dotbracket[fp] in '()':
            return 'X'
        fp += 1
    return 'I'


def compute_structure_array(
    dotbracket: str,
    bp: Dict[int, int]
) -> Tuple[List[str], List[int]]:
    """
    Clasifica cada posición en:
      S: stem
      M: multiloop
      I: internal loop
      B: bulge
      H: hairpin loop
      X: external loop
      E: dangling end
    Marca pseudonudos en k_list.
    """
    n = len(dotbracket)
    s_list = [''] * n
    k_list = [0] * n
    knot_brackets = set(
        '[]{}<>' + ''.join(chr(i) for i in range(65,91)) + ''.join(chr(i) for i in range(97,123))
    )
    # Clasificación inicial: stems y pseudonudos
    for i, ch in enumerate(dotbracket):
        if ch in '()':
            s_list[i] = 'S'
        elif ch in knot_brackets:
            s_list[i] = 'S'
            k_list[i] = 1
    # Clasificación de loops
    for i in range(n):
        if s_list[i] == '':
            fwd_ch, fwd_i = fwd_finder(i, dotbracket, knot_brackets)
            bwd_ch, bwd_i = bwd_finder(i, dotbracket, knot_brackets)
            fwd_pair = bp.get(fwd_i + 1)
            bwd_pair = bp.get(bwd_i + 1)
            if bwd_ch == '(':
                if fwd_ch == '(':  # bulge o loop interno
                    s_list[i] = 'B' if (fwd_pair == bwd_pair - 1) else between(fwd_i, bwd_i, dotbracket)
                elif fwd_ch == ')':
                    s_list[i] = 'H'
                else:
                    s_list[i] = 'E'
            elif bwd_ch == ')':
                if fwd_ch == '(':  # región externa
                    s_list[i] = 'X'
                elif fwd_ch == ')':
                    s_list[i] = 'B' if (fwd_pair == bwd_pair - 1) else between(fwd_i, bwd_i, dotbracket)
                else:
                    s_list[i] = 'E'
            else:
                s_list[i] = 'E'
    # Marcar multiloops M (componentes de tamaño >2 en grafo de stems)
    G = nx.Graph()
    for i, ch in enumerate(dotbracket):
        if ch in '()':
            G.add_node(i)
    for i, j in bp.items():
        G.add_edge(i - 1, j - 1)
    for comp in nx.connected_components(G):
        if len(comp) > 2:
            for idx in comp:
                if s_list[idx] == 'S':
                    s_list[idx] = 'M'
    return s_list, k_list


def build_structure_map(
    seq: str,
    bp: Dict[int, int]
) -> Tuple[str, List[str], List[int], Dict[str, List[str]], int]:
    """
    Pipeline completo:
      1. Extraer stems
      2. Separar pseudonudos
      3. Filtrar bp
      4. Recalcular stems limpios
      5. Generar dotbracket + página
      6. Calcular s_list y k_list
      7. Grafo de segmentos y edges
      8. Poblar structure_types para cada categoría
    """
    all_segments = get_segments(bp)
    clean_segs, pk_segs, warnings = separate_segments(all_segments)
    bp_filt = filter_base_pairs(bp, pk_segs)
    segs2 = get_segments(bp_filt)
    dotb, page_number = compute_dotbracket(segs2, pk_segs, seq)
    s_list, k_list = compute_structure_array(dotb, bp_filt)
    G, edges = build_segment_graph(seq, bp_filt, segs2, pk_segs)
    structure_types: Dict[str, List[str]] = {t: [] for t in ['S','H','B','I','M','X','E','PK']}
    # Stems (S)
    for idx, stem in enumerate(segs2, start=1):
        ranges = ",".join(f"{l}-{r}" for l, r in stem)
        structure_types['S'].append(f"S{idx} {ranges}")
    # Pseudonudos (PK)
    for idx, knot in enumerate(pk_segs, start=1):
        ranges = ",".join(f"{l}-{r}" for l, r in knot)
        structure_types['PK'].append(f"PK{idx} {ranges}")
    # Loops y regiones (H, B, I, X, E, M)
    regions: Dict[str, List[Tuple[int, int]]] = {'H':[], 'B':[], 'I':[], 'X':[], 'E':[], 'M':[]}
    start = 0
    while start < len(s_list):
        t = s_list[start]
        end = start
        while end + 1 < len(s_list) and s_list[end + 1] == t:
            end += 1
        if t in regions:
            regions[t].append((start + 1, end + 1))
        start = end + 1
    for t, segs in regions.items():
        for idx, (a, b) in enumerate(segs, start=1):
            seq_sub = seq[a-1:b]
            structure_types[t].append(f"{t}{idx} {a}..{b} \"{seq_sub}\"")
    return dotb, s_list, k_list, structure_types, page_number


def print_structure_types(
    id_: str,
    seq: str,
    dotbracket: str,
    s_list: List[str],
    k_list: List[int],
    structure_types: Dict[str, List[str]],
    page_number: int,
    warnings: str
) -> None:
    """
    Escribe archivo .st completo.
    """
    with open(f"{id_}.st", 'w') as f:
        f.write(f"#Name: {id}")
        f.write(f"#Length: {len(seq)}")
        f.write(f"#PageNumber: {page_number}")
        if warnings:
            f.write(warnings)
        f.write(seq + "")
        f.write(dotbracket + "")
        f.write(''.join(s_list) + "")
        f.write(''.join('K' if ki else 'N' for ki in k_list) + "")
        for key in ['S','H','B','I','M','X','E','PK']:
            for line in structure_types.get(key, []):
                f.write(line)