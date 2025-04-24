import networkx as nx
from typing import Dict, List, Tuple

class RNAError(Exception): pass

def get_extreme_positions(bp: Dict[int,int]) -> Tuple[int,int]:
    if not bp:
        raise RNAError("No base pairs found.")
    keys = sorted(bp.keys())
    return keys[0], keys[-1]

def in_knot(pos: int, knots: List[List[Tuple[int,int]]] = None) -> bool:
    knots = knots or []
    for knot in knots:
        first5, _ = knot[0]
        _, last5 = knot[-1]
        if (first5 <= pos <= last5) or (first5 >= pos >= last5):
            return True
    return False

def get_next_pair(i: int, bp: Dict[int,int], last_pos: int, knots: List[List[Tuple[int,int]]] = None) -> int:
    knots = knots or []
    for n in range(i+1, last_pos+1):
        if not in_knot(n, knots) and bp.get(n,0):
            return n
    return 0

def get_prev_pair(j: int, bp: Dict[int,int], first_pos: int, knots: List[List[Tuple[int,int]]] = None) -> int:
    knots = knots or []
    for p in range(j-1, first_pos-1, -1):
        if not in_knot(p, knots) and bp.get(p,0):
            return p
    return 0

def get_segments(bp: Dict[int,int], knots: List[List[Tuple[int,int]]] = None) -> List[List[Tuple[int,int]]]:
    first_pos, last_pos = get_extreme_positions(bp)
    knots = knots or []
    segments: List[List[Tuple[int,int]]] = []
    i = get_next_pair(0, bp, last_pos, knots)
    while i:
        j = bp[i]
        if first_pos <= i <= last_pos and i<j:
            seg=[(i,j)]
            while True:
                n = get_next_pair(i,bp,last_pos,knots)
                p = get_prev_pair(j,bp,first_pos,knots)
                if n and p and bp.get(n)==p and n<p:
                    seg.append((n,p)); i,j=n,p
                else:
                    break
            segments.append(seg)
        i = get_next_pair(i,bp,last_pos,knots)
    return segments

def filter_base_pairs(bp: Dict[int,int], knots: List[List[Tuple[int,int]]] = None) -> Dict[int,int]:
    knots = knots or []
    return {i:j for i,j in bp.items() if not in_knot(i,knots) and not in_knot(j,knots)}

def pk_quartet(i,j,k,l) -> bool:
    return (i<k<j<l) or (k<i<l<j)

def separate_segments(segments: List[List[Tuple[int,int]]]) -> Tuple[List,List,str]:
    G=nx.Graph(); n=len(segments)
    G.add_nodes_from(range(n))
    for a in range(n):
        for b in range(a+1,n):
            if pk_quartet(*segments[a][0],*segments[b][0]): G.add_edge(a,b)
    clean_idxs=set(range(n)); pk_idxs=set(); warnings=""
    for comp in nx.connected_components(G):
        if len(comp)>1:
            max_idx=max(comp,key=lambda x:len(segments[x]))
            for idx in comp:
                if idx!=max_idx:
                    pk_idxs.add(idx)
                    warnings+=f"#Warning: segment {idx} assigned to pseudoknot\n"
            clean_idxs-=(comp-{max_idx})
    clean=[segments[i] for i in sorted(clean_idxs)]
    pk=[segments[i] for i in sorted(pk_idxs)]
    return clean, pk, warnings