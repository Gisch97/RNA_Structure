import networkx as nx
from typing import Dict,List,Tuple,Optional
from segments import get_extreme_positions, get_next_pair, get_prev_pair

class GraphError(Exception): pass

def build_segment_graph(
    seq: str, bp: Dict[int,int], segments: List[List[Tuple[int,int]]], knots: Optional[List[List[Tuple[int,int]]]]=None
) -> Tuple[nx.MultiDiGraph,List[Tuple[int,int,int,int,str]]]:
    knots=knots or []
    first_pos,last_pos=get_extreme_positions(bp)
    G=nx.MultiDiGraph(); edges=[]
    for i in range(len(segments)): G.add_node(i)
    for i,seg in enumerate(segments):
        if not seg: continue
        s1_5p,s1_3p=seg[0]
        s1_3p_stop,s1_5p_stop=seg[-1]
        p1_5p_start=get_prev_pair(s1_5p,bp,first_pos,knots)
        n1_3p_stop =get_next_pair(s1_3p_stop,bp,last_pos,knots)
        p1_5p_stop =get_prev_pair(s1_5p_stop,bp,first_pos,knots)
        n1_3p_start=get_next_pair(s1_3p,bp,last_pos,knots)
        for j,seg2 in enumerate(segments):
            if not seg2: continue
            s2_5p,s2_3p=seg2[0]; s2_3p_stop,s2_5p_stop=seg2[-1]
            for lbl,(a,b) in [("1",(n1_3p_start,s2_5p)),
                             ("2",(n1_3p_start,s2_5p_stop)),
                             ("3",(n1_3p_stop,s2_5p)),
                             ("4",(n1_3p_stop,s2_5p_stop))]:
                if a==b:
                    G.add_edge(i,j)
                    edges.append((i,j,a,b,lbl))
    return G,edges