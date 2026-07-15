from itertools import chain, count
from operator import itemgetter
from typing import Iterable, Optional, Tuple

from graph_tool import Graph
from graph_tool.topology import sequential_vertex_coloring


def run_length_encode(text: str) -> str:
    if not text:
        return ''
    encoding = []
    run_length = 1
    prev_char = text[0]
    for char in text[1:]:
        if char == prev_char:
            run_length += 1
        else:
            encoding.append(f'{run_length}{prev_char}')
            prev_char = char
            run_length = 1
    encoding.append(f'{run_length}{prev_char}')
    return ','.join(encoding)


def run_length_decode(encoding: str) -> str:
    return ''.join(int(token[:-1]) * token[-1] for token in encoding.split(',')) if encoding else ''


def get_interval_overlap_graph(intervals: Iterable[Tuple[int, int]], labels: Optional[Iterable] = None, label_type: str = 'string') -> 'Graph':
    # inspired by https://stackoverflow.com/a/19088519
    # build graph of labels where labels are adjacent if their intervals overlap
    if not labels:
        labels = count()
    g = Graph(directed=False)
    g.vp.label = g.new_vertex_property(label_type)
    label_to_vertex = dict()
    active_labels = set()
    boundaries = sorted(
        chain.from_iterable(
            [(a, True, label), (b, False, label)]
            for (a, b), label in zip(intervals, labels)
        ),
        key = itemgetter(0, 1)
    )
    for _, start_of_interval, label in boundaries:
        if start_of_interval:
            if label not in label_to_vertex:
                v = g.add_vertex()
                g.vp.label[v] = label
                label_to_vertex[label] = v
            for other_label in active_labels:
                i = label_to_vertex[label]
                j = label_to_vertex[other_label]
                g.add_edge(i, j)
            active_labels.add(label)
        else:
            active_labels.discard(label)
    return g, label_to_vertex


def generate_subtracks(intervals: Iterable[Tuple[int, int]], labels: Iterable):
    # inspired by https://stackoverflow.com/a/19088519
    # build graph of labels where labels are adjacent if their intervals overlap
    g, vertex_labels = get_interval_overlap_graph(intervals, labels)
    # find vertex coloring of graph
    # all labels w/ same color can be put into same subtrack
    coloring = sequential_vertex_coloring(g)
    label_to_subtrack = dict(zip(vertex_labels, coloring))
    subtracks = max(label_to_subtrack.values(), default=0) + 1
    return label_to_subtrack, subtracks
