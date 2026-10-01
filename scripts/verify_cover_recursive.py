#!/usr/bin/env python3
"""Exact finite checks for the promoted convex-treelike records.

This replays the pentagonal obstruction and validates the compact census
summaries used by the global-minimality record.  It deliberately does not
duplicate the multi-million-candidate census retained in the source archive.
"""
from collections import deque
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1] / "data" / "examples" / "assets"


def drop(word, coordinate):
    low = (1 << coordinate) - 1
    return (word & low) | ((word >> (coordinate + 1)) << coordinate)


def normalize(family, dimension):
    family = tuple(sorted(set(family)))
    varying = [i for i in range(dimension)
               if {word >> i & 1 for word in family} == {0, 1}]
    projected = tuple(sorted({sum(((word >> old) & 1) << new
                                  for new, old in enumerate(varying))
                              for word in family}))
    return projected, len(varying)


def section(family, dimension, coordinate, value):
    return tuple(sorted({drop(word, coordinate) for word in family
                         if word >> coordinate & 1 == value}))


def is_partial_cube(family, dimension):
    family = tuple(family)
    vertices = set(family)
    adjacency = {x: [y for y in vertices if (x ^ y).bit_count() == 1]
                 for x in vertices}
    for source in vertices:
        distance = {source: 0}
        queue = deque([source])
        while queue:
            x = queue.popleft()
            for y in adjacency[x]:
                if y not in distance:
                    distance[y] = distance[x] + 1
                    queue.append(y)
        if len(distance) != len(vertices):
            return False
        if any(distance[target] != (source ^ target).bit_count()
               for target in vertices):
            return False
    return True


@lru_cache(None)
def is_rper(family, dimension):
    family, dimension = normalize(family, dimension)
    if len(family) == 1:
        return True
    if not is_partial_cube(family, dimension):
        return False
    for coordinate in range(dimension):
        zero = set(section(family, dimension, coordinate, 0))
        one = set(section(family, dimension, coordinate, 1))
        if zero <= one:
            base, site = tuple(sorted(one)), tuple(sorted(zero))
        elif one <= zero:
            base, site = tuple(sorted(zero)), tuple(sorted(one))
        else:
            continue
        if is_rper(base, dimension - 1) and is_rper(site, dimension - 1):
            return True
    return False


def shattered(family, dimension, support):
    return {word & support for word in family} == {
        sum(bit << i for i, bit in zip(
            [j for j in range(dimension) if support >> j & 1], values))
        for values in product((0, 1), repeat=support.bit_count())
    }


def strongly_shattered(family, dimension, support):
    outside = ((1 << dimension) - 1) ^ support
    patterns = {word & support for word in family}
    for base in {word & outside for word in family}:
        if all(base | pattern in family for pattern in patterns) and len(patterns) == 1 << support.bit_count():
            return True
    return False


def run():
    p5 = (0, 1, 2, 3, 4, 6, 7, 8, 12, 14, 16, 17, 19, 24, 25, 28)
    assert is_partial_cube(p5, 5)
    shattering = []
    for support in range(1 << 5):
        s = shattered(p5, 5, support)
        ss = strongly_shattered(set(p5), 5, support)
        assert s == ss
        if s:
            shattering.append(support)
    assert max(x.bit_count() for x in shattering) == 2
    assert len(p5) == 1 + 5 + len(list(combinations(range(5), 2))) == 16

    section_sizes = []
    for coordinate in range(5):
        zero = section(p5, 5, coordinate, 0)
        one = section(p5, 5, coordinate, 1)
        section_sizes.append((len(zero), len(one), len(set(zero) & set(one))))
        assert not set(zero) <= set(one) and not set(one) <= set(zero)
        assert is_rper(zero, 4) and is_rper(one, 4)
        contraction = tuple(sorted(set(zero) | set(one)))
        assert is_rper(contraction, 4)
    assert not is_rper(p5, 5)
    assert set(section_sizes) == {(10, 6, 5)}

    small = json.loads((ROOT / 'cover_recursive_census_small.json').read_text())
    dim5 = small['dimension_five_order_21']
    assert dim5 == {
        'intersections': 103668,
        'assignments': 92927,
        'labelled_partial_cubes': 14259,
        'labelled_periphery_free': 3115,
        'additional_types_hex': ['18fb3ab'],
    }
    assert all(not row['critical_hex'] for d, row in small['pointed_through_dimension_four'].items() if int(d) < 4)
    assert small['pointed_through_dimension_four']['4']['critical_hex'] == ['eacf']

    six = json.loads((ROOT / 'cover_recursive_census_six.json').read_text())
    assert six['input_sha256'] == sha256((ROOT / 'cover_recursive_census_small.json').read_bytes()).hexdigest()
    assert six['base_count'] == 806 and six['types_hex'] == ['1fa888a0bb']

    record = json.loads((ROOT / 'cover_recursive_census_record.json').read_text())
    assert record['passed'] and len(record['chains']) == 10
    absent = [row['dimension'] for row in record['chains']
              if row['kind'] == 'periphery_free' and row['output_count'] == 0]
    assert absent == [8, 9, 10, 11, 12]
    audit = json.loads((ROOT / 'cover_recursive_census_audit.json').read_text())
    assert audit['status'] == 'passed'
    assert audit['minimum_additional_cover_recursive_obstruction_order'] == 16
    assert audit['dimension_five_uniqueness_order_bound'] == 21
    assert audit['periphery_free_search_layers_replayed'] == [7, 8, 9, 10, 11, 12]

    result = {
        'status': 'passed',
        'pentagonal': {
            'words': list(p5), 'order': 16, 'dimension': 5,
            'vc_dimension': 2, 'ample': True, 'periphery_free': True,
            'elementary_minors_checked': 15,
            'section_size_triples': [list(x) for x in section_sizes],
        },
        'census': {
            'dimension_five_partial_cubes': 14259,
            'dimension_five_periphery_free': 3115,
            'dimension_five_additional_types': ['18fb3ab'],
            'minimum_additional_order': 16,
            'periphery_free_layers_replayed': [7, 8, 9, 10, 11, 12],
        },
        'scope': 'Exact replay of P5 and integrity checks of compact exhaustive census reports; the source archive owns the full candidate generation.',
        'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    path = ROOT / 'cover_recursive_verification.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    run()
