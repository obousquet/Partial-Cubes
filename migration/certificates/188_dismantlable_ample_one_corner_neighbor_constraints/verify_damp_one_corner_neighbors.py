#!/usr/bin/env python3
"""Check one fixed one-corner obstruction and its adjacent ample additions.

The input is the already certified 291-concept cornerless seed. Verify the
complete placement fibre at gain 2312, then only the adjacent additions to
its fixed extension at 2696. Negative decisions carry complete corner
deletion certificates. This does not test further seeds or expand the
previous pair census, and proves no general repair-inheritance statement.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from explore_damp_cornered_hall_repair_signatures import bounded_peeling
from explore_damp_fixed_291_repair_creation import ample_tip
from verify_damp_singleton_recovery_implications import shadow
from verify_damp_terminal_clause_recovery import corner

ROOT = Path(__file__).resolve().parent


def main():
    source_path = ROOT / 'damp_fixed_291_repair_creation.json'
    source = json.loads(source_path.read_text())
    base = set(source['base'])
    assert len(base) == 291 and len(shadow(base, 4095)) == 291
    assert not any(corner(base, x, 12) for x in base)
    source_model = ROOT / 'damp_hall_carrier_defect_min8.json'
    assert hashlib.sha256(source_model.read_bytes()).hexdigest() == source['source_model_sha256']

    gain, c = 2312, 2696
    tips = [x for x in range(4096) if x not in base and ample_tip(base, x) == gain]
    assert tips == [2184, 2696, 3720, 3818, 3822]
    assert len({x & gain for x in tips}) == 1
    edges = [[x, y] for i, x in enumerate(tips) for y in tips[i + 1:]
             if (x ^ y).bit_count() == 1]
    assert edges == [[2184, 2696], [2696, 3720], [3818, 3822]]
    placements = {x & ~gain for x in tips}
    placement_shadow = shadow(placements, 4095 ^ gain)
    assert len(placement_shadow) > len(placements)
    # Every tip has a neighbor with the same gain, so the proposition applies
    # to all five extensions. Cornerlessness makes each extension nonpeelable.
    for tip in tips:
        family = base | {tip}
        assert [x for x in sorted(family) if corner(family, x, 12)] == [tip]
    family = base | {c}
    checks = []
    for f in range(12):
        tip = c ^ (1 << f)
        if tip in family:
            continue
        support = ample_tip(family, tip)
        if support is None:
            continue
        decision = bounded_peeling(family | {tip}, 10000)
        assert decision['status'] == 'nonpeelable'
        checks.append({'coordinate': f, 'tip': tip, 'support': support, **decision})
    assert [item['tip'] for item in checks] == [2568, 2184, 3720]
    assert [item['visited'] for item in checks] == [3, 5, 4]
    result = {'scope': __doc__.strip(),
              'source_file_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'source_model_sha256': source['source_model_sha256'],
              'base_order': len(base), 'gain': gain, 'complete_gain_fibre': tips,
              'gain_fibre_edges': edges, 'outside_placements': sorted(placements),
              'outside_placement_shadow': sorted(placement_shadow),
              'placement_fibre_is_ample': False, 'all_five_extensions_have_one_corner': True,
              'fixed_source_corner': c, 'source_order': len(family),
              'adjacent_ample_additions': checks, 'all_checks_passed': True}
    (ROOT / 'damp_one_corner_neighbors_verification.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps({'gain_fibre': tips, 'placement_shadow_order': len(placement_shadow),
                      'one_corner_extensions': len(tips),
                      'fixed_source_corner': c,
                      'adjacent_additions': [{k: v for k, v in row.items()
                                               if k in ['tip', 'support', 'status', 'visited']}
                                              for row in checks],
                      'all_checks_passed': True}, indent=2))


if __name__ == '__main__':
    main()
