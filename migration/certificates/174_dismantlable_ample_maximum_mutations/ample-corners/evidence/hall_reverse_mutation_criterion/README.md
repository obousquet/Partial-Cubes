# Hall reverse-mutation check

Replays all40 legal maximum VC3 exchanges of the recorded Hall class.
This reuses the earlier40-exchange result. Twelve neighbors have two corners,
24 have one, and four have none. All two-corner pairs have distance two;
all reverse removed vertices have degree four. Every reverse exchange
passes the independently proved corner-elimination criterion. Every union
with the added word has300 shattered supports and300 words, hence is ample.

This validates a possible final-step configuration at299 vertices; it does
not produce a smaller obstruction. An explicit neighbor replaces Hall word2
by2049 and has corners1026 and2050. Reversing that exchange kills both
corners without new ones. Source snapshots and hashes are retained.
Reproduce from the repository root using the memory/CPU guardrails and
papers/ample-corners/scripts/verify_hall_reverse_mutations.py.
