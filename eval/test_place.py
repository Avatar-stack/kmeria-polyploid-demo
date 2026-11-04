#!/usr/bin/env python3
"""Self-check for evaluate.place(): k-mers in a duplicated segment report both copies,
reverse complements are found, absent k-mers stay unplaced. Run: python eval/test_place.py"""
import os
import numpy as np

src = open(os.path.join(os.path.dirname(__file__), "evaluate.py")).read()
src = src[:src.index("idx = None")]
src = src.replace("W = sys.argv[1]; ", 'W = "/nonexistent"; ').replace(
    'truth = json.load(open(f"{W}/sim/truth.json")); seed = truth["args"]["seed"]', "")
g = {}; exec(src, g)

rng = np.random.default_rng(0)
seq = bytes(rng.choice(list(b"ACGT"), 5000).astype(np.uint8))
seq = seq[:3000] + seq[1000:1100] + seq[3100:]          # identical copy of 1000-1100 at 3000-3100
codes = g["codes"](seq); pos = np.arange(len(codes), dtype=np.uint32)
o = np.lexsort((pos, codes)); idx = (codes[o], pos[o], np.unique(codes))
shared, uniq = seq[1010:1041].decode(), seq[200:231].decode()
rc = uniq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
first, last, npos, in_ref = g["place"]([shared, uniq, rc, "A" * 31], idx)
assert (first[0], last[0], npos[0]) == (1010, 3010, 2)
assert first[1] == last[1] == first[2] == 200 and npos[1] == 1
assert first[3] == -1 and not in_ref[3]
print("place() ok")
