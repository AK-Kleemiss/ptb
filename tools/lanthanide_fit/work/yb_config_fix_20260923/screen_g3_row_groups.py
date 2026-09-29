"""Screen coordinated Yb-row backoffs for a size-feasible, ion-correct branch."""
import itertools
import json
import tempfile
from pathlib import Path

import screen_g3_backoffs as s

OUT = s.HERE / "ion_g3_row_groups_20260929"
OUT.mkdir(exist_ok=True)
RESULT = OUT / "result.json"
ROWS = (7, 8, 9, 10, 12, 13, 14)
SUBSETS = [comb for count in range(1, len(ROWS)+1)
           for comb in itertools.combinations(ROWS, count)]


def main():
    records = []
    with tempfile.TemporaryDirectory(prefix="yb_groups_") as temp:
        base = Path(temp)
        for index, subset in enumerate(SUBSETS, 1):
            work = base / f"case_{index:03d}"
            work.mkdir()
            lines = s.b.copy()
            for row in subset:
                lines[s.start+row] = s.a[s.start+row]
            candidate = work / "candidate.atompara"
            candidate.write_text("\n".join(lines) + "\n", encoding="ascii")
            record = {"reverted_rows": list(subset), "candidate_sha256": s.sha(candidate)}
            try:
                violation, populations = s.gate_violation(s.PTB, candidate, s.BAS, s.TARGETS)
                record["violation"] = violation
                record["populations"] = populations
                if violation == 0.0:
                    record["neutral_radius"] = s.neutral_radius(candidate, work)
                    record["passes_radius"] = record["neutral_radius"] <= s.LIMIT_Q0
                    if record["passes_radius"]:
                        (OUT / f"case_{index:03d}.atompara").write_bytes(candidate.read_bytes())
            except Exception as exc:
                record["error"] = repr(exc)
            records.append(record)
            RESULT.write_text(json.dumps({"g2_sha256": s.sha(s.G2), "g3_sha256": s.sha(s.G3),
                                          "ptb_sha256": s.sha(s.PTB), "basis_sha256": s.sha(s.BAS),
                                          "q0_radius_limit": s.LIMIT_Q0, "completed": len(records),
                                          "total": len(SUBSETS), "records": records}, indent=2) + "\n")
            if index % 10 == 0 or record.get("passes_radius"):
                print(index, len(SUBSETS), subset, record.get("violation"),
                      record.get("neutral_radius"), record.get("passes_radius"), flush=True)


if __name__ == "__main__":
    main()
