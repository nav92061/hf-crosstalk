"""Recover volumes.csv from a Kaggle/Colab run LOG when the output files were lost.

The run echoes each completed row as a line beginning "ROWCSV\t". Paste the whole log
into run_log.txt and this reconstructs the table.

Usage:  python parse_run_log.py run_log.txt volumes_recovered.csv
"""
import sys, csv, io

def parse(log_text):
    rows = [l.split("\t", 1)[1] for l in log_text.splitlines() if l.startswith("ROWCSV\t")]
    if not rows:
        raise SystemExit("no ROWCSV lines found — this log predates the log-recovery change")
    return list(csv.DictReader(io.StringIO("\n".join(rows))))

if __name__ == "__main__":
    src, dst = sys.argv[1], sys.argv[2]
    recs = parse(open(src).read())
    with open(dst, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(recs[0].keys()))
        w.writeheader(); w.writerows(recs)
    nf = sum(1 for r in recs if (r.get("error") or "").strip())
    print(f"recovered {len(recs)} rows ({len(recs)-nf} ok, {nf} failed) -> {dst}")
