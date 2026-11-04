#!/usr/bin/env python2.7
"""Cohort-wide Conpair: all-vs-all concordance and per-sample contamination.

Runs inside the Conpair container (docker://cmopipeline/conpair:v0.3.3) with
its Python 2.7, so type hints are written as comments. Genotype likelihoods
and contamination estimates come from Conpair's own code in the container,
imported at run time. Nothing from Conpair is copied into this repository:
the NYGC license is single-user and does not allow redistribution.

How this differs from running Conpair one pair at a time:

  - Genotype likelihoods are computed once per sample, then concordance is
    computed for every ordered pair (reference, query). The per-pair rule is
    the one in verify_concordance(s).py, so the values are the same.
  - Contamination is computed by the function in
    estimate_tumor_normal_contaminations.py. A sample's "single" estimate
    is Conpair's normal-sample estimate (its own pileup against population
    allele frequencies). Its "paired" estimate is Conpair's tumor-sample
    estimate, using the homozygous markers of the best-matching other
    sample in the cohort (its partner) as the known genotype.

Writes into --out:
  conpair_samples.csv        marker coverage and genotype summary per sample
  conpair_concordance.csv    every ordered pair (reference, query)
  conpair_contamination.csv  single and paired estimate per sample
  conpair_run.csv            version, parameters and timing
"""

from __future__ import division, print_function

import argparse
import csv
import multiprocessing
import os
import sys
import tempfile
import time

import numpy as np

## The container's Python 2.7 has no typing module; the names are only
## needed by type checkers reading the type comments.
try:
    from typing import Any, Dict, List, Optional, Tuple  # noqa: F401
except ImportError:
    pass

CONPAIR_DIR = os.environ.get("CONPAIR_DIR", "/usr/bin/conpair")
sys.path.insert(0, os.path.join(CONPAIR_DIR, "modules"))
sys.path.insert(0, os.path.join(CONPAIR_DIR, "scripts"))

## Conpair's own code, present only inside the container.
import ContaminationMarker  # type: ignore  # noqa: E402
import estimate_tumor_normal_contaminations as conpair_contam  # type: ignore  # noqa: E402

TOOL_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

## Genotype codes: index of the largest Conpair likelihood, as in
## verify_concordance.py. NOCALL marks a marker with no ref or alt reads.
AB = 1
NOCALL = -1

## Set in main() before the worker pool forks, so workers inherit it.
MARKERS = None  # type: Any


def log(message):
    # type: (str) -> None
    """Write a time-stamped progress line to stderr.

    Args:
        message: Text to write.
    """
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    sys.stderr.write("[{0}] {1}\n".format(stamp, message))
    sys.stderr.flush()


def read_version():
    # type: () -> str
    """Read the tool version from conpair/VERSION.

    Returns:
        The version string.
    """
    with open(os.path.join(TOOL_DIR, "VERSION")) as fp:
        return fp.read().strip()


def read_samples(manifest):
    # type: (str) -> List[str]
    """Read sample IDs, in manifest order, from a tab-separated manifest.

    Args:
        manifest: Path to a TSV with a header containing a "sample" column.

    Returns:
        Sample IDs in file order.

    Raises:
        ValueError: If the "sample" column is missing or an ID repeats.
    """
    with open(manifest) as fp:
        header = fp.readline().rstrip("\r\n").split("\t")
        if "sample" not in header:
            raise ValueError("manifest has no 'sample' column: " + manifest)
        col = header.index("sample")
        samples = [line.rstrip("\r\n").split("\t")[col]
                   for line in fp if line.strip()]
    dups = sorted(set(s for s in samples if samples.count(s) > 1))
    if dups:
        raise ValueError("duplicate sample IDs: " + ", ".join(dups))
    return samples


def genotype_calls(job):
    # type: (Tuple[str, str, int, int]) -> Tuple[str, np.ndarray, np.ndarray]
    """Call genotypes for one sample with Conpair's likelihood model.

    Args:
        job: (sample, pileup path, min mapping quality, min base quality).

    Returns:
        (sample, calls, coverage), both arrays in marker-file order. A call
        is the index of the largest Conpair likelihood (0 AA, 1 AB, 2 BB),
        or NOCALL where the marker has no ref or alt reads. Coverage is
        Conpair's depth after the quality filters, 0 at NOCALL markers.
    """
    sample, pileup, mmq, mbq = job
    likelihoods = ContaminationMarker.genotype_likelihoods_for_markers(
        MARKERS, pileup, min_map_quality=mmq, min_base_quality=mbq)
    calls = np.full(len(MARKERS), NOCALL, dtype=np.int8)
    coverage = np.zeros(len(MARKERS), dtype=np.int32)
    for i, marker in enumerate(MARKERS):
        value = likelihoods[marker]
        if value is None:
            continue
        lk = value["likelihoods"]
        calls[i] = lk.index(max(lk))
        coverage[i] = value["coverage"]
    return sample, calls, coverage


def percent(hits, total):
    # type: (int, int) -> Optional[float]
    """Percentage rounded to 2 places as Conpair reports it.

    Args:
        hits: Numerator.
        total: Denominator.

    Returns:
        100 * hits / total rounded to 2 places, or None if total is 0.
    """
    if total == 0:
        return None
    return round(100.0 * hits / total, 2)


def concordance_rows(samples, calls, coverage, min_cov):
    # type: (List[str], np.ndarray, np.ndarray, int) -> List[tuple]
    """Concordance for every ordered pair of samples.

    Follows verify_concordance.py: a marker counts when both samples have
    ref or alt reads and coverage >= min_cov; the pair is concordant there
    when the called genotypes are equal. The "hom" value additionally keeps
    only markers called homozygous in the reference sample, which Conpair
    uses (--normal_homozygous_markers_only) so that LOH and copy number
    changes in the query do not lower concordance.

    Args:
        samples: Sample IDs, the row order of calls and coverage.
        calls: Genotype calls, shape (samples, markers).
        coverage: Coverage, shape (samples, markers).
        min_cov: Minimum coverage in both samples for a marker to count.

    Returns:
        Tuples (reference, query, concordanceHom, nMarkersHom,
        concordanceAll, nMarkersAll), reference != query.
    """
    usable = (calls != NOCALL) & (coverage >= min_cov)
    hom = usable & (calls != AB)
    rows = []
    for r in range(len(samples)):
        for q in range(len(samples)):
            if r == q:
                continue
            both = usable[r] & usable[q]
            same = both & (calls[r] == calls[q])
            n_all = int(both.sum())
            n_hom = int((both & hom[r]).sum())
            rows.append((samples[r], samples[q],
                         percent(int((same & hom[r]).sum()), n_hom), n_hom,
                         percent(int(same.sum()), n_all), n_all))
    return rows


def choose_partners(rows, match, min_markers):
    # type: (List[tuple], float, int) -> Dict[str, Tuple[str, float]]
    """Pick, for each sample, the other sample whose genotype it best matches.

    A candidate partner P for query Q needs concordanceHom(P -> Q) >= match
    over at least min_markers markers. Using P's homozygous markers is what
    Conpair's tumor contamination estimate needs: a marker homozygous in P
    is homozygous in Q when both come from the same individual, even if Q
    has LOH. Among candidates the highest concordance wins, then the most
    markers, then manifest order.

    Args:
        rows: Output of concordance_rows.
        match: Minimum concordanceHom (percent) to count as the same individual.
        min_markers: Minimum markers behind the concordance.

    Returns:
        query -> (partner, concordanceHom(partner -> query)).
    """
    best = {}
    for ref, query, conc_hom, n_hom, _, _ in rows:
        if conc_hom is None or n_hom < min_markers or conc_hom < match:
            continue
        if query not in best or (conc_hom, n_hom) > best[query][1]:
            best[query] = (ref, (conc_hom, n_hom))
    return dict((q, (v[0], v[1][0])) for q, v in best.items())


def parse_estimate(field):
    # type: (str) -> Optional[float]
    """Parse one field of Conpair's contamination result string.

    Args:
        field: "Normal:<pct>", "Tumor:<pct>" or "WARNING:NA".

    Returns:
        The percentage, or None for "WARNING:NA".
    """
    label, value = field.split(":", 1)
    if label == "WARNING" or value == "NA":
        return None
    return float(value)


def estimate_contamination(job):
    # type: (Tuple[str, Optional[str], str, str, str, float, int]) -> tuple
    """Run Conpair's contamination estimate for one (reference, query) job.

    Args:
        job: (reference, query or None, reference pileup, query pileup,
            marker file, grid step, min mapping quality). A None query uses
            an empty pileup, so only the reference's single estimate is made.

    Returns:
        (reference, query, single estimate of reference, paired estimate of
        query, seconds taken). Estimates are percentages or None.
    """
    ref, query, ref_pileup, query_pileup, markers_file, grid, mmq = job
    start = time.time()
    result = conpair_contam.run_estimate_tumor_normal_contamination(
        ref_pileup, query_pileup, markers_file, grid, mmq)
    normal_field, tumor_field = result.split("\t")
    return (ref, query, parse_estimate(normal_field),
            parse_estimate(tumor_field), time.time() - start)


def fmt(value):
    # type: (object) -> str
    """Format a value for TSV output, None as NA.

    Args:
        value: Any value.

    Returns:
        The string form, or "NA" for None.
    """
    return "NA" if value is None else str(value)


def write_csv(path, header, rows):
    # type: (str, List[str], List[tuple]) -> None
    """Write rows to a comma-separated file with a header line.

    Fields are quoted only when they contain a comma, a quote or a line
    break.

    Args:
        path: Output file.
        header: Column names.
        rows: Row tuples, None written as NA.
    """
    with open(path, "w") as fp:
        writer = csv.writer(fp, lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow([fmt(v) for v in row])


def parse_args(argv):
    # type: (List[str]) -> argparse.Namespace
    """Parse command line arguments.

    Args:
        argv: Arguments without the program name.

    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        prog="conpairCohort.py",
        description="All-vs-all Conpair concordance and per-sample "
                    "contamination for a cohort of pileups.")
    parser.add_argument("--manifest", required=True,
                        help="TSV with a 'sample' column")
    parser.add_argument("--pileup-dir", required=True,
                        help="directory holding <sample>.pileup")
    parser.add_argument("--markers", required=True,
                        help="Conpair marker .txt file")
    parser.add_argument("--out", required=True, help="output directory")
    parser.add_argument("--min-cov", type=int, default=10,
                        help="min coverage to call a genotype [10]")
    parser.add_argument("--mmq", type=int, default=10,
                        help="min mapping quality [10]")
    parser.add_argument("--mbq", type=int, default=20,
                        help="min base quality, concordance only [20]")
    parser.add_argument("--grid", type=float, default=0.01,
                        help="contamination grid step [0.01]")
    parser.add_argument("--match", type=float, default=95.0,
                        help="min concordanceHom (percent) for a partner [95]")
    parser.add_argument("--min-markers", type=int, default=100,
                        help="min markers behind a concordance value [100]")
    parser.add_argument("--threads", type=int, default=1,
                        help="worker processes [1]")
    return parser.parse_args(argv)


def main(argv):
    # type: (List[str]) -> int
    """Run the cohort analysis.

    Args:
        argv: Arguments without the program name.

    Returns:
        Process exit status.
    """
    global MARKERS
    args = parse_args(argv)
    started = time.time()
    version = read_version()
    log("conpairCohort.py {0}".format(version))

    samples = read_samples(args.manifest)
    pileups = dict((s, os.path.join(args.pileup_dir, s + ".pileup"))
                   for s in samples)
    missing = [s for s in samples if not os.path.isfile(pileups[s])]
    if missing:
        log("ERROR: missing pileups: " + ", ".join(missing))
        return 1

    MARKERS = ContaminationMarker.get_markers(args.markers)
    log("{0} samples, {1} markers, {2} threads".format(
        len(samples), len(MARKERS), args.threads))
    pool = multiprocessing.Pool(args.threads)

    ##
    ## Genotype calls, once per sample.
    ##
    jobs = [(s, pileups[s], args.mmq, args.mbq) for s in samples]
    called = dict((s, (c, v)) for s, c, v in pool.imap_unordered(genotype_calls, jobs))
    calls = np.vstack([called[s][0] for s in samples])
    coverage = np.vstack([called[s][1] for s in samples])
    log("genotypes called")

    usable = (calls != NOCALL) & (coverage >= args.min_cov)
    sample_rows = []
    for i, s in enumerate(samples):
        n_usable = int(usable[i].sum())
        n_het = int((usable[i] & (calls[i] == AB)).sum())
        sample_rows.append((s, int((calls[i] != NOCALL).sum()), n_usable,
                            percent(n_het, n_usable),
                            float(np.median(coverage[i]))))
    write_csv(os.path.join(args.out, "conpair_samples.csv"),
              ["sample", "markersWithReads", "markersCovered", "pctHet",
               "medianCoverage"],
              sample_rows)

    ##
    ## Concordance, every ordered pair.
    ##
    conc_rows = concordance_rows(samples, calls, coverage, args.min_cov)
    write_csv(os.path.join(args.out, "conpair_concordance.csv"),
              ["reference", "query", "concordanceHom", "nMarkersHom",
               "concordanceAll", "nMarkersAll"],
              conc_rows)
    log("concordance done, {0} ordered pairs".format(len(conc_rows)))

    ##
    ## Contamination. One job per sample with a partner: Conpair's function
    ## returns the partner's single estimate and the sample's paired
    ## estimate together. Samples that are no one's partner get a job of
    ## their own against an empty pileup, for their single estimate.
    ##
    partners = choose_partners(conc_rows, args.match, args.min_markers)
    fd, empty_pileup = tempfile.mkstemp(suffix=".pileup", dir=args.out)
    os.close(fd)
    jobs = [(partners[q][0], q, pileups[partners[q][0]], pileups[q],
             args.markers, args.grid, args.mmq)
            for q in samples if q in partners]
    as_partner = set(p for p, _ in partners.values())
    jobs += [(s, None, pileups[s], empty_pileup, args.markers, args.grid, args.mmq)
             for s in samples if s not in as_partner]
    log("contamination: {0} jobs, {1} samples with a partner".format(
        len(jobs), len(partners)))

    single = {}
    paired = {}
    for k, (ref, query, est_single, est_paired, secs) in enumerate(
            pool.imap_unordered(estimate_contamination, jobs), 1):
        single[ref] = est_single
        if query is not None:
            paired[query] = est_paired
        log("contamination {0}/{1}: {2} -> {3} ({4:.0f} s)".format(
            k, len(jobs), ref, query, secs))
    pool.close()
    pool.join()
    os.remove(empty_pileup)

    contam_rows = []
    for s in samples:
        partner, partner_conc = partners.get(s, (None, None))
        contam_rows.append((s, single.get(s), partner, partner_conc,
                            paired.get(s)))
    write_csv(os.path.join(args.out, "conpair_contamination.csv"),
              ["sample", "contaminationSingle", "partner",
               "partnerConcordance", "contaminationPaired"],
              contam_rows)

    finished = time.time()
    run_rows = [
        ("version", version),
        ("conpairDir", CONPAIR_DIR),
        ("markers", args.markers),
        ("nSamples", len(samples)),
        ("nMarkers", len(MARKERS)),
        ("minCov", args.min_cov),
        ("mmq", args.mmq),
        ("mbq", args.mbq),
        ("grid", args.grid),
        ("match", args.match),
        ("minMarkers", args.min_markers),
        ("started", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(started))),
        ("seconds", int(finished - started)),
    ]
    write_csv(os.path.join(args.out, "conpair_run.csv"), ["key", "value"], run_rows)
    log("done in {0:.0f} s".format(finished - started))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
