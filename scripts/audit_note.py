"""
Check every quantitative claim in note/note.tex against a saved log.

An earlier version pasted whole log files into one blob and asked whether each
number in the note appeared somewhere in it. That is too weak: short strings such
as "10" or "14" match incidentally, and it let the note claim agreement "to
1e-14" when the worst measured value was 9.37e-13, off by more than an order of
magnitude.

This version is a whitelist with arithmetic. Every numeric quantity the note
asserts is declared below together with the value in the log it comes from and
the relation it should satisfy:

  exact   the note quotes the logged value verbatim
  round   the note rounds it, and the rounding is checked numerically
  bound   the note claims agreement better than some figure, and the logged
          value must actually be smaller

Anything numeric in the note that is neither declared here nor structural fails
the audit.
"""
import os, re, sys

D = "out/verify"

TEXT = [
    ("s2^-1 s1 s2^-1 s1", "headline_reverify_2026-07-29.log", "vineyard braid word"),
    ("t - 3 + 1/t",       "alexander_2026-07-29.log",         "Alexander polynomial"),
    ("[1, -3, 1]",        "alexander_2026-07-29.log",         "Alexander coefficients"),
    ("0.33333",           "rigor.log",                        "T/3 set-return"),
    ("0.66667",           "rigor.log",                        "2T/3 set-return"),
    ("order 3",           "full_period_corrected.log",        "order over T/3"),
    ("order 1",           "full_period_corrected.log",        "order over T"),
    ("isosceles configurations in one period: 12",
                          "what_it_means.log",                "twelve walls"),
    ("HOLDS",             "wall_proof.log",                   "hypothesis of the proposition"),
    ("(C) holds and the conclusion holds",
                          "concavity_lemma.log",              "concavity closes the argument"),
    ("858 of 858",        "deaths_and_symmetry.log",          "deaths are the MST"),
    ("vines MEET",        "wall_classification.log",          "six vine collisions"),
    ("6.7e-16",           "frame_test.log",                   "distances across frames"),
    ("0.00e+00",          "crosscheck_gudhi.log",             "gudhi circle agreement"),
    ("m004",              "knot_id_independent.log",          "snappy identifies 4_1"),
    ("ORDER of the births equals ORDER of the opposite sides: 858 of 858",
                          "completeness_test.log",            "strand order is the side order"),
    ("cauchy s=0.35",     "kernel_diagrams.log",              "Cauchy kernel, real diagram"),
    ("exponential cusp",  "kernel_diagrams.log",              "cusped exponential kernel"),
    ("student t3",        "kernel_diagrams.log",              "Student t kernel"),
    ("conjugate to s1 s2^-1 in 18 of 18 projections",
                          "trajectory_braid.log",             "trajectory braid class"),
    ("conjugate as braids?      yes, by s1",
                          "trajectory_braid.log",             "vineyard is the square"),
    ("= s2^-1 s1 s2^-1 s1", "trajectory_braid.log",           "conjugation returns the word"),
    ("once-punctured-torus bundle",
                          "knot_meaning_and_N.log",           "the bundle interpretation"),
    ("ratio 2.000000",    "trajectory_braid.log",             "entropy ratio is exactly two"),
]

NUMERIC = [
    ("1.924847", "1.924847", "knot_meaning_and_N.log",  "exact",
     "vineyard topological entropy"),
    ("0.962424", "0.962424", "knot_meaning_and_N.log",  "exact",
     "trajectory topological entropy"),
    ("0.06",     "0.0598",   "sigma_window.log",        "round",
     "least persistent feature over the loop"),
    ("0.24",     "0.24",     "sigma_window.log",        "exact",
     "lower end of the bandwidth window"),
    ("0.30",     "0.30",     "sigma_window.log",        "exact",
     "upper end of the bandwidth window"),
]

# indices, exponents, subscripts, and counts the note states in words
STRUCTURAL = {"0", "1", "2", "3", "4", "10", "12", "13", "23"}


def note_numbers():
    tex = open("note/note.tex").read()
    tex = re.sub(r"\\(documentclass|usepackage|setlength|linespread|pagestyle)[^\n]*", "", tex)
    tex = re.sub(r"\[[^\]]*(mm|em|pt)[^\]]*\]", "", tex)
    tex = re.sub(r"\\vspace\{[^}]*\}", "", tex)
    return sorted(set(re.findall(r"\d+\.\d+(?:e[-+]?\d+)?|\d+", tex)))


def main():
    ok = bad = 0
    print("text claims")
    print("-" * 78)
    for needle, log, what in TEXT:
        p = os.path.join(D, log)
        if os.path.exists(p) and needle in open(p).read():
            print(f"  ok    {what:40s} {log}"); ok += 1
        else:
            print(f"  FAIL  {what:40s} {log}"); bad += 1

    print("\nnumeric claims, arithmetic checked")
    print("-" * 78)
    accounted = set(STRUCTURAL)
    for shown, actual, log, rel, what in NUMERIC:
        p = os.path.join(D, log)
        if not os.path.exists(p) or actual not in open(p).read():
            print(f"  FAIL  {what:40s} {actual} not in {log}"); bad += 1; continue
        a, b = float(shown), float(actual)
        good = (a == b if rel == "exact"
                else b <= a if rel == "bound"
                else abs(a - b) <= 0.02 * abs(b))
        print(f"  {'ok  ' if good else 'FAIL'}  {what:40s} note {shown} {rel} log {actual}")
        if good: ok += 1
        else: bad += 1
        for tok in re.findall(r"\d+\.\d+|\d+", shown):
            accounted.add(tok)

    loose = [n for n in note_numbers() if n not in accounted]
    print("\n" + "-" * 78)
    print(f"{ok} verified, {bad} failed")
    print(f"undeclared numbers in the note: {loose if loose else 'none'}")
    sys.exit(1 if (bad or loose) else 0)


if __name__ == "__main__":
    main()
