"""
cefe_lab.tests -- a minimal framework for theory verification scripts.

A scientist subclasses TheoryTest, implements `measure()` returning a dict
of measured quantities, and declares `expectations`: tuples of
    (key, target, tolerance, label)
The framework produces pass/fail verdicts and a report.  A theory is
"confirmed" only when every expectation passes; a single decisive failure
falsifies it (as it should).

Example
-------
class AreaLaw(TheoryTest):
    name = "Srednicki area law"
    def measure(self):
        ...
        return {"alpha": 2.09, "kappa": 0.314}
    expectations = [("alpha", 2.0, 0.15, "area exponent"),
                    ("kappa", 0.30, 0.03, "Srednicki coefficient")]

print(AreaLaw().run().report())
"""

from dataclasses import dataclass, field
import numpy as np


@dataclass
class TestResult:
    key: str
    label: str
    measured: float
    target: float
    tolerance: float
    passed: bool

    def line(self):
        mark = "PASS" if self.passed else "FAIL"
        return (f"  [{mark}] {self.label}: measured {self.measured:.6g} "
                f"vs target {self.target:.6g} +/- {self.tolerance:.6g}")


@dataclass
class TheoryReport:
    name: str
    results: list = field(default_factory=list)

    @property
    def passed(self):
        return all(r.passed for r in self.results)

    def report(self):
        lines = [f"TheoryTest: {self.name}", "-" * 60]
        lines += [r.line() for r in self.results]
        verdict = "CONFIRMED (all expectations within tolerance)" if self.passed \
                  else "FALSIFIED or INCONCLUSIVE (see failed expectations)"
        lines += ["-" * 60, f"verdict: {verdict}"]
        return "\n".join(lines)


class TheoryTest:
    name = "unnamed theory test"
    expectations = []   # (key, target, tolerance, label)

    def measure(self):
        """Return a dict of measured floats.  Implemented by the user."""
        raise NotImplementedError

    def run(self):
        measured = self.measure()
        results = []
        for key, target, tol, label in self.expectations:
            if key not in measured:
                results.append(TestResult(key, label, float("nan"),
                                          target, tol, False))
                continue
            val = float(measured[key])
            passed = abs(val - target) <= tol
            results.append(TestResult(key, label, val, target, tol, passed))
        return TheoryReport(self.name, results)
