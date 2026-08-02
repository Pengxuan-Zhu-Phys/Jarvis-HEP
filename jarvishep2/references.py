#!/usr/bin/env python3
"""Citable references for Jarvis2 and its bundled samplers."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Reference:
    """One human-readable literature reference."""

    label: str
    title: str
    authors: str
    arxiv: str | None = None
    doi: str | None = None
    url: str | None = None

    def format(self, number: int) -> str:
        lines = [f"[{number}] {self.label}", f"    {self.title}", f"    {self.authors}"]
        if self.arxiv:
            lines.append(f"    arXiv: {self.arxiv}")
        if self.doi:
            lines.append(f"    DOI: {self.doi}")
        if self.url:
            lines.append(f"    URL: {self.url}")
        return "\n".join(lines)


FRAMEWORK_REFERENCES: tuple[Reference, ...] = (
    Reference(
        label="Jarvis-HEP framework",
        title="Jarvis-HEP: A lightweight Python framework for workflow composition and parameter scans in high-energy physics",
        authors="Erdong Guo, Paul Jackson, Jin-Min Yang, Pengxuan Zhu",
        arxiv="2604.25557",
    ),
)

SAMPLER_REFERENCES: tuple[Reference, ...] = (
    Reference(
        label="Bridson / AdaptiveBridson",
        title="Fast Poisson Disk Sampling in Arbitrary Dimensions",
        authors="Robert Bridson",
        doi="10.1145/1278780.1278807",
        url="https://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf",
    ),
    Reference(
        label="MultiNest",
        title="MULTINEST: an efficient and robust Bayesian inference tool for cosmology and particle physics",
        authors="F. Feroz, M. P. Hobson, M. Bridges",
        arxiv="0809.3437",
        doi="10.1111/j.1365-2966.2009.14548.x",
    ),
    Reference(
        label="Nested sampling",
        title="Nested Sampling for General Bayesian Computation",
        authors="John Skilling",
        doi="10.1214/06-BA127",
    ),
    Reference(
        label="Dynesty",
        title="dynesty: A Dynamic Nested Sampling Package for Estimating Bayesian Posteriors and Evidences",
        authors="Joshua S. Speagle",
        arxiv="1904.02180",
        doi="10.1093/mnras/staa278",
    ),
    Reference(
        label="MCMC family",
        title="Equation of State Calculations by Fast Computing Machines",
        authors="Nicholas Metropolis et al.",
        doi="10.1063/1.1699114",
    ),
    Reference(
        label="MCMC family",
        title="Monte Carlo Sampling Methods Using Markov Chains and Their Applications",
        authors="W. K. Hastings",
        doi="10.1093/biomet/57.1.97",
    ),
)


def render_references() -> str:
    """Return the stable text emitted by ``Jarvis2 refs``.

    Random, Grid, and CSV are native proposal/input strategies rather than
    bundled external algorithms, so they deliberately have no separate paper
    in this list.
    """
    sections = ["Jarvis-HEP V2 references", "", "Framework"]
    number = 1
    for reference in FRAMEWORK_REFERENCES:
        sections.extend((reference.format(number), ""))
        number += 1
    sections.append("Built-in sampler methods")
    for reference in SAMPLER_REFERENCES:
        sections.extend((reference.format(number), ""))
        number += 1
    sections.append(
        "Random, Grid, and CSV are native Jarvis2 proposal/input strategies; no external sampler citation is required."
    )
    return "\n".join(sections).rstrip() + "\n"


__all__ = ["FRAMEWORK_REFERENCES", "SAMPLER_REFERENCES", "Reference", "render_references"]
