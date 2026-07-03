"""Dataset discovery.

Resolves a *data root* into a :class:`~stride_analysis.models.Dataset` spanning
both data levels, without hardcoding any absolute path or local layout. Two
layouts are supported and auto-detected:

**Nested layout** (real STRIDE tree)::

    <root>/
      1st_run/DENV1/analysis_output/DENV1_correlations_v5.csv
      2nd_run/DENV1/analysis_output/DENV1_correlations_v5.csv
      3rd_run/DENV1/analysis_output/DENV1_correlations_v5.csv
      summaries/DENV1_profile.csv          (or profiles alongside root)
      summaries/DENV1_mechanism.json

**Flat layout** (convenient for examples/tests)::

    <root>/
      replicates/<run>/<serotype>/analysis_output/<serotype>_correlations_v5.csv
      summaries/<serotype>_profile.csv
      summaries/<serotype>_mechanism.json

Discovery is permissive about *where* the summaries live (root, ``summaries/``,
or next to the replicates) but strict about *consistency*: it fails loudly on
missing components, inconsistent replicate counts, or duplicate serotypes.
"""
from __future__ import annotations

from pathlib import Path

from ..models import (
    Dataset,
    ReplicateInput,
    SerotypeDataset,
    SummaryInput,
)
from ..models.errors import DiscoveryError
from ..models.schema import (
    ANALYSIS_OUTPUT_DIRNAME,
    CORRELATIONS_SUFFIX,
    MECHANISM_SUFFIX,
    PROFILE_SUFFIX,
)


def _run_dirs(root: Path) -> list[Path]:
    """Return candidate run directories under ``root`` (nested or ./replicates)."""
    replicates_root = root / "replicates"
    base = replicates_root if replicates_root.is_dir() else root
    runs = [
        d
        for d in sorted(base.iterdir())
        if d.is_dir() and _has_analysis_children(d)
    ]
    return runs


def _has_analysis_children(run_dir: Path) -> bool:
    """A run dir looks like a run dir if any child has an analysis_output/."""
    for child in run_dir.iterdir():
        if child.is_dir() and (child / ANALYSIS_OUTPUT_DIRNAME).is_dir():
            return True
    return False


def _find_correlations(analysis_dir: Path, serotype: str) -> Path | None:
    """Locate the replicate table inside an analysis_output directory."""
    exact = analysis_dir / f"{serotype}{CORRELATIONS_SUFFIX}"
    if exact.is_file():
        return exact
    # fall back to any *_correlations_v5.csv (single-file convention)
    hits = sorted(analysis_dir.glob(f"*{CORRELATIONS_SUFFIX}"))
    return hits[0] if hits else None


def _summary_search_dirs(root: Path) -> list[Path]:
    """Directories to search for Level-2 summaries, in priority order."""
    candidates = [root / "summaries", root]
    return [d for d in candidates if d.is_dir()]


def discover_replicates(root: Path) -> dict[str, list[ReplicateInput]]:
    """Discover Level-1 replicate tables, grouped by serotype.

    Replicate index is 1-based in sorted run-dir order (so ``1st_run`` → 1).
    """
    runs = _run_dirs(root)
    by_serotype: dict[str, list[ReplicateInput]] = {}
    for idx, run_dir in enumerate(runs, start=1):
        for sero_dir in sorted(run_dir.iterdir()):
            if not sero_dir.is_dir():
                continue
            analysis_dir = sero_dir / ANALYSIS_OUTPUT_DIRNAME
            if not analysis_dir.is_dir():
                continue
            serotype = sero_dir.name
            corr = _find_correlations(analysis_dir, serotype)
            if corr is None:
                raise DiscoveryError(
                    f"run {run_dir.name!r} serotype {serotype!r}: no "
                    f"'*{CORRELATIONS_SUFFIX}' in {analysis_dir}"
                )
            by_serotype.setdefault(serotype, []).append(
                ReplicateInput(
                    serotype=serotype,
                    run_dir=run_dir.name,
                    replicate_index=idx,
                    correlations_path=corr,
                )
            )
    return by_serotype


def discover_summaries(root: Path) -> dict[str, SummaryInput]:
    """Discover Level-2 profile/mechanism pairs, grouped by serotype.

    Fails on any profile without a mechanism (or vice versa) or on a duplicate
    profile for the same serotype across search directories.
    """
    profiles: dict[str, Path] = {}
    mechanisms: dict[str, Path] = {}
    for sdir in _summary_search_dirs(root):
        for p in sorted(sdir.glob(f"*{PROFILE_SUFFIX}")):
            sero = p.name[: -len(PROFILE_SUFFIX)]
            if sero in profiles:
                raise DiscoveryError(
                    f"duplicate profile for serotype {sero!r}: "
                    f"{profiles[sero]} and {p}"
                )
            profiles[sero] = p
        for p in sorted(sdir.glob(f"*{MECHANISM_SUFFIX}")):
            sero = p.name[: -len(MECHANISM_SUFFIX)]
            if sero in mechanisms:
                raise DiscoveryError(
                    f"duplicate mechanism for serotype {sero!r}: "
                    f"{mechanisms[sero]} and {p}"
                )
            mechanisms[sero] = p

    prof_only = sorted(set(profiles) - set(mechanisms))
    mech_only = sorted(set(mechanisms) - set(profiles))
    problems = []
    if prof_only:
        problems.append(f"profiles without mechanism: {prof_only}")
    if mech_only:
        problems.append(f"mechanisms without profile: {mech_only}")
    if problems:
        raise DiscoveryError("unmatched Level-2 summaries: " + "; ".join(problems))

    return {
        sero: SummaryInput(
            serotype=sero,
            profile_path=profiles[sero],
            mechanism_path=mechanisms[sero],
        )
        for sero in sorted(profiles)
    }


def discover_dataset(
    root: str | Path,
    *,
    require_replicates: bool = True,
    require_summaries: bool = True,
    enforce_equal_replicate_counts: bool = True,
) -> Dataset:
    """Discover a full dataset under ``root``.

    Parameters
    ----------
    require_replicates, require_summaries
        If True, every serotype must provide that data level; otherwise the
        level is optional (partial datasets are allowed).
    enforce_equal_replicate_counts
        If True (default), all serotypes must have the same replicate count;
        a mismatch raises :class:`DiscoveryError`.

    Raises :class:`DiscoveryError` on any structural problem.
    """
    root_path = Path(root)
    if not root_path.is_dir():
        raise DiscoveryError(f"data root does not exist: {root_path}")

    replicates = discover_replicates(root_path)
    summaries = discover_summaries(root_path)

    all_serotypes = sorted(set(replicates) | set(summaries))
    if not all_serotypes:
        raise DiscoveryError(
            f"no STRIDE data found under {root_path} "
            f"(expected replicate '*{CORRELATIONS_SUFFIX}' and/or "
            f"Level-2 '*{PROFILE_SUFFIX}' files)"
        )

    # requirement checks -----------------------------------------------------
    if require_replicates:
        missing = [s for s in all_serotypes if s not in replicates]
        if missing:
            raise DiscoveryError(
                f"serotype(s) missing replicate data: {missing}"
            )
    if require_summaries:
        missing = [s for s in all_serotypes if s not in summaries]
        if missing:
            raise DiscoveryError(
                f"serotype(s) missing Level-2 summaries: {missing}"
            )

    # equal replicate counts -------------------------------------------------
    if enforce_equal_replicate_counts and replicates:
        counts = {s: len(replicates.get(s, [])) for s in replicates}
        distinct = set(counts.values())
        if len(distinct) > 1:
            raise DiscoveryError(
                f"inconsistent replicate counts across serotypes: {counts}"
            )

    serotype_datasets = tuple(
        SerotypeDataset(
            serotype=s,
            replicates=tuple(replicates.get(s, ())),
            summary=summaries.get(s),
        )
        for s in all_serotypes
    )
    return Dataset(serotypes=serotype_datasets)
