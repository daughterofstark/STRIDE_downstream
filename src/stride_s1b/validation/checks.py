"""Structural validation of the S1B annotation layer.

Structural checks only — no statistical assertions. Each check appends a
:class:`~stride_s1b.models.ValidationCheck` to the report and raises
:class:`~stride_s1b.models.errors.ConsistencyError` on failure, so later stages
can trust the annotation layer.

Checks:

- every residue has exactly one biological annotation;
- hierarchy paths are unique (one residue per (serotype, hierarchy_path));
- every domain membership is internally consistent (residue counts agree, and
  every annotated residue's (serotype, chain, domain) exists in the domain
  annotation);
- no orphan annotations (every annotated residue exists in the S1A canonical
  residues);
- every serotype references existing canonical residues;
- referential integrity between every generated table.
"""
from __future__ import annotations

import pandas as pd

from ..models import S1BReport
from ..models.errors import ConsistencyError


def validate_one_annotation_per_residue(
    residue_annotation: pd.DataFrame,
    canonical_residues: pd.DataFrame,
    report: S1BReport,
) -> None:
    """Every residue has exactly one biological annotation, and no orphans.

    The residue annotation key set must equal the S1A canonical residue key set
    (a bijection): no residue annotated twice, none missing, none orphaned.
    """
    ann_keys_list = list(
        zip(
            residue_annotation["serotype"],
            residue_annotation["canon_label"],
            strict=True,
        )
    )
    ann_keys = set(ann_keys_list)
    if len(ann_keys_list) != len(ann_keys):
        raise ConsistencyError(
            "residue annotation has duplicate (serotype, canon_label) rows"
        )
    canon_keys = set(
        zip(
            canonical_residues["serotype"],
            canonical_residues["canon_label"],
            strict=True,
        )
    )
    missing = canon_keys - ann_keys  # residue without an annotation
    orphan = ann_keys - canon_keys   # annotation without a residue
    if missing:
        raise ConsistencyError(
            f"{len(missing)} canonical residue(s) have no annotation, e.g. "
            f"{sorted(missing)[:3]}"
        )
    if orphan:
        raise ConsistencyError(
            f"{len(orphan)} annotation(s) reference no canonical residue "
            f"(orphan), e.g. {sorted(orphan)[:3]}"
        )
    report.add(
        "every residue has exactly one biological annotation",
        "global",
        True,
        f"{len(ann_keys)} residues annotated 1:1",
    )


def validate_unique_hierarchy_paths(
    hierarchy_annotation: pd.DataFrame, report: S1BReport
) -> None:
    """Hierarchy paths remain unique: one row per (serotype, hierarchy_path)."""
    if hierarchy_annotation.empty:
        report.add("hierarchy paths unique", "global", True, "empty")
        return
    dup = hierarchy_annotation.duplicated(
        ["serotype", "hierarchy_path"], keep=False
    )
    if dup.any():
        examples = (
            hierarchy_annotation.loc[dup, ["serotype", "hierarchy_path"]]
            .drop_duplicates()
            .head(3)
            .to_dict("records")
        )
        raise ConsistencyError(
            f"hierarchy annotation key not unique; examples: {examples}"
        )
    report.add(
        "hierarchy paths remain unique",
        "global",
        True,
        f"{len(hierarchy_annotation)} unique paths",
    )


def validate_domain_membership(
    domain_annotation: pd.DataFrame,
    residue_annotation: pd.DataFrame,
    report: S1BReport,
) -> None:
    """Every domain membership is internally consistent.

    - the domain annotation key (serotype, chain, domain) is unique;
    - each annotated residue's (serotype, chain, domain) exists in the domain
      annotation (no residue belongs to an unknown domain);
    - each domain's ``n_residues`` equals the number of residue-annotation rows
      that map to it.
    """
    dup = domain_annotation.duplicated(["serotype", "chain", "domain"], keep=False)
    if dup.any():
        raise ConsistencyError(
            "domain annotation key (serotype, chain, domain) not unique"
        )

    domain_keys = set(
        zip(
            domain_annotation["serotype"],
            domain_annotation["chain"],
            domain_annotation["domain"],
            strict=True,
        )
    )
    if not residue_annotation.empty:
        res_domain_keys = set(
            zip(
                residue_annotation["serotype"],
                residue_annotation["chain"],
                residue_annotation["domain"],
                strict=True,
            )
        )
        unknown = res_domain_keys - domain_keys
        if unknown:
            raise ConsistencyError(
                f"{len(unknown)} residue domain membership(s) reference a domain "
                f"absent from the domain annotation, e.g. {sorted(unknown)[:3]}"
            )

        # counts must agree
        counts = (
            residue_annotation.groupby(["serotype", "chain", "domain"])
            .size()
            .to_dict()
        )
        for row in domain_annotation.itertuples(index=False):
            key = (row.serotype, row.chain, row.domain)
            expected = int(counts.get(key, 0))
            if int(row.n_residues) != expected:
                raise ConsistencyError(
                    f"domain {key} declares n_residues={row.n_residues} but "
                    f"{expected} residue annotation(s) map to it"
                )
    report.add(
        "every domain membership is internally consistent",
        "global",
        True,
        f"{len(domain_keys)} domains, counts agree with residue annotation",
    )


def validate_serotype_references(
    serotype_annotation: pd.DataFrame,
    canonical_residues: pd.DataFrame,
    residue_annotation: pd.DataFrame,
    report: S1BReport,
) -> None:
    """Every serotype references existing canonical residues.

    - the serotype annotation key is unique;
    - every serotype in the annotation exists in the S1A canonical residues;
    - each serotype's ``n_residues`` equals its residue-annotation row count.
    """
    if serotype_annotation.duplicated(["serotype"]).any():
        raise ConsistencyError("serotype annotation has duplicate serotype rows")

    annot_serotypes = set(serotype_annotation["serotype"])
    canon_serotypes = set(canonical_residues["serotype"])
    unknown = annot_serotypes - canon_serotypes
    if unknown:
        raise ConsistencyError(
            f"serotype annotation references serotype(s) with no canonical "
            f"residues: {sorted(unknown)}"
        )

    res_counts = (
        residue_annotation.groupby("serotype").size().to_dict()
        if not residue_annotation.empty
        else {}
    )
    for row in serotype_annotation.itertuples(index=False):
        expected = int(res_counts.get(row.serotype, 0))
        if int(row.n_residues) != expected:
            raise ConsistencyError(
                f"serotype {row.serotype!r} declares n_residues={row.n_residues} "
                f"but {expected} residue annotation(s) exist"
            )
    report.add(
        "every serotype references existing canonical residues",
        "global",
        True,
        f"{len(annot_serotypes)} serotypes",
    )


def validate_referential_integrity(
    residue_annotation: pd.DataFrame,
    domain_annotation: pd.DataFrame,
    hierarchy_annotation: pd.DataFrame,
    serotype_annotation: pd.DataFrame,
    report: S1BReport,
) -> None:
    """Referential integrity between every generated table.

    - every residue annotation has a matching hierarchy annotation row
      (by (serotype, hierarchy_path)) and vice versa;
    - every (serotype, chain, domain) in the residue annotation appears in the
      domain annotation;
    - every serotype appearing in any table appears in the serotype annotation.
    """
    if residue_annotation.empty:
        report.add("referential integrity between tables", "global", True, "empty")
        return

    res_paths = set(
        zip(
            residue_annotation["serotype"],
            residue_annotation["hierarchy_path"],
            strict=True,
        )
    )
    hier_paths = set(
        zip(
            hierarchy_annotation["serotype"],
            hierarchy_annotation["hierarchy_path"],
            strict=True,
        )
    )
    if res_paths != hier_paths:
        only_res = sorted(res_paths - hier_paths)[:3]
        only_hier = sorted(hier_paths - res_paths)[:3]
        raise ConsistencyError(
            "residue and hierarchy annotations disagree on (serotype, "
            f"hierarchy_path); residue-only={only_res}, hierarchy-only={only_hier}"
        )

    res_domain_keys = set(
        zip(
            residue_annotation["serotype"],
            residue_annotation["chain"],
            residue_annotation["domain"],
            strict=True,
        )
    )
    domain_keys = set(
        zip(
            domain_annotation["serotype"],
            domain_annotation["chain"],
            domain_annotation["domain"],
            strict=True,
        )
    )
    if not res_domain_keys <= domain_keys:
        missing = sorted(res_domain_keys - domain_keys)[:3]
        raise ConsistencyError(
            f"residue annotation references domains absent from the domain "
            f"annotation, e.g. {missing}"
        )

    all_serotypes = (
        set(residue_annotation["serotype"])
        | set(domain_annotation["serotype"])
        | set(hierarchy_annotation["serotype"])
    )
    annot_serotypes = set(serotype_annotation["serotype"])
    if all_serotypes != annot_serotypes:
        raise ConsistencyError(
            "serotype annotation does not cover exactly the serotypes present in "
            f"the other tables; tables={sorted(all_serotypes)}, "
            f"serotype_annotation={sorted(annot_serotypes)}"
        )
    report.add(
        "referential integrity between every generated table",
        "global",
        True,
        f"{len(res_paths)} residues cross-linked across 4 tables",
    )
