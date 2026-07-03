# Contributing

> **Placeholder.** Adapt before making the repository public.

Thanks for your interest in contributing to `stride-dengue-analysis`.

## Development setup

```bash
git clone https://github.com/OWNER/stride-dengue-analysis
cd stride-dengue-analysis
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Before opening a pull request

```bash
pytest                       # all tests must pass (synthetic fixtures only)
pytest --cov=stride_analysis # keep coverage high
ruff check src tests         # lint
mypy                         # type-check
```

## Conventions

- **No real data in git.** Tests use the synthetic generator
  (`stride_analysis._synthetic`) and the committed example dataset only.
- **No generated outputs in git.** Everything under `outputs/` is regenerable.
- **Fail loudly.** New validation raises a typed `StrideAnalysisError` subclass
  with an actionable message; never silently drop or coerce data.
- **Keep the two data levels separate.** Replicate observations (Level 1) and
  STRIDE summaries (Level 2) are distinct canonical tables and must not be
  merged.
- **Reusable subpackages.** New functionality goes into `io/`, `validation/`,
  `canonical/`, or `models/`; the stage orchestrators only compose them.
- **This is S0.** Do not add biology, statistics, or figures here — those are
  later stages.

## Reporting issues

Please include the framework version, a minimal reproduction (ideally against a
synthetic dataset), and the full error message.
