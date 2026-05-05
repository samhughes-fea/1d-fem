# Benchmark ladder standard

This document standardizes how validation cases are promoted across the repository.

## Benchmark levels

### Level 0 — Smoke

- artifact existence
- basic positivity/finiteness
- dispatch and runner wiring

### Level 1 — Pinned acceptance

- repository-pinned outputs
- stable artifact contract
- explicit benchmark job and README

### Level 2 — Calibrated repository reference

- pinned scalar or field target
- explicit tolerance policy
- documented assumptions and quantity compared

### Level 3 — External-reference validated

- benchmark tied to handbook, textbook, or independent FE reference
- tolerance backed by external source or documented engineering rationale

## Required fields for every benchmark doc

Each benchmark document must define:

- purpose
- geometry / mesh
- BC convention
- load / forcing convention
- quantity compared
- source of truth
- tolerance policy
- scope limitations

## Promotion rule

A benchmark should only move up one level at a time:

- Smoke → Pinned acceptance
- Pinned acceptance → Calibrated repository reference
- Calibrated repository reference → External-reference validated

## Test expectations by level

| Level | Typical test expectation |
|---|---|
| Smoke | files/logs exist, run completes |
| Pinned acceptance | artifacts + stable benchmark-specific metadata |
| Calibrated repository reference | scalar/path/field target within tolerance |
| External-reference validated | external-source comparison within documented tolerance |
