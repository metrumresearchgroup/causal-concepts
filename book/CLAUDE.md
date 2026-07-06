# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A Quarto book titled *Causal Concepts for Pharmacometricians*, jointly developed by MetrumRG and Genentech/Roche. 
It covers causal DAGs, estimands, and related causal inference topics aimed at pharmacometrics practitioners.

## Build Commands

All commands run from the `book/` directory (the working directory).

```bash
# Render the full book
quarto render

# Live preview with hot reload
quarto preview

# Restore R package environment (run once after cloning or after renv.lock changes)
Rscript -e "renv::restore()"
```

There is no lint or test step — the primary validation is a successful `quarto render`.

## Architecture

### Book Structure

Defined in `_quarto.yml`. Chapters are individual `.qmd` files:
- `index.qmd` — Preface
- `Chapter0/intro.qmd` — Introduction (Why causal inference, motivating examples)
- `Chapter1/dags.qmd` — Causal DAGs (conventions, examples)
- `Chapter2/estimands.qmd` — Estimands, potential outcomes, ATE/CATE, dose-exposure-response
- `references.qmd` — Bibliography (sourced from `references.bib`)

### TikZ DAG Rendering

DAGs are rendered inside R code chunks using the knitr `tikz` engine, not as standalone `.tex` files. 
Each chapter that renders DAGs sets the knitr engine options in a setup chunk:

```r
knitr::opts_chunk$set(engine.opts = list(
  template = "../tikz_template.tex",
  convert = "dvisvgm"
))
```

`tikz_template.tex` provides the LaTeX preamble (TikZ libraries, `amsmath`, etc.) and 
the `%% TIKZ_CODE %%` placeholder that knitr fills in. 
The `dvisvgm` converter is required to preserve math italics in SVG output.

### DAG Visual Conventions

- **Blue edges**: the causal pathway of interest
- **Red edges**: biasing/backdoor paths
- **Gray edges**: other causal dependencies
- **Dashed edges**: speculative/uncertain causal dependencies
- **Rectangular nodes** (`\node[draw]`): adjusted/conditioned-on variables

### Output and Ignored Paths

The `.gitignore` excludes `/.quarto/` (build cache) and `/book/` (rendered HTML output). The `dags/` directory contains a legacy standalone LaTeX file — it is not part of the Quarto build pipeline.

### R Package Management

Uses `renv` (lockfile at `renv/renv.lock` in the repo root). Key packages: `knitr`, `quarto`, `dplyr`, `ggplot2`, `dagitty`. R version: 4.5.1.
