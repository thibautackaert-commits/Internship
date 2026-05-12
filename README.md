# RiboHub

A command-line tool for building [UCSC Genome Browser](https://genome.ucsc.edu/) track hubs from P-shifted bigWig files produced by a RiboSeq pipeline.

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Expected File Layout](#expected-file-layout)
- [Usage](#usage)
  - [generate](#generate)
  - [Options](#options)
  - [Sample Selection](#sample-selection)
  - [Include Modes](#include-modes)
  - [Colors](#colors)
- [Examples](#examples)
- [Output](#output)
- [Environment Variables](#environment-variables)
- [Troubleshooting](#troubleshooting)

---

## Overview

RiboHub takes a set of SRR sample IDs, locates their sorted bigWig files on disk, and builds a fully staged UCSC track hub. Each sample gets its own SuperTrack container. Depending on the `--include` mode you choose, RiboHub will generate:

- Plain unique-read tracks per strand (minimal)
- A transparent strand overlay aggregate (aggregate)
- Unique vs. multimapped composite tracks per strand (composite)
- All of the above (full)

---

## Requirements

- Python ≥ 3.10
- [`trackhub`](https://github.com/daler/trackhub)
- [`click`](https://click.palletsprojects.com/)

Install dependencies:

```bash
pip install trackhub click
```

---

## Installation

```bash
git clone https://github.com/your-org/ribohub.git
cd ribohub
pip install -e .
```

Or run directly:

```bash
python ribohub.py --help
```

---

## Expected File Layout

RiboHub expects bigWig files to be pre-sorted into the following directory structure (matching the output of `sorting.sh`):

```
{data_dir}/
└── {SRR[:6]}/
    └── {SRR[6:8]}/
        ├── {SRR}_pshifted_forward.bigWig          # all reads, forward
        ├── {SRR}_pshifted_reverse.bigWig          # all reads, reverse
        ├── {SRR}_pshifted_unique_forward.bigWig
        ├── {SRR}_pshifted_unique_reverse.bigWig
        ├── {SRR}_pshifted_multimapped_forward.bigWig
        └── {SRR}_pshifted_multimapped_reverse.bigWig
```

**Example** for `SRR1234567`:

```
data/
└── SRR123/
    └── 45/
        ├── SRR1234567_pshifted_forward.bigWig
        ├── SRR1234567_pshifted_reverse.bigWig
        ├── SRR1234567_pshifted_unique_forward.bigWig
        └── ...
```

Files not matching this naming scheme are skipped with a debug warning.

---

## Usage

### generate

```bash
python ribohub.py generate [OPTIONS]
```

Builds and stages a UCSC track hub for the specified samples.

### Options

| Option | Required | Default | Description |
|---|---|---|---|
| `--samples` | ✅ | — | Sample selection (see below) |
| `--data-dir` | ✅ | `$RIBOHUB_DATA_DIR` | Root directory of sorted bigWig files |
| `--output-dir` | ✅ | `$RIBOHUB_OUTPUT_DIR` | Directory where the hub will be written |
| `--base-url` | ✅ | `$RIBOHUB_BASE_URL` | Public URL serving `output-dir` |
| `--genome` | — | `hg38` | UCSC assembly (e.g. `mm10`, `dm6`) |
| `--hub-name` | — | `RiboSeqHub` | Hub directory name (no spaces) |
| `--email` | — | `your@email.com` | Contact email written into `hub.txt` |
| `--include` | — | `minimal` | Track layout mode (see below) |
| `--strict` | — | off | Exit with error if any sample is missing or partial |
| `--dry-run` | — | off | Report what would be built without writing files |
| `--color-fwd` | — | `#FF0000` | Forward strand color (hex) |
| `--color-rev` | — | `#0000FF` | Reverse strand color (hex) |
| `--color-fwd-multi` | — | `#FF9696` | Forward multimapped color (hex) |
| `--color-rev-multi` | — | `#6496FF` | Reverse multimapped color (hex) |
| `--verbose` | — | off | Enable debug logging |

### Sample Selection

`--samples` accepts four formats:

| Format | Example |
|---|---|
| Single SRR ID | `--samples SRR1234567` |
| Comma-separated list | `--samples SRR1234567,SRR7654321` |
| Plain text file (one ID per line) | `--samples samples.txt` |
| CSV file (first column used) | `--samples metadata.csv` |

### Include Modes

| Mode | Tracks built per sample | Files required |
|---|---|---|
| `minimal` | 1 unique-forward track + 1 unique-reverse track | `*_unique_forward.bigWig`, `*_unique_reverse.bigWig` |
| `aggregate` | 1 transparent strand overlay | `*_forward.bigWig`, `*_reverse.bigWig` (bare) |
| `composite` | FW composite (unique vs multi) + REV composite (unique vs multi) | all 4 strand+kind files |
| `full` | aggregate + both composites | all 6 files |

Samples that don't have enough files for the chosen mode are skipped with a warning. Use `--strict` to treat this as an error.

### Colors

All color options accept standard hex format (`#RRGGBB`). Colors are automatically converted to the `R,G,B` format required by UCSC.

---

## Examples

**Minimal hub for two samples:**

```bash
python ribohub.py generate \
  --samples SRR1234567,SRR7654321 \
  --data-dir /data/riboseq \
  --output-dir /var/www/hubs \
  --base-url https://myserver.example.com/hubs \
  --genome hg38 \
  --include minimal
```

**Full hub from a CSV, with custom colors, dry run first:**

```bash
python ribohub.py generate \
  --samples metadata.csv \
  --data-dir /data/riboseq \
  --output-dir /var/www/hubs \
  --base-url https://myserver.example.com/hubs \
  --genome mm10 \
  --include full \
  --color-fwd "#E63946" \
  --color-rev "#457B9D" \
  --dry-run
```

**Strict mode (fail if any sample is missing):**

```bash
python ribohub.py generate \
  --samples samples.txt \
  --data-dir /data/riboseq \
  --output-dir /var/www/hubs \
  --base-url https://myserver.example.com/hubs \
  --strict
```

---

## Output

After a successful run, the hub is staged under `{output-dir}/{hub-name}/`:

```
/var/www/hubs/RiboSeqHub/
├── hub.txt
├── genomes.txt
└── hg38/
    ├── trackDb.txt
    └── ... (symlinks or copies of bigWig files)
```

Load into UCSC by pasting the hub URL into **My Hubs**:

```
https://myserver.example.com/hubs/RiboSeqHub/hub.txt
```

---

## Environment Variables

You can set these instead of passing the flags every time:

```bash
export RIBOHUB_DATA_DIR=/data/riboseq
export RIBOHUB_OUTPUT_DIR=/var/www/hubs
export RIBOHUB_BASE_URL=https://myserver.example.com/hubs
```

---

## Troubleshooting

**"No samples found. Nothing to build."**
The directory structure doesn't match what RiboHub expects. Run with `--verbose` to see which paths were checked.

**"No tracks built (insufficient files for --include mode)"**
The sample directory exists but is missing files required by the chosen mode. Switch to `--include minimal` or check that your pipeline produced all expected bigWig variants.

**Hub loads in UCSC but shows no tracks**
Verify that `--base-url` is publicly reachable by UCSC's servers (not `localhost` or an internal IP). The bigWig files must be accessible at `{base-url}/{relative-path}`.

**Hex color error**
Colors must be exactly 6 hex digits with an optional `#` prefix, e.g. `#FF0000` or `FF0000`.
