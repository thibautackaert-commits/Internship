# Internship

Code for my 2026 Internship. Builds GWIPS-viz track hubs from p-shifted RiboSeq bigWigs.

## How it works
<img width="4891" height="6100" alt="Untitled diagram-2026-05-18-112200" src="https://github.com/user-attachments/assets/c30604a3-0050-481a-8009-70da60536ec1" />

## Files

- `hub_builder.py` - the CLI tool
- `sorting.sh` - organizes raw bigWigs into the expected layout

## Use

```bash
pip install trackhub click
python hub_builder.py generate --help
```

## Status

Work in progress.
