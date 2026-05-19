# Internship

Code for my 2026 Internship. Builds GWIPS-viz track hubs from p-shifted RiboSeq bigWigs.

## How it works
<img width="5341" height="6159" alt="Command Flow for Sample-2026-05-18-124244" src="https://github.com/user-attachments/assets/c80026fc-a2c7-4b5d-98fa-5d820d00aa05" />

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
