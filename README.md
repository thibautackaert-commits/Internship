# Internship

Code for my 2026 Internship. Builds GWIPS-viz track hubs from p-shifted RiboSeq bigWigs.

## How it works
<img width="4891" height="6100" alt="Untitled diagram-2026-05-18-112114" src="https://github.com/user-attachments/assets/21d6c421-f903-4d80-8fde-9680ff7859f1" />
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
