# Internship

Code for my 2026 Internship. Builds GWIPS-viz track hubs from p-shifted RiboSeq bigWigs.

## How it works
```mermaid
flowchart TD
    %%{init: {"layout": "elk"}}%%
    classDef input stroke:#38bdf8,fill:#f0f9ff,color:#0c4a6e
    classDef process stroke:#4ade80,fill:#f0fdf4,color:#14532d
    classDef decision stroke:#facc15,fill:#fefce8,color:#713f12
    classDef output stroke:#818cf8,fill:#eef2ff,color:#312e81
    classDef optional stroke:#fb923c,fill:#fff7ed,color:#7c2d12,stroke-dasharray:4 3

    subgraph "🚀 Command Flow"
        Start([Generate]):::input --> Parse[Parse samples<br/>and build context]:::process
        Parse --> Discover[Discover bigWigs<br/>on disk]:::process
        Discover --> Init[Initialize default hub]:::process
        Init --> Loop[Loop over SRRs<br/>build supertracks]:::process
        Loop --> Mode{--include mode}:::decision

        Mode -->|minimal| Min[Unique reads only<br/>forward + reverse]:::process
        Mode -->|aggregate| Agg[Strand overlay<br/>fwd + rev on one track]:::process
        Mode -->|composite| Comp[Unique vs multimapped<br/>grouped per strand]:::process
        Mode -->|full| Full[Overlay + composites<br/>everything]:::process

        Min --> Write[Write hub]:::output
        Agg --> Write
        Comp --> Write
        Full --> Write
        Write --> End([Print hub URL]):::output
    end

    subgraph "⚙️ Optional Flags"
        Meta[Load metadata CSV<br/>enrich sample labels]:::optional
        Dry[--dry-run<br/>print plan & exit]:::optional
        Strict[--strict<br/>abort on missing samples]:::optional
        Single[--output-format single-file<br/>write one .hub.txt]:::optional
    end

    Parse -.-> Meta
    Discover -.-> Dry
    Discover -.-> Strict
    Write -.-> Single
```

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
