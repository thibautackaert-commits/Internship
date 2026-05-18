# Internship

Code for my 2026 Internship. Builds GWIPS-viz track hubs from p-shifted RiboSeq bigWigs.

## How it works
```mermaid
flowchart TD
    Start([generate command]) --> Parse[parse_samplesvalidate colorsbuild BuildContext]
    Parse --> Discover[discover_sampleswalk srr_to_dirclassify bigWigs]
    Discover --> LoadMeta{--metadataprovided?}
    LoadMeta -->|yes| Meta[load_metadataSRR-keyed rows]
    LoadMeta -->|no| FoundCheck
    Meta --> FoundCheck{any samplesfound?}

    FoundCheck -->|no| Exit1[exit 1:nothing to build]
    FoundCheck -->|yes| DryRun{--dry-run?}

    DryRun -->|yes| PrintPlan[print planreturn]
    DryRun -->|no| Strict{--strictand gaps?}

    Strict -->|yes| Exit2[exit 1:refuse partial hub]
    Strict -->|no| DefaultHub[trackhub.default_hubcreate hub, genome, trackDb]

    DefaultHub --> Loop[for each SRR:build_sample_supertrack]

    Loop --> Include{--include mode}
    Include -->|minimal| Minimal[plain unique tracksforward + reverse]
    Include -->|aggregate| Agg[build_aggregatestrand overlay]
    Include -->|composite| Comp[build_kind_compositeunique vs multimapped]
    Include -->|full| Both[aggregate + composite]

    Minimal --> AnyBuilt{samples_built> 0?}
    Agg --> AnyBuilt
    Comp --> AnyBuilt
    Both --> AnyBuilt

    AnyBuilt -->|no| Exit3[exit 1:no valid samples]
    AnyBuilt -->|yes| Format{--output-format}

    Format -->|directory| DirHub[write_directory_hubstage_hub + symlinks]
    Format -->|single-file| SingleHub[write_single_file_hubuseOneFile .hub.txt]

    DirHub --> Done([echo hub file +public URL])
    SingleHub --> Done
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
