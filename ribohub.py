"""
Build a Gwips RiboSeq track hub from sorted bigWig files.

Layout expected on disk:
    {data_dir}/{srr[:6]}/{srr[6:8]}/{filename}

Filenames follow:
    {SRR}_pshifted_{forward,reverse}.{bw,bigWig}            (all reads)
    {SRR}_pshifted_unique_{forward,reverse}.{bw,bigWig}     (unique only)
    {SRR}_pshifted_multimapped_{forward,reverse}.{bw,bigWig} (multimapped only)

Output formats:
    directory   : standard multi-file hub (hub.txt + genomes.txt + hg38/trackDb.txt + symlinks)
    single-file : useOneFile hub (everything in one .hub.txt, no genomes.txt or trackDb.txt)
"""
import csv
import logging
import shutil
import sys
from dataclasses import dataclass  
from pathlib import Path
from typing import Literal, NamedTuple  


import click
import trackhub  # type: ignore


# ----- type aliases --------------------------------------------------------
Strand = Literal["forward", "reverse"]
Kind = Literal["unique", "multimapped"]  # None separately means "all reads / bare file"

class FileKey(NamedTuple):
    """Dict key identifying one bigWig variant for a sample."""
    strand: Strand
    kind: Kind | None  # None means "all reads" (no kind suffix in filename)
#Named tuples are used so the discovery dict is keyed by something readable (e.g. FileKey("forward", "unique")) 
#instead of a raw tuple ("forward", "unique") 
#and to allow type annotations on the components.

class BigwigEntry(NamedTuple):
    """One classified bigWig file: which strand, which kind, where it lives."""
    strand: Strand
    kind: Kind | None
    path: Path


@dataclass(frozen=True)
class BuildContext:
    """Run-wide configuration passed to every per-sample builder."""
    base_url: str
    colors: dict[str, str]
    label_fields: tuple[str, ...]
    include_mode: str
# frozen=True so the context can't be mutated mid-build, 
# which makes it so there can be no changes mid build.

# ----- defaults (informative: runtime values come from CLI) ----------------
#Look into color brewer or similar for better defaults if desired. These are just a starting point.
DEFAULT_COLORS: dict[str, str] = {
    "fwd":       "#FF0000",   # red
    "rev":       "#0000FF",   # blue
    "fwd_multi": "#FF9696",   # light red
    "rev_multi": "#6496FF",   # light blue
}
VIEW_LIMITS: str = "-127:127"

# Default columns to surface in a sample's long label in render order.
# Chosen from RiboSeqOrg metadata based on fill rate + biological signal. --> #TODO: add numbers
DEFAULT_LABEL_FIELDS: tuple[str, ...] = (
    "ScientificName",
    "LIBRARYTYPE",
    "TISSUE",
    "CELL_LINE",
    "CONDITION",
    "INHIBITOR",
    "TIMEPOINT",
    #"REPLICATE", might not be needed
)

# Values that RiboSeqOrg uses to mean "no data". Case-insensitive.
NULL_VALUES: frozenset[str] = frozenset({
    "", "nana", "na", "n/a", "none", "null", "nan", "-", "unknown",
})

log = logging.getLogger("ribohub")


# ----- input parsing -------------------------------------------------------
#This is so the main command function can assume it's working with a clean set of SRR IDs
#and not have to worry about the various ways users might specify them.
#Overall give you more options
def parse_samples(value: str) -> set[str]:
    """Parse --samples into a non-empty set of SRR IDs.

    Accepts: single ID, comma-separated list, .txt file (one per line),
    or .csv file (first column).
    """
    path = Path(value)
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as f:
            if path.suffix.lower() == ".csv":
                ids = {row[0].strip() for row in csv.reader(f)
                       if row and row[0].strip()}
            else:
                ids = {line.strip() for line in f if line.strip()}
    else:
        ids = {s.strip() for s in value.split(",") if s.strip()}

    if not ids:
        raise click.BadParameter("No sample IDs found in --samples value.")
    return ids

# UCSC/Gwips trackDb expects "R,G,B" decimal strings (e.g. "255,0,0"), not hex.
def hex_to_trackhub_rgb(value: str) -> str:
    """Validate hex color and convert to trackhub RGB format ('R,G,B')."""
    v = value.lstrip("#")
    if len(v) != 6 or not all(c in "0123456789abcdefABCDEF" for c in v):
        raise click.BadParameter(f"Invalid hex color: {value!r}. Expected #RRGGBB.")
    return trackhub.helpers.hex2rgb(value)
#Validate first so a bad --color-fwd value fails with a clear CLI
#error, not a cryptic trackhub crash mid-build.


def parse_label_fields(value: str) -> tuple[str, ...]:
    """Parse --label-fields into an ordered, de-duplicated tuple."""
    fields: list[str] = []
    seen: set[str] = set()
    for raw in value.split(","):
        col = raw.strip()
        if not col or col in seen:
            continue
        seen.add(col)
        fields.append(col)
    if not fields:
        raise click.BadParameter("--label-fields cannot be empty.")
    return tuple(fields)
#Manual de-dup loop (rather than `set(value.split(","))`) preserves
#the order the user typed, which determines render order in labels.

# Mirrors sorting.sh: organizes files by first 6 + next 2 chars of the SRR ID.
def srr_to_dir(srr_id: str, data_dir: Path) -> Path:
    """Compute on-disk directory for a sample. Mirrors sorting.sh."""
    return data_dir / srr_id[:6] / srr_id[6:8]
#For better understanding look at the sorting.sh script in the same repo, 
#which organizes files into this structure.

def classify_bigwig(path: Path) -> BigwigEntry | None:
    """Recognize a bigWig file's strand and kind. Returns None if unrecognized."""
    fname = path.name
    #Looks for the expected patterns in the filename to determine strand and kind.
    if "_forward" in fname:
        strand: Strand = "forward"
    elif "_reverse" in fname:
        strand = "reverse"
    else:
        return None #No strand info in filename → skip this file
    #Look for kind info, but it's optional.
    if "_unique" in fname:
        kind: Kind | None = "unique"
    elif "_multimapped" in fname:
        kind = "multimapped"
    else:
        kind = None

    return BigwigEntry(strand=strand, kind=kind, path=path)


# ----- metadata ------------------------------------------------------------
def _has_value(value: str | None) -> bool:
    """Return True if value is non-empty and not a known null sentinel.

    `value` is typed `str | None` because metadata dicts may legitimately
    not have a given column at all (returning None from .get), and many
    columns in the RiboSeqOrg CSV are blank strings to mean "missing".
    Both cases collapse to "no real data here".
    """
    if value is None:
        return False
    return value.strip().lower() not in NULL_VALUES


def _clean_value(col: str, val: str) -> str: #Doesnt really matter and could be removed if REPLICATE doesnt get used
    """Normalize known cosmetic quirks in RiboSeqOrg metadata.

    Currently:
      - REPLICATE arrives as '1.0', '2.0' (CSV float coercion) → '1', '2'.
    """
    v = val.strip()
    if not v:
        return v
    if col == "REPLICATE" and v.endswith(".0") and v[:-2].isdigit():
        return v[:-2]
    return v

#This is where the metadata CSV is loaded and cleaned up a bit. 
#The main thing is to key it by SRR ID so it's easy to look up when building each sample's tracks. 
#Also does some basic validation of the expected columns and logs warnings about duplicates or missing metadata.
def load_metadata(
    path: Path,
    srr_column: str,
    label_fields: tuple[str, ...] | None = None,
) -> dict[str, dict[str, str]]:
    """Load metadata CSV into an SRR-keyed dict of full rows."""
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames or []
        if srr_column not in header:
            raise click.BadParameter(
                f"Column {srr_column!r} not found in metadata CSV. "
                f"Available columns: {', '.join(header)}"
            )
        if label_fields:
            unknown = [c for c in label_fields if c not in header]
            if unknown:
                raise click.BadParameter(
                    f"--label-fields contains column(s) not in CSV: "
                    f"{', '.join(unknown)}. Available: {', '.join(header)}"
                )
        #Overall QC of the metadata CSV.
        rows: dict[str, dict[str, str]] = {}
        duplicates: list[str] = []
        for raw_row in reader:
            srr = (raw_row.get(srr_column) or "").strip()
            if not srr:
                continue
            row = {k: _clean_value(k, v or "") for k, v in raw_row.items()}
            if srr in rows:
                duplicates.append(srr)
            rows[srr] = row

    if duplicates:
        log.warning(
            f"Metadata CSV has {len(duplicates)} duplicate SRR(s); "
            f"keeping last occurrence. Examples: {', '.join(duplicates[:5])}"
        )
    log.info(f"Loaded metadata for {len(rows)} sample(s).")
    return rows


def long_label_for_sample(
    srr_id: str,
    metadata_row: dict | None,
    label_fields: tuple[str, ...],
) -> str:
    """Build the rich long label used at the supertrack level."""
    if not metadata_row: #No metadata just gives the SRR ID as label
        return srr_id
    parts = [srr_id]
    for col in label_fields:
        val = metadata_row.get(col) #Skips missing columns gracefully (e.g. if user specified --label-fields that aren't in the CSV)
        if _has_value(val):
            parts.append(val.strip())
    return " | ".join(parts)


# ----- discovery -----------------------------------------------------------
#This is kind of the mapbuilder of the whole operation:
#It takes the set of requested SRR IDs and looks on disk for files that match the expected patterns.
#It classifies them by strand and kind, and returns a nested dict of what it found,
#along with sets of which SRR IDs were missing entirely or had no valid bigWigs.
def discover_samples(
    sample_set: set[str],
    data_dir: Path,
) -> tuple[dict[str, dict[FileKey, str]], set[str], set[str]]:
    """For each requested SRR, find files on disk and classify them."""
    found: dict[str, dict[FileKey, str]] = {}
    dir_missing: set[str] = set()
    dir_empty: set[str] = set()

    for srr_id in sorted(sample_set):
        sample_dir = srr_to_dir(srr_id, data_dir)

        if not sample_dir.is_dir():
            dir_missing.add(srr_id)
            continue

        bws = sorted(
            list(sample_dir.glob(f"{srr_id}_*.bw")) +
            list(sample_dir.glob(f"{srr_id}_*.bigWig"))
        )
        if not bws:
            dir_empty.add(srr_id)
            continue

        classified: dict[FileKey, str] = {}
        for path in bws:
            entry = classify_bigwig(path)
            if entry is None:
                log.debug(f"Skipping unrecognized file: {path}")
                continue
            rel = path.relative_to(data_dir).as_posix()
            classified[FileKey(entry.strand, entry.kind)] = rel

        if classified:
            found[srr_id] = classified
        else:
            dir_empty.add(srr_id)

    return found, dir_missing, dir_empty


def report_discovery(
    found: dict[str, dict],
    dir_missing: set[str],
    dir_empty: set[str],
    sample_set: set[str],
) -> None:
    """Print a summary of what was found vs missing/empty."""
    click.echo(f"Requested: {len(sample_set)} samples")
    click.echo(f"  Found:           {len(found)}")
    if dir_missing:
        click.echo(f"  Dir missing:     {len(dir_missing)} → "
                   f"{', '.join(sorted(dir_missing))}")
    if dir_empty:
        click.echo(f"  Dir empty:       {len(dir_empty)} → "
                   f"{', '.join(sorted(dir_empty))}")
    click.echo("")


# ----- track building ------------------------------------------------------
#These functions take the classified files for one sample and build the appropriate trackhub.Track, trackhub.CompositeTrack 
#or trackhub.AggregateTrack objects according to the requested --include mode.

#kwargs is used to build the trackhub.Track objects in a more compact way, since they have a lot of parameters.
#and to add the negateValues parameter conditionally without repeating a lot of code.

def _make_plain_track(
    srr_id: str,
    strand: Strand,
    kind: Kind | None,
    rel_path: str,
    ctx: BuildContext,
    negate: bool = False,
    metadata_row: dict | None = None,
) -> trackhub.Track:
    """Build a single plain Track for a strand/kind combo."""
    color = ctx.colors["fwd"] if strand == "forward" else ctx.colors["rev"]
    name_suffix = f"_{kind}" if kind else ""
    rich = long_label_for_sample(srr_id, metadata_row, ctx.label_fields)
    strand_suffix = f" ({strand}{', ' + kind if kind else ''})"
    kwargs = dict(
        name=trackhub.helpers.sanitize(f"{srr_id}{name_suffix}_{strand}"),
        url=f"{ctx.base_url}/{rel_path}",
        tracktype="bigWig",
        short_label=srr_id,
        long_label=rich + strand_suffix,
        color=color,
        visibility="full",
    )
    if negate:
        kwargs["negateValues"] = "on"
    return trackhub.Track(**kwargs)

def build_aggregate(
    srr_id: str,
    files: dict[FileKey, str],
    ctx: BuildContext,
    metadata_row: dict | None = None,
) -> trackhub.AggregateTrack | None:
    """Strand-overlay aggregate using bare (all-reads) files."""
    fwd = files.get(FileKey("forward", None))
    rev = files.get(FileKey("reverse", None))
    if not (fwd and rev):
        return None

    #Forms an overall with a long label that includes metadata, 
    #and two subtracks for forward and reverse strands colored differently 
    #and with negateValues on the reverse strand to flip it in the browser.
    rich = long_label_for_sample(srr_id, metadata_row, ctx.label_fields)
    agg = trackhub.AggregateTrack(
        name=trackhub.helpers.sanitize(f"{srr_id}_agg"),
        tracktype="bigWig",
        short_label="strand overlay",
        long_label=f"{rich} (strand overlay)",
        aggregate="transparentOverlay",
        showSubtrackColorOnUi="on",
        visibility="full",
        viewLimits=VIEW_LIMITS,
        autoScale="off",
        maxHeightPixels="100:50:8",
    )
    agg.add_subtrack(trackhub.Track(
        name=trackhub.helpers.sanitize(f"{srr_id}_agg_forward"),
        url=f"{ctx.base_url}/{fwd}",
        tracktype="bigWig",
        short_label="forward",
        long_label=f"{srr_id} forward",
        color=ctx.colors["fwd"],
    ))
    agg.add_subtrack(trackhub.Track(
        name=trackhub.helpers.sanitize(f"{srr_id}_agg_reverse"),
        url=f"{ctx.base_url}/{rev}",
        tracktype="bigWig",
        short_label="reverse",
        long_label=f"{srr_id} reverse",
        color=ctx.colors["rev"],
        negateValues="on",
    ))
    return agg

#This needs both unique and multimapped bigWigs for a given strand to be present, 
#but it gives a nice direct comparison of them colored differently 
#and with a long label that includes metadata.
def build_kind_composite(
    srr_id: str,
    strand: Strand,
    files: dict[FileKey, str],
    ctx: BuildContext,
    negate: bool = False,
    metadata_row: dict | None = None,
) -> trackhub.CompositeTrack | None:
    """Unique-vs-multimapped composite for one strand."""
    uniq = files.get(FileKey(strand, "unique"))
    multi = files.get(FileKey(strand, "multimapped"))
    if not (uniq and multi):
        return None

    if strand == "forward":
        base_color, light_color = ctx.colors["fwd"], ctx.colors["fwd_multi"]
        suffix = "FW"
    else:
        base_color, light_color = ctx.colors["rev"], ctx.colors["rev_multi"]
        suffix = "REV"

    rich = long_label_for_sample(srr_id, metadata_row, ctx.label_fields)
    comp = trackhub.CompositeTrack(
        name=trackhub.helpers.sanitize(f"{srr_id}_{suffix}"),
        tracktype="bigWig",
        short_label=f"{strand} kinds",
        long_label=f"{rich} ({strand}: unique vs multimapped)",
        visibility="full",
        viewLimits=VIEW_LIMITS,
        autoScale="off",
    )
    uniq_kwargs = dict(
        name=trackhub.helpers.sanitize(f"{srr_id}_{suffix}_unique"),
        url=f"{ctx.base_url}/{uniq}",
        tracktype="bigWig",
        short_label="unique",
        long_label=f"{srr_id} {strand} unique",
        color=base_color,
    )
    multi_kwargs = dict(
        name=trackhub.helpers.sanitize(f"{srr_id}_{suffix}_multimapped"),
        url=f"{ctx.base_url}/{multi}",
        tracktype="bigWig",
        short_label="multi",
        long_label=f"{srr_id} {strand} multimapped",
        color=light_color,
    )
    if negate:
        uniq_kwargs["negateValues"] = "on"
        multi_kwargs["negateValues"] = "on"
    comp.add_subtrack(trackhub.Track(**uniq_kwargs))
    comp.add_subtrack(trackhub.Track(**multi_kwargs))
    return comp

#This is the highest parent track for one sample, which will go under the main trackDb.
#Whether it contains just the minimal unique bigWigs or the aggregates or the composites
#depends on the --include mode, but it always has the same long label format that can include metadata.
def build_sample_supertrack(
    srr_id: str,
    files: dict[FileKey, str],
    metadata_row: dict | None,
    ctx: BuildContext,
) -> trackhub.SuperTrack | None:
    """Per-sample SuperTrack containing the tracks for ctx.include_mode."""
    sample_super = trackhub.SuperTrack(
        name=trackhub.helpers.sanitize(srr_id),
        short_label=srr_id,
        long_label=long_label_for_sample(srr_id, metadata_row, ctx.label_fields),
    )

    added_anything = False

    if ctx.include_mode == "minimal":
        for strand, negate in (("forward", False), ("reverse", True)):
            uniq = files.get(FileKey(strand, "unique"))
            if uniq:
                sample_super.add_tracks(
                    _make_plain_track(
                        srr_id, strand, "unique", uniq,
                        ctx=ctx, negate=negate,
                        metadata_row=metadata_row,
                    )
                )
                added_anything = True

    if ctx.include_mode in ("aggregate", "full"):
        agg = build_aggregate(srr_id, files, ctx, metadata_row=metadata_row)
        if agg:
            sample_super.add_tracks(agg)
            added_anything = True

    if ctx.include_mode in ("composite", "full"):
        fw = build_kind_composite(
            srr_id, "forward", files, ctx,
            negate=False, metadata_row=metadata_row,
        )
        if fw:
            sample_super.add_tracks(fw)
            added_anything = True
        rev = build_kind_composite(
            srr_id, "reverse", files, ctx,
            negate=True, metadata_row=metadata_row,
        )
        if rev:
            sample_super.add_tracks(rev)
            added_anything = True

    return sample_super if added_anything else None


# ----- output writers ------------------------------------------------------
#This might get changed later on, because it literally deletes the old output directory if it exists.
def write_directory_hub(hub, hub_name: str, output_dir: Path) -> Path:
    """Write the standard multi-file hub. Returns the path to the hub.txt file."""
    staging = output_dir / hub_name
    if staging.exists():
        shutil.rmtree(staging)
    trackhub.upload.stage_hub(hub, staging=str(staging))
    return staging / f"{hub_name}.hub.txt"

#This is the single file option which is favorable for Galaxy integration.
#It puts everything in one .hub.txt file with "useOneFile on", 
#which Gwips can read directly without needing the usual directory structure or separate trackDb.txt.
def write_single_file_hub(
    hub_name: str,
    short_label: str,
    long_label: str,
    email: str,
    genome: str,
    trackdb,
    output_dir: Path,
) -> Path:
    """Write a useOneFile hub."""
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{hub_name}.hub.txt"

    lines: list[str] = [
        f"hub {hub_name}",
        f"shortLabel {short_label}",
        f"longLabel {long_label}",
        "useOneFile on",
        f"email {email}",
        f"genome {genome}",
        "",
    ]
    for track in trackdb.children:
        lines.append(str(track))
        lines.append("")

    out_path.write_text("\n".join(lines))
    return out_path


# ----- CLI -----------------------------------------------------------------

@click.group()
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
@click.pass_context
def cli(ctx: click.Context, verbose: bool) -> None:
    """RiboHub: Build a Gwips RiboSeq track hub from organized bigWig files."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose


@cli.command()
# Required inputs
@click.option("--samples", required=True,
              help="Sample selection. Single ID, comma list, .txt, or .csv (first column).")
@click.option("--data-dir", required=True, envvar="RIBOHUB_DATA_DIR",
              type=click.Path(exists=True, file_okay=False, path_type=Path),
              help="Root directory containing sorted bigWig files.")
@click.option("--output-dir", required=True, envvar="RIBOHUB_OUTPUT_DIR",
              type=click.Path(file_okay=False, path_type=Path),
              help="Directory where the hub will be written.")
@click.option("--base-url", required=True, envvar="RIBOHUB_BASE_URL",
              help="Public URL where bigWig files are served (must be reachable by Gwips).")
# Metadata
@click.option("--metadata", "metadata_path", default=None,
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help="Optional RiboSeqOrg-style metadata CSV. Rows matched to "
                   "samples by --srr-column.")
@click.option("--srr-column", default="Run", show_default=True,
              help="Column name in metadata CSV that holds the SRR ID.")
@click.option("--label-fields",
              default=",".join(DEFAULT_LABEL_FIELDS), show_default=True,
              help="Comma-separated metadata columns to include in each sample's "
                   "long label.")
# Hub configuration
@click.option("--genome", default="hg38", show_default=True,
              help="Gwips genome assembly.")
@click.option("--hub-name", default="RiboSeqHub", show_default=True,
              help="Hub directory name (no spaces).")
@click.option("--email", default="your@email.com", show_default=True,
              help="Contact email written into hub.txt.")
# Output
@click.option("--output-format", default="directory", show_default=True,
              type=click.Choice(["directory", "single-file"]),
              help="directory: multi-file hub; single-file: useOneFile hub.")
# Build behavior
@click.option("--include", "include_mode", default="minimal", show_default=True,
              type=click.Choice(["minimal", "aggregate", "composite", "full"]),
              help="What to build per sample.")
@click.option("--strict", is_flag=True,
              help="Exit non-zero if any requested sample is missing or partial.")
@click.option("--dry-run", is_flag=True,
              help="Report what would be built without writing the hub.")
# Colors
@click.option("--color-fwd", default=DEFAULT_COLORS["fwd"], show_default=True,
              help="Forward strand color (hex).")
@click.option("--color-rev", default=DEFAULT_COLORS["rev"], show_default=True,
              help="Reverse strand color (hex).")
@click.option("--color-fwd-multi", default=DEFAULT_COLORS["fwd_multi"], show_default=True,
              help="Forward multimapped color (hex).")
@click.option("--color-rev-multi", default=DEFAULT_COLORS["rev_multi"], show_default=True,
              help="Reverse multimapped color (hex).")
def generate(
    samples: str,
    data_dir: Path,
    output_dir: Path,
    base_url: str,
    metadata_path: Path | None,
    srr_column: str,
    label_fields: str,
    genome: str,
    hub_name: str,
    email: str,
    output_format: str,
    include_mode: str,
    strict: bool,
    dry_run: bool,
    color_fwd: str,
    color_rev: str,
    color_fwd_multi: str,
    color_rev_multi: str,
) -> None:
    """Generate a Gwips track hub for the given samples."""
    base_url = base_url.rstrip("/")

    #Validate and convert colors to trackhub format.
    colors: dict[str, str] = {
        "fwd":       hex_to_trackhub_rgb(color_fwd),
        "rev":       hex_to_trackhub_rgb(color_rev),
        "fwd_multi": hex_to_trackhub_rgb(color_fwd_multi),
        "rev_multi": hex_to_trackhub_rgb(color_rev_multi),
    }

    parsed_label_fields = parse_label_fields(label_fields)

    #Freeze the run config into one immutable context object passed to
    #every builder, instead of threading 4-5 args through every call.
    ctx = BuildContext(
        base_url=base_url,
        colors=colors,
        label_fields=parsed_label_fields,
        include_mode=include_mode,
    )

    sample_set = parse_samples(samples)
    log.info(f"Resolved {len(sample_set)} requested samples.")

    found, dir_missing, dir_empty = discover_samples(sample_set, data_dir)
    report_discovery(found, dir_missing, dir_empty, sample_set)

    metadata: dict[str, dict[str, str]] = {}
    if metadata_path:
        metadata = load_metadata(metadata_path, srr_column, parsed_label_fields)
        missing_meta = sorted(set(found) - set(metadata))
        if missing_meta:
            log.warning(
                f"{len(missing_meta)} discovered sample(s) have no metadata row. "
                f"Examples: {', '.join(missing_meta[:5])}"
            )

    if not found:
        click.echo("No samples found. Nothing to build.", err=True)
        sys.exit(1)


    #Dry-run intentionally runs AFTER discovery and metadata loading
    #so the user sees real numbers and column errors before committing to a build.
    if dry_run:
        click.echo(f"Dry run: would build {len(found)} sample(s) "
                   f"with --include {include_mode} and "
                   f"--output-format {output_format}.")
        if strict and (dir_missing or dir_empty):
            click.echo(f"  Note: --strict would reject this run "
                       f"({len(dir_missing)} missing, {len(dir_empty)} empty).",
                       err=True)
        if metadata:
            with_meta = sum(1 for s in found if s in metadata)
            click.echo(f"Metadata matched: {with_meta}/{len(found)} sample(s).")
        return

    if strict and (dir_missing or dir_empty):
        click.echo("--strict: refusing to build a partial hub.", err=True)
        sys.exit(1)

    hub, _, _, trackdb = trackhub.default_hub(
        hub_name=hub_name,
        short_label=hub_name,
        long_label=f"{hub_name} Track Hub",
        genome=genome,
        email=email,
    )

    samples_built = 0
    for srr_id in sorted(found.keys()):
        sample_super = build_sample_supertrack(
            srr_id,
            found[srr_id],
            metadata.get(srr_id),
            ctx,
        )
        if sample_super:
            trackdb.add_tracks(sample_super)
            samples_built += 1
        else:
            log.warning(f"{srr_id}: no tracks built (insufficient files for "
                        f"--include {include_mode})")

    if samples_built == 0:
        click.echo("No samples had enough files for the chosen --include mode.",
                   err=True)
        sys.exit(1)

    if output_format == "directory":
        hub_file = write_directory_hub(hub, hub_name, output_dir)
        public_url = f"{base_url}/{hub_name}/{hub_name}.hub.txt"
    else:
        hub_file = write_single_file_hub(
            hub_name=hub_name,
            short_label=hub_name,
            long_label=f"{hub_name} Track Hub",
            email=email,
            genome=genome,
            trackdb=trackdb,
            output_dir=output_dir,
        )
        public_url = f"{base_url}/{hub_name}.hub.txt"

    click.echo(f"Built hub for {samples_built} sample(s) "
               f"with --include {include_mode} and --output-format {output_format}.")
    click.echo(f"Wrote: {hub_file}")
    click.echo(f"Hub URL: {public_url}")


if __name__ == "__main__":
    cli()
