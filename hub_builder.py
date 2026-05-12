"""
Build a UCSC RiboSeq track hub from sorted bigWig files.

Layout expected on disk:
    {data_dir}/{srr[:6]}/{srr[6:8]}/{srr}_pshifted_[unique_|multimapped_]{forward|reverse}.bigWig
    {data_dir}/{srr[:6]}/{srr[6:8]}/{srr}_pshifted_{forward|reverse}.bigWig    (bare = all reads)
"""
import csv
import logging
import shutil
import sys
from pathlib import Path

import click
import trackhub # type: ignore


# ----- defaults (informative; runtime values come from CLI) ----------------

DEFAULT_COLORS: dict[str, str] = {
    "fwd":       "#FF0000",   # red
    "rev":       "#0000FF",   # blue
    "fwd_multi": "#FF9696",   # light red
    "rev_multi": "#6496FF",   # light blue
}

VIEW_LIMITS: str = "-127:127"

log = logging.getLogger("ribohub")


# ----- small helpers -------------------------------------------------------

def parse_samples(value: str) -> set[str]:
    """Parse --samples into a non-empty set of SRR IDs.

    Accepts: single ID, comma-separated list, .txt file (one per line),
    or .csv file (first column).
    """
    path = Path(value)
    if path.is_file():
        with path.open(newline="") as f:
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


def parse_hex_color(value: str) -> str:
    """Validate hex color and convert to UCSC RGB format ('R,G,B')."""
    v = value.lstrip("#")
    if len(v) != 6 or not all(c in "0123456789abcdefABCDEF" for c in v):
        raise click.BadParameter(f"Invalid hex color: {value!r}. Expected #RRGGBB.")
    return trackhub.helpers.hex2rgb(value)


def srr_to_dir(srr_id: str, data_dir: Path) -> Path:
    """Compute on-disk directory for a sample. Mirrors sorting.sh."""
    return data_dir / srr_id[:6] / srr_id[6:8]


def classify_bigwig(path: Path) -> dict | None:
    """Recognize a bigWig file's strand and kind. Returns None if unrecognized.

    Returns {strand, kind, path}. kind=None means a 'bare' file (no kind suffix
    in filename — used for the aggregate view).
    """
    fname = path.name

    if "_forward" in fname:
        strand = "forward"
    elif "_reverse" in fname:
        strand = "reverse"
    else:
        return None

    if "_unique" in fname:
        kind: str | None = "unique"
    elif "_multimapped" in fname:
        kind = "multimapped"
    else:
        kind = None

    return {"strand": strand, "kind": kind, "path": path}


def long_label_for_sample(srr_id: str, metadata_row: dict | None = None) -> str:
    """Build the rich long label used at the supertrack level.
    Stub for V3 — returns just the SRR ID. Wire metadata in here later.
    """
    if not metadata_row:
        return srr_id
    parts = [srr_id]
    for col in ("ScientificName", "TISSUE", "CELL_LINE", "CONDITION"):
        val = metadata_row.get(col, "").strip()
        if val:
            parts.append(val)
    return " | ".join(parts)


# ----- discovery -----------------------------------------------------------

def discover_samples(
    sample_set: set[str],
    data_dir: Path,
) -> tuple[dict[str, dict[tuple[str, str | None], str]], set[str], set[str]]:
    """For each requested SRR, find files on disk and classify them.

    Returns:
        found:        {srr_id: {(strand, kind): rel_posix_path}}
        dir_missing:  set of SRR IDs with no directory on disk
        dir_empty:    set of SRR IDs with empty / unrecognized contents
    """
    found: dict[str, dict[tuple[str, str | None], str]] = {}
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

        classified: dict[tuple[str, str | None], str] = {}
        for path in bws:
            entry = classify_bigwig(path)
            if entry is None:
                log.debug(f"Skipping unrecognized file: {path}")
                continue
            rel = path.relative_to(data_dir).as_posix()
            classified[(entry["strand"], entry["kind"])] = rel

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

def _make_plain_track(
    srr_id: str,
    strand: str,
    kind: str | None,
    rel_path: str,
    base_url: str,
    colors: dict[str, str],
    negate: bool = False,
) -> trackhub.Track:
    """Build a single plain Track for a strand/kind combo."""
    color = colors["fwd"] if strand == "forward" else colors["rev"]
    name_suffix = f"_{kind}" if kind else ""
    return trackhub.Track(
        name=trackhub.helpers.sanitize(f"{srr_id}{name_suffix}_{strand}"),
        url=f"{base_url}/{rel_path}",
        tracktype="bigWig",
        short_label=srr_id,
        long_label=srr_id,
        color=color,
        visibility="full",
        negateValues="on" if negate else "off",
    )


def build_aggregate(
    srr_id: str,
    files: dict[tuple[str, str | None], str],
    base_url: str,
    colors: dict[str, str],
) -> trackhub.AggregateTrack | None:
    """Build the strand-overlay aggregate using bare (all-reads) files.
    Returns None if either bare file is missing.
    """
    fwd = files.get(("forward", None))
    rev = files.get(("reverse", None))
    if not (fwd and rev):
        return None

    agg = trackhub.AggregateTrack(
        name=trackhub.helpers.sanitize(f"{srr_id}_agg"),
        tracktype="bigWig",
        short_label="strand overlay",
        long_label=f"{srr_id} strand overlay (all reads)",
        aggregate="transparentOverlay",
        showSubtrackColorOnUi="on",
        visibility="full",
        viewLimits=VIEW_LIMITS,
        autoScale="off",
        maxHeightPixels="100:50:8",
    )
    agg.add_subtrack(trackhub.Track(
        name=trackhub.helpers.sanitize(f"{srr_id}_agg_forward"),
        url=f"{base_url}/{fwd}",
        tracktype="bigWig",
        short_label="forward",
        long_label=f"{srr_id} forward",
        color=colors["fwd"],
    ))
    agg.add_subtrack(trackhub.Track(
        name=trackhub.helpers.sanitize(f"{srr_id}_agg_reverse"),
        url=f"{base_url}/{rev}",
        tracktype="bigWig",
        short_label="reverse",
        long_label=f"{srr_id} reverse",
        color=colors["rev"],
        negateValues="on",
    ))
    return agg


def build_kind_composite(
    srr_id: str,
    strand: str,
    files: dict[tuple[str, str | None], str],
    base_url: str,
    colors: dict[str, str],
    negate: bool = False,
) -> trackhub.CompositeTrack | None:
    """Build the unique-vs-multimapped composite for one strand.
    Returns None unless both unique and multimapped exist for that strand.
    """
    uniq = files.get((strand, "unique"))
    multi = files.get((strand, "multimapped"))
    if not (uniq and multi):
        return None

    if strand == "forward":
        base_color, light_color = colors["fwd"], colors["fwd_multi"]
        suffix = "FW"
    else:
        base_color, light_color = colors["rev"], colors["rev_multi"]
        suffix = "REV"

    comp = trackhub.CompositeTrack(
        name=trackhub.helpers.sanitize(f"{srr_id}_{suffix}"),
        tracktype="bigWig",
        short_label=f"{strand} kinds",
        long_label=f"{srr_id} {strand} unique vs multimapped",
        visibility="full",
        viewLimits=VIEW_LIMITS,
        autoScale="off",
    )
    comp.add_subtrack(trackhub.Track(
        name=trackhub.helpers.sanitize(f"{srr_id}_{suffix}_unique"),
        url=f"{base_url}/{uniq}",
        tracktype="bigWig",
        short_label="unique",
        long_label=f"{srr_id} {strand} unique",
        color=base_color,
        negateValues="on" if negate else "off",
    ))
    comp.add_subtrack(trackhub.Track(
        name=trackhub.helpers.sanitize(f"{srr_id}_{suffix}_multimapped"),
        url=f"{base_url}/{multi}",
        tracktype="bigWig",
        short_label="multi",
        long_label=f"{srr_id} {strand} multimapped",
        color=light_color,
        negateValues="on" if negate else "off",
    ))
    return comp


def build_sample_supertrack(
    srr_id: str,
    files: dict[tuple[str, str | None], str],
    base_url: str,
    include_mode: str,
    colors: dict[str, str],
    metadata_row: dict | None = None,
) -> trackhub.SuperTrack | None:
    """Build the per-sample container with the right tracks for include_mode.

    include_mode:
      - 'minimal'   : just unique forward + unique reverse, plain tracks
      - 'aggregate' : aggregate only
      - 'composite' : FW composite + REV composite
      - 'full'      : aggregate + FW composite + REV composite
    """
    sample_super = trackhub.SuperTrack(
        name=trackhub.helpers.sanitize(srr_id),
        short_label=srr_id,
        long_label=long_label_for_sample(srr_id, metadata_row),
    )

    added_anything = False

    if include_mode == "minimal":
        for strand, negate in (("forward", False), ("reverse", True)):
            uniq = files.get((strand, "unique"))
            if uniq:
                sample_super.add_tracks(
                    _make_plain_track(srr_id, strand, "unique", uniq, base_url,
                                      colors, negate)
                )
                added_anything = True

    if include_mode in ("aggregate", "full"):
        agg = build_aggregate(srr_id, files, base_url, colors)
        if agg:
            sample_super.add_tracks(agg)
            added_anything = True

    if include_mode in ("composite", "full"):
        fw = build_kind_composite(srr_id, "forward", files, base_url, colors,
                                  negate=False)
        if fw:
            sample_super.add_tracks(fw)
            added_anything = True
        rev = build_kind_composite(srr_id, "reverse", files, base_url, colors,
                                   negate=True)
        if rev:
            sample_super.add_tracks(rev)
            added_anything = True

    return sample_super if added_anything else None


# ----- CLI -----------------------------------------------------------------

@click.group()
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
@click.pass_context
def cli(ctx: click.Context, verbose: bool) -> None:
    """RiboHub: build UCSC RiboSeq track hubs from sorted bigWig files."""
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
              help="Directory where the staged hub will be written.")
@click.option("--base-url", required=True, envvar="RIBOHUB_BASE_URL",
              help="Public URL where output-dir is served (must be reachable by UCSC).")
# Hub configuration
@click.option("--genome", default="hg38", show_default=True,
              help="UCSC genome assembly (e.g. hg38, mm10).")
@click.option("--hub-name", default="RiboSeqHub", show_default=True,
              help="Hub directory name (no spaces).")
@click.option("--email", default="your@email.com", show_default=True,
              help="Contact email written into hub.txt.")
# Build behavior
@click.option("--include", "include_mode", default="minimal", show_default=True,
              type=click.Choice(["minimal", "aggregate", "composite", "full"]),
              help="What to build per sample. minimal=just unique strands; "
                   "aggregate=overlay only; composite=kind comparisons only; "
                   "full=aggregate + both kind composites.")
@click.option("--strict", is_flag=True,
              help="Exit non-zero if any requested sample is missing or partial.")
@click.option("--dry-run", is_flag=True,
              help="Report what would be built without writing the hub.")
# Colors (hex format, e.g. #FF0000)
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
    genome: str,
    hub_name: str,
    email: str,
    include_mode: str,
    strict: bool,
    dry_run: bool,
    color_fwd: str,
    color_rev: str,
    color_fwd_multi: str,
    color_rev_multi: str,
) -> None:
    """Generate a UCSC track hub for the given samples."""
    base_url = base_url.rstrip("/")

    colors: dict[str, str] = {
        "fwd":       parse_hex_color(color_fwd),
        "rev":       parse_hex_color(color_rev),
        "fwd_multi": parse_hex_color(color_fwd_multi),
        "rev_multi": parse_hex_color(color_rev_multi),
    }

    sample_set = parse_samples(samples)
    log.info(f"Resolved {len(sample_set)} requested samples.")

    found, dir_missing, dir_empty = discover_samples(sample_set, data_dir)
    report_discovery(found, dir_missing, dir_empty, sample_set)

    if strict and (dir_missing or dir_empty):
        click.echo("--strict: refusing to build a partial hub.", err=True)
        sys.exit(1)
    if not found:
        click.echo("No samples found. Nothing to build.", err=True)
        sys.exit(1)

    if dry_run:
        click.echo(f"Dry run: would build {len(found)} sample(s) "
                   f"with --include {include_mode}.")
        return

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
            srr_id, found[srr_id], base_url, include_mode,
            colors=colors,
            metadata_row=None,
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

    staging = output_dir / hub_name
    if staging.exists():
        shutil.rmtree(staging)
    trackhub.upload.stage_hub(hub, staging=str(staging))

    click.echo(f"Built hub for {samples_built} sample(s) "
               f"with --include {include_mode}.")
    click.echo(f"Hub URL: {base_url}/{hub_name}/hub.txt")


if __name__ == "__main__":
    cli()
