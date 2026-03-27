"""Command line interface for the DepMap database application."""

import functools
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from .config import get_logger, get_settings, reload_settings
from .database import create_tables, get_current_schema_version, get_db_manager
from .database.connection import reset_db_manager
from .database.queries import get_available_views
from .downloader import FileManager, RefreshPlanner, ReleaseTracker
from .etl.processors import (
    GeneEffectWideProcessor,
    GeneExpressionWideProcessor,
    GeneProcessor,
    ModelProcessor,
)
from .etl.processors.base import BaseProcessor
from .utils.constants import (
    PHASE_1_DATASETS,
    SUPPORTED_DATASETS,
)

logger = get_logger(__name__)
console = Console()


def _get_processor(dataset_name: str) -> BaseProcessor | None:
    """Get the appropriate processor for a dataset.

    All data is stored in wide format (canonical storage).
    """
    processors: dict[str, type] = {
        "CRISPRGeneEffect": GeneEffectWideProcessor,
        "GeneExpression": GeneExpressionWideProcessor,
        "Gene": GeneProcessor,
        "Model": ModelProcessor,
    }

    processor_class = processors.get(dataset_name)
    return processor_class() if processor_class else None


def _detect_dataset_from_filename(filename: str) -> str | None:
    """Detect dataset type from filename."""
    filename_lower = filename.lower()

    if (
        "crisprgeneeffect" in filename_lower
        or "crispr_gene_effect" in filename_lower
    ):
        return "CRISPRGeneEffect"
    elif (
        "omicsexpression" in filename_lower
        or "expression" in filename_lower
        and "tpm" in filename_lower
    ):
        return "GeneExpression"
    elif (
        "model" in filename_lower and "data" in filename_lower
    ) or filename_lower.startswith("model"):
        return "Model"
    elif filename_lower.startswith("gene") and filename_lower.endswith(".csv"):
        return "Gene"
    elif (
        "modelcondition" in filename_lower
        or "model_condition" in filename_lower
    ):
        return "ModelCondition"

    return None


def handle_exceptions(func: Callable[..., Any]) -> Callable[..., Any]:
    """Decorator to handle and display exceptions gracefully."""

    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        try:
            return func(*args, **kwargs)
        except (OSError, ValueError, RuntimeError) as e:
            logger.error("Command failed: %s", e)
            console.print(f"[red]Error: {e}[/red]")
            sys.exit(1)

    return wrapper


def _resolve_requested_datasets(datasets: tuple[str, ...]) -> list[str]:
    """Resolve and validate requested dataset names."""
    if not datasets:
        return PHASE_1_DATASETS

    invalid = set(datasets) - set(SUPPORTED_DATASETS)
    if invalid:
        invalid_list = ", ".join(sorted(invalid))
        supported_list = ", ".join(SUPPORTED_DATASETS)
        raise ValueError(
            f"Invalid datasets: {invalid_list}. Available datasets: {supported_list}"
        )

    return list(datasets)


def _load_downloaded_files(
    downloaded_files: dict[str, Path], force: bool = False
) -> None:
    """Load downloaded dataset files into the database."""
    if not get_current_schema_version():
        console.print(
            "[red]Database not initialized. Run 'depmap-db init' first.[/red]"
        )
        return

    for dataset_name, file_path in downloaded_files.items():
        processor = _get_processor(dataset_name)
        if not processor:
            console.print(
                f"[yellow]Skipping {dataset_name} - processor not yet implemented[/yellow]"
            )
            continue

        console.print(f"[blue]Processing {dataset_name}...[/blue]")
        result = processor.process_file(file_path, force_reload=force)

        if result.status == "success":
            console.print(
                f"[green]✓ {dataset_name}: {result.records_inserted:,} records loaded[/green]"
            )
        elif result.status == "skipped":
            console.print(
                f"[yellow]- {dataset_name}: {result.records_skipped:,} records already exist[/yellow]"
            )
        else:
            console.print(
                f"[red]✗ {dataset_name}: {result.error_message}[/red]"
            )


@click.group()
@click.option(
    "--config-file",
    type=click.Path(exists=True, path_type=Path),
    help="Path to configuration file",
)
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
@click.option("--debug", is_flag=True, help="Enable debug mode")
@click.option(
    "--memory", is_flag=True, help="Use in-memory database for all operations"
)
@click.pass_context
def cli(
    ctx: click.Context,
    config_file: Path | None,
    verbose: bool,
    debug: bool,
    memory: bool,
) -> None:
    """DepMap Database - Build and query local DepMap databases."""
    # Ensure ctx.obj exists
    ctx.ensure_object(dict)

    # Store options in context
    ctx.obj["config_file"] = config_file
    ctx.obj["verbose"] = verbose
    ctx.obj["debug"] = debug
    ctx.obj["memory"] = memory

    # Configure global memory setting if requested
    if memory:
        import os

        os.environ["DEPMAP_DATABASE__MEMORY"] = "true"
        # Force reload of settings and reset database manager
        reload_settings()
        reset_db_manager()

    # Configure logging level if verbose or debug
    if debug:
        import logging

        logging.getLogger().setLevel(logging.DEBUG)
        logger.info("Debug mode enabled")
    elif verbose:
        import logging

        logging.getLogger().setLevel(logging.INFO)
        logger.info("Verbose mode enabled")


@cli.command()
@click.option(
    "--database-path",
    type=click.Path(path_type=Path),
    help="Path to database file",
)
@click.option("--memory", is_flag=True, help="Use in-memory database")
@handle_exceptions
def init(database_path: Path | None, memory: bool) -> None:
    """Initialize the database schema."""
    console.print("[blue]Initializing DepMap database...[/blue]")

    settings = get_settings()

    # Override database settings if provided
    if database_path:
        settings.database.path = database_path
    if memory:
        settings.database.memory = memory

    # Show configuration
    config_table = Table(title="Database Configuration")
    config_table.add_column("Setting", style="cyan")
    config_table.add_column("Value", style="green")

    if memory:
        config_table.add_row("Mode", "In-Memory")
    else:
        config_table.add_row("Database Path", str(settings.database.path))
    config_table.add_row("Max Memory", settings.database.max_memory)

    console.print(config_table)

    # Check if database already exists
    current_version = get_current_schema_version()
    if current_version:
        console.print(
            f"[yellow]Database already exists (version: {current_version})[/yellow]"
        )
        if not click.confirm(
            "Do you want to continue and recreate the schema?"
        ):
            console.print("[yellow]Initialization cancelled.[/yellow]")
            return

    # Create tables
    console.print("[blue]Creating database tables...[/blue]")
    create_tables()

    # Create views (temporarily disabled for wide format)
    # console.print("[blue]Creating predefined views...[/blue]")
    # create_all_views()

    console.print(
        "[green]✓ Database initialization completed successfully![/green]"
    )


@cli.command()
@handle_exceptions
def status() -> None:
    """Show database status and information."""
    settings = get_settings()
    db_manager = get_db_manager()

    # Database info panel
    if settings.database.memory:
        db_info = "[yellow]In-Memory Database[/yellow]"
    else:
        db_path = settings.database.path
        db_exists = db_path.exists() if db_path else False
        db_size = db_path.stat().st_size if db_exists else 0
        db_info = f"""Path: {db_path}
Exists: {"Yes" if db_exists else "No"}
Size: {db_size:,} bytes"""

    console.print(
        Panel(db_info, title="Database Information", title_align="left")
    )

    # Schema version
    schema_version = get_current_schema_version()
    if schema_version:
        console.print(f"[green]Schema Version: {schema_version}[/green]")
    else:
        console.print(
            "[red]No schema found - run 'depmap-db init' first[/red]"
        )
        return

    # Tables info
    try:
        # Get table names and row counts
        tables_query = """
        SELECT
            table_name,
            (SELECT COUNT(*) FROM information_schema.tables t2 WHERE t2.table_name = t1.table_name) as exists
        FROM information_schema.tables t1
        WHERE table_schema = 'main'
        ORDER BY table_name
        """

        result = db_manager.execute(tables_query)
        tables_data = result.fetchall()

        if tables_data:
            table_info = Table(title="Database Tables")
            table_info.add_column("Table", style="cyan")
            table_info.add_column("Row Count", style="green", justify="right")

            for table_name, _ in tables_data:
                try:
                    count_result = db_manager.execute(
                        f"SELECT COUNT(*) FROM {table_name}"
                    )
                    row_count = count_result.fetchone()[0]
                    table_info.add_row(table_name, f"{row_count:,}")
                except (OSError, RuntimeError):
                    table_info.add_row(table_name, "Error")

            console.print(table_info)
        else:
            console.print("[yellow]No tables found[/yellow]")

    except (OSError, RuntimeError) as e:
        console.print(f"[red]Error getting table information: {e}[/red]")


@cli.command()
@click.option(
    "--view",
    type=click.Choice(get_available_views()),
    help="Specific view to query",
)
@click.option("--limit", type=int, default=10, help="Number of rows to return")
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["table", "csv", "json"]),
    default="table",
    help="Output format",
)
@handle_exceptions
def query(view: str | None, limit: int, output_format: str) -> None:
    """Run predefined queries on the database."""
    db_manager = get_db_manager()

    # Check if database is initialized
    if not get_current_schema_version():
        console.print(
            "[red]Database not initialized. Run 'depmap-db init' first.[/red]"
        )
        return

    if not view:
        # Show available views
        console.print("[blue]Available predefined views:[/blue]")
        views_table = Table()
        views_table.add_column("View Name", style="cyan")
        views_table.add_column("Description", style="white")

        # TODO: Add view descriptions from ViewManager
        available_views = get_available_views()
        for view_name in available_views:
            views_table.add_row(view_name, "Predefined query view")

        console.print(views_table)
        return

    # Execute the view query
    console.print(f"[blue]Querying view: {view}[/blue]")

    query_sql = f"SELECT * FROM {view} LIMIT {limit}"

    try:
        if output_format == "table":
            df = db_manager.fetch_df(query_sql)

            if df.empty:
                console.print("[yellow]No results found[/yellow]")
                return

            # Create rich table
            result_table = Table(title=f"Results from {view} (limit {limit})")

            # Add columns
            for col in df.columns:
                result_table.add_column(str(col), style="cyan")

            # Add rows
            for _, row in df.head(limit).iterrows():
                result_table.add_row(*[str(val) for val in row])

            console.print(result_table)

        elif output_format == "csv":
            df = db_manager.fetch_df(query_sql)
            console.print(df.to_csv(index=False))

        elif output_format == "json":
            df = db_manager.fetch_df(query_sql)
            console.print(df.to_json(orient="records", indent=2))

    except (OSError, RuntimeError) as e:
        console.print(f"[red]Query failed: {e}[/red]")


@cli.command()
@click.option(
    "--datasets",
    multiple=True,
    help="Specific datasets to download (default: Phase 1 datasets)",
)
@click.option(
    "--cache-dir",
    type=click.Path(path_type=Path),
    help="Directory for caching downloads",
)
@click.option(
    "--force", is_flag=True, help="Force re-download even if files exist"
)
@click.option(
    "--load-data", is_flag=True, help="Also load downloaded data into database"
)
@handle_exceptions
def download(
    datasets: tuple[str, ...],
    cache_dir: Path | None,
    force: bool,
    load_data: bool,
) -> None:
    """Download DepMap data files."""
    from .downloader.depmap_client import download_datasets_sync

    requested_datasets = _resolve_requested_datasets(datasets)

    console.print(
        f"[blue]Downloading datasets: {', '.join(requested_datasets)}[/blue]"
    )

    # Set up file manager
    file_manager = FileManager(cache_dir)

    # Show cache status
    cache_stats = file_manager.get_cache_stats()
    if cache_stats["cached_datasets"] > 0:
        console.print(
            f"[green]Cache: {cache_stats['cached_datasets']} files ({cache_stats['total_size_formatted']})[/green]"
        )

    try:
        # Download files
        console.print("[blue]Starting downloads...[/blue]")
        downloaded_files = download_datasets_sync(
            requested_datasets, cache_dir
        )

        # Register downloads with file manager
        for dataset_name, file_path in downloaded_files.items():
            # Construct source URL (this is a simplified version)
            source_url = (
                f"https://depmap.org/portal/api/download/file/{file_path.name}"
            )
            file_manager.register_download(dataset_name, file_path, source_url)

        # Show results
        results_table = Table(title="Download Results")
        results_table.add_column("Dataset", style="cyan")
        results_table.add_column("File", style="white")
        results_table.add_column("Size", style="green", justify="right")

        for dataset_name, file_path in downloaded_files.items():
            file_size = file_path.stat().st_size
            from .utils.helpers import format_file_size

            results_table.add_row(
                dataset_name, file_path.name, format_file_size(file_size)
            )

        console.print(results_table)
        console.print(
            f"[green]✓ Successfully downloaded {len(downloaded_files)} files[/green]"
        )

        if load_data:
            console.print("[blue]Loading data into database...[/blue]")
            _load_downloaded_files(downloaded_files, force=force)
            console.print("[green]✓ Data loading completed[/green]")

    except (OSError, RuntimeError) as e:
        console.print(f"[red]Download failed: {e}[/red]")
        logger.error("Download failed: %s", e)
        raise


@cli.command()
@click.option(
    "--folder",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    required=True,
    help="Folder containing DepMap CSV files",
)
@click.option(
    "--force", is_flag=True, help="Force reload even if data already exists"
)
@click.option(
    "--pattern", default="*.csv", help="File pattern to match (default: *.csv)"
)
@handle_exceptions
def load_folder(folder: Path, force: bool, pattern: str) -> None:
    """Load all DepMap files from a folder into the database."""
    console.print(f"[blue]Loading DepMap files from: {folder}[/blue]")

    # Check database is initialized
    if not get_current_schema_version():
        console.print(
            "[red]Database not initialized. Run 'depmap-db init' first.[/red]"
        )
        return

    # Find CSV files in the folder
    csv_files = list(folder.glob(pattern))
    if not csv_files:
        console.print(
            f"[yellow]No files found matching pattern '{pattern}' in {folder}[/yellow]"
        )
        return

    console.print(f"[green]Found {len(csv_files)} files to process[/green]")

    # Process each file
    results = []
    processed_count = 0

    for file_path in csv_files:
        # Detect dataset type from filename
        dataset_name = _detect_dataset_from_filename(file_path.name)

        if not dataset_name:
            console.print(
                f"[yellow]Skipping {file_path.name} - cannot detect dataset type[/yellow]"
            )
            continue

        # Get processor
        processor = _get_processor(dataset_name)

        if not processor:
            console.print(
                f"[yellow]Skipping {file_path.name} - no processor for {dataset_name}[/yellow]"
            )
            continue

        # Show file info
        file_size = file_path.stat().st_size
        from .utils.helpers import format_file_size

        console.print(
            f"[cyan]Processing {file_path.name} → {dataset_name} ({format_file_size(file_size)})[/cyan]"
        )

        try:
            result = processor.process_file(file_path, force_reload=force)
            results.append((dataset_name, result))

            if result.status == "success":
                console.print(
                    f"  [green]✓ {result.records_inserted:,} records loaded[/green]"
                )
                processed_count += 1
            elif result.status == "skipped":
                console.print(
                    f"  [yellow]- {result.records_skipped:,} records already exist (use --force to reload)[/yellow]"
                )
            else:
                console.print(f"  [red]✗ Failed: {result.error_message}[/red]")

            if result.warnings:
                console.print(
                    f"  [yellow]Warnings: {len(result.warnings)}[/yellow]"
                )

        except (OSError, ValueError, RuntimeError) as e:
            console.print(f"  [red]✗ Error: {e}[/red]")
            logger.error("Failed to process %s: %s", file_path, e)

    # Show summary
    console.print("\n[blue]Summary:[/blue]")
    console.print(f"  Files found: {len(csv_files)}")
    console.print(f"  Files processed: {processed_count}")

    if results:
        summary_table = Table(title="Loading Results")
        summary_table.add_column("Dataset", style="cyan")
        summary_table.add_column("Status", style="white")
        summary_table.add_column("Records", style="green", justify="right")
        summary_table.add_column("Time", style="blue", justify="right")

        for dataset_name, result in results:
            status_color = (
                "green"
                if result.status == "success"
                else "yellow"
                if result.status == "skipped"
                else "red"
            )
            status_text = (
                f"[{status_color}]{result.status.title()}[/{status_color}]"
            )

            records_text = (
                f"{result.records_inserted:,}"
                if result.status == "success"
                else f"{result.records_skipped:,}"
            )
            time_text = f"{result.processing_time_seconds:.1f}s"

            summary_table.add_row(
                dataset_name, status_text, records_text, time_text
            )

        console.print(summary_table)

    # Show database status
    console.print("\n[green]✓ Folder loading completed[/green]")


@cli.command()
@click.option(
    "--datasets",
    multiple=True,
    help="Specific datasets to refresh (default: Phase 1 datasets)",
)
@click.option(
    "--cache-dir",
    type=click.Path(path_type=Path),
    help="Directory for cached downloads and release tracking",
)
@click.option(
    "--apply/--plan-only",
    default=False,
    help="Apply the refresh plan instead of only showing it",
)
@click.option(
    "--force", is_flag=True, help="Force reload into the database when applying"
)
@click.option(
    "--load-data",
    is_flag=True,
    help="Also load refreshed datasets into the database when applying",
)
@handle_exceptions
def refresh(
    datasets: tuple[str, ...],
    cache_dir: Path | None,
    apply: bool,
    force: bool,
    load_data: bool,
) -> None:
    """Plan or apply a refresh for the configured DepMap release."""
    requested_datasets = _resolve_requested_datasets(datasets)
    file_manager = FileManager(cache_dir)
    tracker = None
    if cache_dir:
        tracker = ReleaseTracker(cache_dir / "release_state.json")
    planner = RefreshPlanner(file_manager=file_manager, tracker=tracker)
    plan = planner.build_plan(requested_datasets)

    console.print(
        f"[blue]Refresh plan for release '{plan.snapshot.release_label}'[/blue]"
    )
    console.print(f"[cyan]Reason:[/cyan] {plan.reason}")
    console.print(
        f"[cyan]Cached:[/cyan] {', '.join(plan.cached_datasets) if plan.cached_datasets else 'none'}"
    )
    console.print(
        f"[cyan]To download:[/cyan] {', '.join(plan.datasets_to_download) if plan.datasets_to_download else 'none'}"
    )

    if not apply:
        console.print(
            "[yellow]Plan only. Re-run with --apply to download and record this release.[/yellow]"
        )
        return

    downloaded_files: dict[str, Path] = {}
    if plan.datasets_to_download:
        from .downloader.depmap_client import download_datasets_sync

        console.print("[blue]Downloading refresh datasets...[/blue]")
        downloaded_files = download_datasets_sync(
            plan.datasets_to_download, cache_dir
        )

        for dataset_name, file_path in downloaded_files.items():
            source_url = (
                f"https://depmap.org/portal/api/download/file/{file_path.name}"
            )
            file_manager.register_download(dataset_name, file_path, source_url)
    else:
        console.print("[green]No downloads needed for this release.[/green]")

    planner.mark_applied(plan.snapshot)
    console.print("[green]✓ Recorded applied release state[/green]")

    if load_data and downloaded_files:
        console.print("[blue]Loading refreshed data into database...[/blue]")
        _load_downloaded_files(downloaded_files, force=force)
        console.print("[green]✓ Refresh data loading completed[/green]")


@cli.command()
@click.option("--sql", help="Custom SQL query to execute")
@click.option(
    "--file",
    "sql_file",
    type=click.Path(exists=True, path_type=Path),
    help="File containing SQL query",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["table", "csv", "json"]),
    default="table",
    help="Output format",
)
@handle_exceptions
def sql(sql: str | None, sql_file: Path | None, output_format: str) -> None:
    """Execute custom SQL queries."""
    db_manager = get_db_manager()

    # Check if database is initialized
    if not get_current_schema_version():
        console.print(
            "[red]Database not initialized. Run 'depmap-db init' first.[/red]"
        )
        return

    # Get SQL query
    if sql_file:
        query_sql = sql_file.read_text()
    elif sql:
        query_sql = sql
    else:
        console.print(
            "[red]Please provide either --sql or --file option[/red]"
        )
        return

    console.print("[blue]Executing SQL query...[/blue]")

    try:
        if output_format == "table":
            df = db_manager.fetch_df(query_sql)

            if df.empty:
                console.print("[yellow]No results returned[/yellow]")
                return

            # Create rich table
            result_table = Table(title="Query Results")

            # Add columns
            for col in df.columns:
                result_table.add_column(str(col), style="cyan")

            # Add rows (limit to first 50 for display)
            for _, row in df.head(50).iterrows():
                result_table.add_row(*[str(val) for val in row])

            if len(df) > 50:
                console.print(
                    f"[yellow]Showing first 50 of {len(df)} results[/yellow]"
                )

            console.print(result_table)

        elif output_format == "csv":
            df = db_manager.fetch_df(query_sql)
            console.print(df.to_csv(index=False))

        elif output_format == "json":
            df = db_manager.fetch_df(query_sql)
            console.print(df.to_json(orient="records", indent=2))

    except (OSError, RuntimeError) as e:
        console.print(f"[red]Query failed: {e}[/red]")


@cli.command()
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    required=True,
    help="Output CSV file path",
)
@click.option(
    "--data-type",
    type=click.Choice(
        ["expression", "crispr", "integrated", "matched", "joint"]
    ),
    required=True,
    help="Type of data to export: expression, crispr, integrated, matched, or joint (perfectly matched CRISPR + expression)",
)
@click.option(
    "--filter",
    "models_filter",
    help="SQL WHERE clause to filter models (e.g., \"oncotree_lineage = 'Lung'\")",
)
@click.option(
    "--genes",
    help="Comma-separated list of genes to include (for integrated export only)",
)
@handle_exceptions
def export(
    output: Path,
    data_type: str,
    models_filter: str | None,
    genes: str | None,
) -> None:
    """Export data from database to CSV file.

    Examples:

      # Export all gene expression data
      depmap-db export -o expression.csv --data-type expression

      # Export CRISPR data for lung cancer only
      depmap-db export -o lung_crispr.csv --data-type crispr --filter "oncotree_lineage = 'Lung'"

      # Export integrated expression + dependency for specific genes
      depmap-db export -o tp53_myc.csv --data-type integrated --genes TP53,MYC

      # Export expression data matched to CRISPR (only models and genes in CRISPR)
      depmap-db export -o matched_expression.csv --data-type matched

      # Export perfectly matched CRISPR + expression for correlation analysis
      depmap-db export -o matched --data-type joint
    """
    from .export import (
        export_gene_effects_csv,
        export_gene_expression_csv,
        export_integrated_data_csv,
        export_joint_crispr_expression_csv,
        export_matched_expression_csv,
    )

    # Check if database is initialized
    if not get_current_schema_version():
        console.print(
            "[red]Database not initialized. Run 'depmap-db init' first.[/red]"
        )
        return

    console.print(f"[blue]Exporting {data_type} data...[/blue]")

    try:
        if data_type == "expression":
            num_models = export_gene_expression_csv(
                output_path=output, models_filter=models_filter
            )
            console.print(
                f"[green]✓ Exported {num_models:,} models to {output}[/green]"
            )

        elif data_type == "crispr":
            num_models = export_gene_effects_csv(
                output_path=output, models_filter=models_filter
            )
            console.print(
                f"[green]✓ Exported {num_models:,} models to {output}[/green]"
            )

        elif data_type == "integrated":
            genes_list = None
            if genes:
                genes_list = [g.strip() for g in genes.split(",")]
                console.print(
                    f"[cyan]Including genes: {', '.join(genes_list)}[/cyan]"
                )

            num_models = export_integrated_data_csv(
                output_path=output,
                genes_to_include=genes_list,
                models_filter=models_filter,
            )
            console.print(
                f"[green]✓ Exported {num_models:,} models to {output}[/green]"
            )

        elif data_type == "matched":
            console.print(
                "[cyan]Matching expression data to CRISPR data...[/cyan]"
            )
            result = export_matched_expression_csv(
                output_path=output,
                models_filter=models_filter,
                report_missing=True,
            )

            if result["models_exported"] > 0:
                console.print(
                    f"[green]✓ Exported {result['models_exported']:,} models with {result['genes_exported']:,} genes to {output}[/green]"
                )
                console.print(
                    f"[cyan]  Models in both datasets: {result['models_in_both']:,}/{result['total_crispr_models']:,} CRISPR models[/cyan]"
                )
                console.print(
                    f"[cyan]  Genes in both datasets: {result['genes_in_both']:,}/{result['total_crispr_genes']:,} CRISPR genes[/cyan]"
                )

                if result["models_in_crispr_not_expression"] > 0:
                    console.print(
                        f"[yellow]  ⚠️  {result['models_in_crispr_not_expression']:,} CRISPR models missing in expression data[/yellow]"
                    )
                if result["genes_in_crispr_not_expression"] > 0:
                    console.print(
                        f"[yellow]  ⚠️  {result['genes_in_crispr_not_expression']:,} CRISPR genes missing in expression data[/yellow]"
                    )
            else:
                console.print(
                    "[red]✗ No overlapping data found between CRISPR and expression[/red]"
                )
                return

        elif data_type == "joint":
            # For joint export, we need two output files
            # Generate filenames based on the provided output path
            output_stem = output.stem
            output_dir = output.parent

            expression_output = output_dir / f"{output_stem}_expression.csv"
            crispr_output = output_dir / f"{output_stem}_crispr.csv"

            console.print(
                "[cyan]Exporting perfectly matched CRISPR + expression datasets...[/cyan]"
            )
            console.print(f"[cyan]  Expression → {expression_output}[/cyan]")
            console.print(f"[cyan]  CRISPR     → {crispr_output}[/cyan]")

            result = export_joint_crispr_expression_csv(
                expression_output_path=expression_output,
                crispr_output_path=crispr_output,
                models_filter=models_filter,
                report_missing=True,
            )

            if result["models_exported"] > 0:
                console.print(
                    "\n[green]✓ Successfully exported perfectly matched datasets![/green]"
                )
                console.print(
                    f"[green]  {result['models_exported']:,} models × {result['genes_exported']:,} genes[/green]"
                )
                console.print(
                    f"[cyan]  Expression: {expression_output.name} ({result['expression_file_size_bytes']:,} bytes)[/cyan]"
                )
                console.print(
                    f"[cyan]  CRISPR:     {crispr_output.name} ({result['crispr_file_size_bytes']:,} bytes)[/cyan]"
                )

                # Show what was excluded
                excluded = []
                if result["models_in_crispr_not_expression"] > 0:
                    excluded.append(
                        f"{result['models_in_crispr_not_expression']:,} CRISPR-only models"
                    )
                if result["genes_in_crispr_not_expression"] > 0:
                    excluded.append(
                        f"{result['genes_in_crispr_not_expression']:,} CRISPR-only genes"
                    )
                if result["models_in_expression_not_crispr"] > 0:
                    excluded.append(
                        f"{result['models_in_expression_not_crispr']:,} expression-only models"
                    )
                if result["genes_in_expression_not_crispr"] > 0:
                    excluded.append(
                        f"{result['genes_in_expression_not_crispr']:,} expression-only genes"
                    )

                if excluded:
                    console.print(
                        f"\n[yellow]  Excluded: {', '.join(excluded)}[/yellow]"
                    )

                console.print(
                    "\n[bold green]Both files contain identical models and genes - ready for correlation analysis![/bold green]"
                )
            else:
                console.print(
                    "[red]✗ No overlapping data found between CRISPR and expression[/red]"
                )
                return

        # Show file info (skip for joint export since we show it inline)
        if data_type != "joint":
            file_size = output.stat().st_size
            from .utils.helpers import format_file_size

            console.print(
                f"[cyan]File size: {format_file_size(file_size)}[/cyan]"
            )

    except (OSError, RuntimeError) as e:
        console.print(f"[red]Export failed: {e}[/red]")
        logger.error("Export failed: %s", e)
        raise


if __name__ == "__main__":
    cli()
