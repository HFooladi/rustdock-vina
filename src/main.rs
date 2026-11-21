//! Main executable for rustdock-vina

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};
use log::{info, warn};
use nalgebra::Vector3;
use rayon::prelude::*;
use rustdock_vina::forcefield::{
    ad4::AD4ForceField, vina::VinaForceField, vinardo::VinardoForceField, ForceField,
};
use rustdock_vina::io::{parse_flex_pdbqt, parse_pdbqt, write_docking_results, write_pose};
use rustdock_vina::optimization::monte_carlo::{MonteCarlo, MonteCarloParams, Target};
use rustdock_vina::optimization::DockingResult;
use rustdock_vina::scoring::{ligand_internal_energy, score_pose_with_flex};
use rustdock_vina::{FlexResidues, Molecule};
use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

/// Command-line arguments for the application
#[derive(Parser, Debug)]
#[clap(
    name = "rustdock",
    version = rustdock_vina::VERSION,
    about = "A Rust implementation of AutoDock Vina for molecular docking"
)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

// Parsed once per run, so the size of `Dock` does not matter.
#[allow(clippy::large_enum_variant)]
#[derive(Subcommand, Debug)]
enum Commands {
    /// Dock a ligand to a receptor
    Dock(DockArgs),

    /// Score an existing pose
    Score {
        /// PDBQT file containing the receptor
        #[clap(long, value_parser)]
        receptor: PathBuf,

        /// PDBQT file containing the ligand pose to score
        #[clap(long, value_parser)]
        ligand: PathBuf,

        /// Flexible side chains (PDBQT of BEGIN_RES/END_RES blocks, as for
        /// Vina's --flex); --receptor is then the rigid part
        #[clap(long, value_parser)]
        flex: Option<PathBuf>,

        /// Scoring function to use (vina, ad4, vinardo)
        #[clap(long, default_value = "vina")]
        scoring: String,
    },
}

/// Options of `rustdock dock`. Those with a default are `Option`s so that
/// one left off the command line can come from `--config`.
#[derive(Args, Debug)]
struct DockArgs {
    /// PDBQT file containing the receptor
    #[clap(long, value_parser)]
    receptor: Option<PathBuf>,

    /// PDBQT file(s) containing the ligand to dock; each is docked on its
    /// own
    #[clap(long, value_parser, num_args = 1..)]
    ligand: Vec<PathBuf>,

    /// Flexible side chains (PDBQT of BEGIN_RES/END_RES blocks, as for
    /// Vina's --flex); --receptor is then the rigid part. Their torsions
    /// are searched with the ligand's and written after it in each mode
    #[clap(long, value_parser)]
    flex: Option<PathBuf>,

    /// Vina configuration file: `key = value` lines naming any option
    /// in Vina's spelling. Options on the command line take precedence
    #[clap(long, short, value_parser)]
    config: Option<PathBuf>,

    /// Center of the search box (x,y,z)
    #[clap(
        long,
        value_parser,
        value_delimiter = ',',
        allow_negative_numbers = true
    )]
    center: Option<Vec<f64>>,

    /// Center of the search box, X (Angstrom)
    #[clap(
        long,
        visible_alias = "center_x",
        allow_negative_numbers = true,
        conflicts_with = "center"
    )]
    center_x: Option<f64>,

    /// Center of the search box, Y (Angstrom)
    #[clap(
        long,
        visible_alias = "center_y",
        allow_negative_numbers = true,
        conflicts_with = "center"
    )]
    center_y: Option<f64>,

    /// Center of the search box, Z (Angstrom)
    #[clap(
        long,
        visible_alias = "center_z",
        allow_negative_numbers = true,
        conflicts_with = "center"
    )]
    center_z: Option<f64>,

    /// Size of the search box (x,y,z)
    #[clap(long, value_parser, value_delimiter = ',')]
    size: Option<Vec<f64>>,

    /// Size of the search box, X (Angstrom)
    #[clap(long, visible_alias = "size_x", conflicts_with = "size")]
    size_x: Option<f64>,

    /// Size of the search box, Y (Angstrom)
    #[clap(long, visible_alias = "size_y", conflicts_with = "size")]
    size_y: Option<f64>,

    /// Size of the search box, Z (Angstrom)
    #[clap(long, visible_alias = "size_z", conflicts_with = "size")]
    size_z: Option<f64>,

    /// Scoring function to use (vina, ad4, vinardo) [default: vina]
    #[clap(long)]
    scoring: Option<String>,

    /// Exhaustiveness of the search (higher values increase accuracy but
    /// take longer) [default: 8]
    #[clap(long)]
    exhaustiveness: Option<usize>,

    /// Number of binding modes to generate [default: 9]
    #[clap(long, visible_alias = "num_modes")]
    num_modes: Option<usize>,

    /// Minimum heavy-atom RMSD between output modes (Angstrom) [default: 1]
    #[clap(long, visible_alias = "min_rmsd")]
    min_rmsd: Option<f64>,

    /// Output file for the docking results
    #[clap(long, short, value_parser)]
    out: Option<PathBuf>,

    /// Ligands to screen, as PDBQT files or directories of them. They are
    /// docked in parallel against one set of receptor grids; each one's
    /// poses go to `--dir`, and `--dir`/scores.tsv lists every best
    /// affinity. A ligand that fails is reported there and skipped.
    #[clap(long, value_parser, num_args = 1..)]
    batch: Vec<PathBuf>,

    /// Output directory for the docking results (required with --batch)
    #[clap(long, value_parser)]
    dir: Option<PathBuf>,

    /// Number of threads; 0 uses all cores [default: 0]
    #[clap(long)]
    cpu: Option<usize>,

    /// Energy range for output poses (kcal/mol) [default: 3]
    #[clap(long, visible_alias = "energy_range")]
    energy_range: Option<f64>,

    /// Spacing of the precomputed affinity grids (Angstroms) [default: 0.375]
    #[clap(long, visible_alias = "spacing")]
    grid_spacing: Option<f64>,

    /// Search on exact pairwise energies instead of precomputed grids (slow)
    #[clap(long)]
    no_grid: bool,

    /// Score the ligand where it stands and print the breakdown, without
    /// searching; the search box may be left out
    #[clap(long, visible_alias = "score_only")]
    score_only: bool,

    /// Minimize the ligand from where it stands and write the one pose,
    /// without searching
    #[clap(long, visible_alias = "local_only")]
    local_only: bool,

    /// Take the search box from the ligand's own extent, with Vina's 4 A of
    /// padding (only with --score_only or --local_only)
    #[clap(long)]
    autobox: bool,

    /// Random seed; the same seed and inputs reproduce the same poses.
    /// 0 draws a random seed and prints it, as Vina does [default: 0]
    #[clap(long)]
    seed: Option<u64>,
}

fn main() -> Result<()> {
    // Default to `info` so a run is not silent. Previously every progress
    // message was suppressed unless the user happened to know to set RUST_LOG,
    // which is documented nowhere.
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp(None)
        .init();

    // Parse command-line arguments
    let cli = Cli::parse();

    match cli.command {
        Commands::Dock(args) => {
            let config = match &args.config {
                Some(path) => {
                    let text = std::fs::read_to_string(path).with_context(|| {
                        format!("Failed to read config file: {}", path.display())
                    })?;
                    parse_config(&text).with_context(|| {
                        format!("Failed to parse config file: {}", path.display())
                    })?
                }
                None => VinaConfig::default(),
            };
            let DockSettings {
                receptor,
                ligands,
                batch,
                flex,
                out,
                dir,
                scoring,
                center,
                size,
                score_only,
                local_only,
                autobox,
                exhaustiveness,
                num_modes,
                min_rmsd,
                energy_range,
                grid_spacing,
                no_grid,
                cpu,
                seed,
            } = resolve(args, config)?;

            // Determine scoring function
            let forcefield = select_forcefield(&scoring)?;

            info!("Using {} scoring function", forcefield.name());

            if let Some(threads) = cpu {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build_global()
                    .context("failed to set the thread count")?;
            }

            let batch = expand_batch(&batch)?;
            if ligands.is_empty() && batch.is_empty() {
                return Err(anyhow::anyhow!("--batch found no .pdbqt ligands"));
            }

            // Parse receptor
            info!("Loading receptor: {}", receptor.display());
            let receptor_molecule = parse_pdbqt(&receptor).with_context(|| {
                format!("Failed to parse receptor file: {}", receptor.display())
            })?;

            // Scoring a pose needs no search, no grids and no box.
            if score_only {
                let flex_residues = flex
                    .as_ref()
                    .map(|path| load_flex(path, &receptor_molecule))
                    .transpose()?;
                if flex_residues.is_some() {
                    // Vina's own --score_only takes the unbound system's
                    // energy from the grids while reporting the bound
                    // internal energy exactly, so the two do not cancel and
                    // its affinity sits a little apart from this one. Every
                    // term below still matches.
                    warn!(
                        "with --flex, Vina's --score_only affinity differs from this one: it \
                         takes the unbound system's energy from its grids, so (2) and (4) do \
                         not cancel. Pass Vina --unbound_energy with the internal energy it \
                         prints to compare them."
                    );
                }
                for ligand_path in ligands.iter().chain(&batch) {
                    info!("Loading ligand: {}", ligand_path.display());
                    let ligand = parse_pdbqt(ligand_path).with_context(|| {
                        format!("Failed to parse ligand file: {}", ligand_path.display())
                    })?;
                    if autobox {
                        let (center, size) = autobox_of(&ligand);
                        info!("Grid center: ligand center (autobox)");
                        info!(
                            "Grid size  : X {} Y {} Z {} (autobox)",
                            size.x, size.y, size.z
                        );
                        let _ = center;
                    }
                    if ligands.len() + batch.len() > 1 {
                        println!("\n{}", ligand_path.display());
                    }
                    show_score(
                        &ligand,
                        &receptor_molecule,
                        flex_residues.as_ref(),
                        forcefield.as_ref(),
                    )?;
                }
                return Ok(());
            }

            // Minimizing or docking needs a box. `resolve` has insisted on
            // one unless --autobox is to derive it from the ligand.
            let (center, box_size) = match (center, size) {
                (Some(center), Some(size)) => (center, size),
                _ if autobox => {
                    let first = ligands.first().or_else(|| batch.first()).ok_or_else(|| {
                        anyhow::anyhow!("--autobox needs a ligand to take the box from")
                    })?;
                    let ligand = parse_pdbqt(first).with_context(|| {
                        format!("Failed to parse ligand file: {}", first.display())
                    })?;
                    let (center, size) = autobox_of(&ligand);
                    info!("Grid center: ligand center (autobox)");
                    info!(
                        "Grid size  : X {} Y {} Z {} (autobox)",
                        size.x, size.y, size.z
                    );
                    (center, size)
                }
                _ => return Err(anyhow::anyhow!("docking needs --center and --size")),
            };

            // As Vina, draw a seed when none is given and print it, so that
            // any run can be reproduced.
            let seed = seed.unwrap_or_else(|| rand::random::<u64>().max(1));
            info!("Seed: {seed}");

            info!(
                "Search box center: ({}, {}, {})",
                center.x, center.y, center.z
            );
            info!(
                "Search box size: ({}, {}, {})",
                box_size.x, box_size.y, box_size.z
            );

            // Create optimizer
            if !no_grid && !(grid_spacing.is_finite() && grid_spacing > 0.0) {
                return Err(anyhow::anyhow!(
                    "--grid-spacing must be a positive number, got {grid_spacing}"
                ));
            }
            let optimizer = MonteCarlo::with_params(MonteCarloParams {
                exhaustiveness,
                use_grid: !no_grid,
                grid_spacing,
                seed: Some(seed),
                min_rmsd,
                ..MonteCarloParams::default()
            });
            // Every ligand shares one set of receptor grids.
            let target = match &flex {
                Some(path) => {
                    info!("Loading flexible residues: {}", path.display());
                    let flex_residues =
                        parse_flex_pdbqt(path, &receptor_molecule).with_context(|| {
                            format!("Failed to parse flex residue file: {}", path.display())
                        })?;
                    warn_flex_outside_box(&flex_residues, &center, &box_size);
                    optimizer.target_with_flex(
                        receptor_molecule,
                        flex_residues,
                        forcefield.as_ref(),
                        center,
                        box_size,
                    )?
                }
                None => {
                    optimizer.target(receptor_molecule, forcefield.as_ref(), center, box_size)?
                }
            };
            // Minimizing takes the pose as it is, so there is no search and
            // only ever one mode to write.
            if local_only {
                for ligand_path in ligands.iter().chain(&batch) {
                    info!("Loading ligand: {}", ligand_path.display());
                    let ligand = parse_pdbqt(ligand_path).with_context(|| {
                        format!("Failed to parse ligand file: {}", ligand_path.display())
                    })?;
                    info!("Minimizing {}", ligand_path.display());
                    let local = optimizer.local_optimize(&target, &ligand, forcefield.as_ref())?;
                    let output_path = match (&out, &dir) {
                        (Some(out_path), _) => out_path.clone(),
                        (None, Some(dir_path)) => {
                            let stem = ligand_path
                                .file_stem()
                                .unwrap_or_default()
                                .to_string_lossy();
                            dir_path.join(format!("{stem}_out.pdbqt"))
                        }
                        (None, None) => default_output(ligand_path),
                    };
                    if ligands.len() + batch.len() > 1 {
                        println!("\n{}", ligand_path.display());
                    }
                    show_score(
                        &local.pose.molecule,
                        target.receptor(),
                        local.flex.as_ref(),
                        forcefield.as_ref(),
                    )?;
                    info!("Writing the minimized pose to {}", output_path.display());
                    // As Vina's own --local_only: a plain PDBQT, which Vina
                    // can read back (it rejects MODEL tags on one pose).
                    write_pose(&local.pose, &output_path)
                        .with_context(|| format!("Failed to write {}", output_path.display()))?;
                }
                return Ok(());
            }

            let job = DockJob {
                optimizer: &optimizer,
                target: &target,
                forcefield: forcefield.as_ref(),
                num_modes,
                energy_range,
            };

            // Process each ligand
            for ligand_path in ligands {
                // Output path: --out, else --dir/<stem>_out.pdbqt, else
                // Vina's default beside the ligand.
                let output_path = match (&out, &dir) {
                    (Some(out_path), _) => out_path.clone(),
                    (None, Some(dir_path)) => {
                        let stem = ligand_path
                            .file_stem()
                            .unwrap_or_default()
                            .to_string_lossy();
                        dir_path.join(format!("{stem}_out.pdbqt"))
                    }
                    (None, None) => default_output(&ligand_path),
                };

                info!(
                    "Docking ligand {} with exhaustiveness {}",
                    ligand_path.display(),
                    exhaustiveness
                );
                let filtered_results = job.dock(&ligand_path)?;

                // Print the familiar Vina mode/affinity/RMSD table. Without
                // this the command produced no terminal output at all and the
                // user had to open the output file to learn anything.
                println!("\nmode |   affinity | rmsd l.b.| rmsd u.b.");
                println!("     | (kcal/mol) | (A)      | (A)");
                println!("-----+------------+----------+----------");
                for (i, r) in filtered_results.iter().enumerate() {
                    println!(
                        "{:4} | {:10.3} | {:8.3} | {:8.3}",
                        i + 1,
                        r.energy,
                        r.rmsd_lb.unwrap_or(0.0),
                        r.rmsd_ub.unwrap_or(0.0)
                    );
                }
                println!();

                // Write results
                info!(
                    "Writing {} docking poses to {}",
                    filtered_results.len(),
                    output_path.display()
                );
                write_docking_results(&filtered_results, &output_path).with_context(|| {
                    format!(
                        "Failed to write docking results to {}",
                        output_path.display()
                    )
                })?;
            }

            if let Some(dir) = dir.as_ref().filter(|_| !batch.is_empty()) {
                run_batch(&job, &batch, dir)?;
            }

            info!("Docking completed successfully");
        }

        Commands::Score {
            receptor,
            ligand,
            flex,
            scoring,
        } => {
            // Determine scoring function
            let forcefield = select_forcefield(&scoring)?;

            info!("Using {} scoring function", forcefield.name());

            // Parse receptor
            info!("Loading receptor: {}", receptor.display());
            let receptor_molecule = parse_pdbqt(&receptor).with_context(|| {
                format!("Failed to parse receptor file: {}", receptor.display())
            })?;

            // Parse ligand
            info!("Loading ligand: {}", ligand.display());
            let ligand_molecule = parse_pdbqt(&ligand)
                .with_context(|| format!("Failed to parse ligand file: {}", ligand.display()))?;
            let flex_residues = flex
                .as_ref()
                .map(|path| {
                    info!("Loading flexible residues: {}", path.display());
                    parse_flex_pdbqt(path, &receptor_molecule).with_context(|| {
                        format!("Failed to parse flex residue file: {}", path.display())
                    })
                })
                .transpose()?;

            println!("\nLigand atoms: {}", ligand_molecule.atoms.len());
            println!("Ligand bonds: {}", ligand_molecule.bonds.len());
            println!("Ligand torsions: {}", ligand_molecule.torsions.len());
            println!("Receptor atoms: {}", receptor_molecule.atoms.len());
            if let Some(flex) = &flex_residues {
                println!(
                    "Flexible residues: {} ({} atoms, {} torsions)",
                    flex.num_residues(),
                    flex.molecule.atoms.len(),
                    flex.molecule.torsions.len()
                );
            }

            // Vina formula: affinity = intermolecular / (1 + weight_rot * num_tors)
            // The internal energy cancels with unbound energy in the full formula
            let score = score_pose_with_flex(
                &ligand_molecule,
                &receptor_molecule,
                flex_residues.as_ref(),
                forcefield.as_ref(),
            )?;
            let torsional_effect = score.intermolecular - score.affinity;

            println!("\nScoring breakdown:");
            println!("  Atom pairs within 8A: {}", score.pairs);
            println!("  Positive energy (repulsion): {:.3}", score.positive);
            println!("  Negative energy (attraction): {:.3}", score.negative);
            println!("  Intermolecular energy: {:.3}", score.intermolecular);
            if flex_residues.is_some() {
                println!(
                    "    Ligand - Receptor: {:.3}",
                    score.intermolecular - score.ligand_flex
                );
                println!("    Ligand - Flex side chains: {:.3}", score.ligand_flex);
                println!(
                    "  Flex - Receptor (internal, not in the affinity): {:.3}",
                    score.flex_receptor
                );
                println!(
                    "  Flex - Flex side chains (internal, not in the affinity): {:.3}",
                    score.flex_flex
                );
            }
            println!(
                "  Torsional effect: {:.3} ({} torsions, weighted count {})",
                torsional_effect,
                ligand_molecule.torsions.len(),
                score.num_tors
            );
            println!("\nAffinity: {:.3} (kcal/mol)", score.affinity);
            info!("Scoring completed successfully");
        }
    }

    Ok(())
}

/// Options read from a Vina configuration file. Every field is optional:
/// the command line takes precedence, then this file, then Vina's defaults.
#[derive(Debug, Default, PartialEq)]
struct VinaConfig {
    receptor: Option<PathBuf>,
    ligand: Vec<PathBuf>,
    batch: Vec<PathBuf>,
    flex: Option<PathBuf>,
    out: Option<PathBuf>,
    dir: Option<PathBuf>,
    scoring: Option<String>,
    center: [Option<f64>; 3],
    size: [Option<f64>; 3],
    exhaustiveness: Option<usize>,
    num_modes: Option<usize>,
    min_rmsd: Option<f64>,
    energy_range: Option<f64>,
    spacing: Option<f64>,
    cpu: Option<usize>,
    seed: Option<i64>,
    score_only: Option<bool>,
    local_only: Option<bool>,
    autobox: Option<bool>,
}

/// Vina options rustdock does not implement, with Vina 1.2.7's default.
/// A config file may set one to its default; any other value is an error,
/// because ignoring it would change the result.
const UNSUPPORTED_WITH_DEFAULT: &[(&str, f64)] = &[
    ("max_evals", 0.0),
    ("weight_gauss1", -0.035579),
    ("weight_gauss2", -0.005156),
    ("weight_repulsion", 0.840245),
    ("weight_hydrophobic", -0.035069),
    ("weight_hydrogen", -0.587439),
    ("weight_rot", 0.05846),
    ("weight_vinardo_gauss1", -0.045),
    ("weight_vinardo_repulsion", 0.8),
    ("weight_vinardo_hydrophobic", -0.035),
    ("weight_vinardo_hydrogen", -0.6),
    ("weight_vinardo_rot", 0.05846),
    ("weight_ad4_vdw", 0.1662),
    ("weight_ad4_hb", 0.1209),
    ("weight_ad4_elec", 0.1406),
    ("weight_ad4_dsolv", 0.1322),
    ("weight_ad4_rot", 0.2983),
    ("weight_glue", 50.0),
];

/// Vina switches rustdock does not implement; `false` is accepted.
const UNSUPPORTED_SWITCHES: &[&str] = &["randomize_only", "no_refine", "force_even_voxels"];

/// Read a Vina configuration file as Vina 1.2.7 (Boost.Program_options)
/// does: `key = value` lines, `#` to the end of a line is a comment, an
/// unknown key or a key given twice is an error, and `ligand` / `batch`
/// lines accumulate.
fn parse_config(config_str: &str) -> Result<VinaConfig> {
    fn num<T: std::str::FromStr>(key: &str, value: &str) -> Result<T> {
        value
            .parse()
            .map_err(|_| anyhow::anyhow!("{key}: {value:?} is not a valid number"))
    }
    fn switch(key: &str, value: &str) -> Result<bool> {
        match value.to_ascii_lowercase().as_str() {
            "" | "1" | "true" | "yes" | "on" => Ok(true),
            "0" | "false" | "no" | "off" => Ok(false),
            _ => Err(anyhow::anyhow!("{key}: {value:?} is not true or false")),
        }
    }

    let mut config = VinaConfig::default();
    let mut seen = HashSet::new();
    for (index, raw) in config_str.lines().enumerate() {
        let line = raw.split('#').next().unwrap_or_default().trim();
        if line.is_empty() {
            continue;
        }
        let at = || format!("config line {}", index + 1);
        let (key, value) = line
            .split_once('=')
            .map(|(k, v)| (k.trim(), v.trim()))
            .ok_or_else(|| anyhow::anyhow!("{}: expected `key = value`", at()))?;
        if key != "ligand" && key != "batch" && !seen.insert(key.to_string()) {
            return Err(anyhow::anyhow!("{}: {key} is set twice", at()));
        }
        let path = || PathBuf::from(value);
        let result: Result<()> = (|| {
            match key {
                "receptor" => config.receptor = Some(path()),
                "ligand" => config.ligand.push(path()),
                "batch" => config.batch.push(path()),
                "flex" => config.flex = Some(path()),
                "out" => config.out = Some(path()),
                "dir" => config.dir = Some(path()),
                "scoring" => config.scoring = Some(value.to_string()),
                "center_x" => config.center[0] = Some(num(key, value)?),
                "center_y" => config.center[1] = Some(num(key, value)?),
                "center_z" => config.center[2] = Some(num(key, value)?),
                "size_x" => config.size[0] = Some(num(key, value)?),
                "size_y" => config.size[1] = Some(num(key, value)?),
                "size_z" => config.size[2] = Some(num(key, value)?),
                "exhaustiveness" => config.exhaustiveness = Some(num(key, value)?),
                "num_modes" => config.num_modes = Some(num(key, value)?),
                "min_rmsd" => config.min_rmsd = Some(num(key, value)?),
                "energy_range" => config.energy_range = Some(num(key, value)?),
                "spacing" => config.spacing = Some(num(key, value)?),
                "cpu" => config.cpu = Some(num(key, value)?),
                "seed" => config.seed = Some(num(key, value)?),
                "score_only" => config.score_only = Some(switch(key, value)?),
                "local_only" => config.local_only = Some(switch(key, value)?),
                "autobox" => config.autobox = Some(switch(key, value)?),
                "verbosity" => {
                    if value != "1" {
                        warn!("ignoring verbosity = {value}; set RUST_LOG instead");
                    }
                }
                _ => {
                    if let Some((_, default)) =
                        UNSUPPORTED_WITH_DEFAULT.iter().find(|(k, _)| *k == key)
                    {
                        if num::<f64>(key, value)? != *default {
                            return Err(anyhow::anyhow!(
                                "{key} is not supported by rustdock; only Vina's \
                                 default ({default}) is accepted"
                            ));
                        }
                    } else if UNSUPPORTED_SWITCHES.contains(&key) {
                        if switch(key, value)? {
                            return Err(anyhow::anyhow!("{key} is not supported by rustdock"));
                        }
                    } else if matches!(key, "maps" | "write_maps" | "unbound_energy") {
                        return Err(anyhow::anyhow!("{key} is not supported by rustdock"));
                    } else {
                        return Err(anyhow::anyhow!("unknown option {key:?}"));
                    }
                }
            }
            Ok(())
        })();
        result.with_context(at)?;
    }
    Ok(config)
}

/// Every `dock` option after merging the command line, the config file and
/// Vina 1.2.7's defaults.
#[derive(Debug, PartialEq)]
struct DockSettings {
    receptor: PathBuf,
    ligands: Vec<PathBuf>,
    batch: Vec<PathBuf>,
    flex: Option<PathBuf>,
    out: Option<PathBuf>,
    dir: Option<PathBuf>,
    scoring: String,
    /// `None` only with `score_only`, which needs no box.
    center: Option<Vector3<f64>>,
    size: Option<Vector3<f64>>,
    /// Score the ligand as given, without searching.
    score_only: bool,
    /// Minimize the ligand from where it stands, without searching.
    local_only: bool,
    /// Take the box from the ligand (Vina's `--autobox`).
    autobox: bool,
    exhaustiveness: usize,
    num_modes: usize,
    min_rmsd: f64,
    energy_range: f64,
    grid_spacing: f64,
    no_grid: bool,
    /// `None` uses every core (Vina's `cpu = 0`).
    cpu: Option<usize>,
    /// `None` draws a random seed (Vina's `seed = 0`).
    seed: Option<u64>,
}

/// One box vector (`--center` or `--size`): the comma form, else the
/// per-axis flags, each axis falling back to the config file.
fn box_vector(
    name: &str,
    joined: Option<Vec<f64>>,
    axes: [Option<f64>; 3],
    config: [Option<f64>; 3],
) -> Result<Vector3<f64>> {
    let axes = match joined {
        Some(v) if v.len() == 3 => [Some(v[0]), Some(v[1]), Some(v[2])],
        Some(v) => {
            return Err(anyhow::anyhow!(
                "--{name} needs exactly 3 comma-separated values, got {}",
                v.len()
            ))
        }
        None => axes,
    };
    let mut out = Vector3::zeros();
    let mut missing = Vec::new();
    for (i, axis) in ["x", "y", "z"].iter().enumerate() {
        match axes[i].or(config[i]) {
            Some(value) => out[i] = value,
            None => missing.push(format!("{name}_{axis}")),
        }
    }
    if !missing.is_empty() {
        return Err(anyhow::anyhow!(
            "search box incomplete, missing {}: pass --center and --size, or \
             set center_x/y/z and size_x/y/z in a --config file",
            missing.join(", ")
        ));
    }
    Ok(out)
}

/// Merge `args` with `config`; the command line wins, as in Vina.
fn resolve(args: DockArgs, config: VinaConfig) -> Result<DockSettings> {
    let receptor = args.receptor.or(config.receptor).ok_or_else(|| {
        anyhow::anyhow!("no receptor: pass --receptor FILE or set receptor in the --config file")
    })?;
    // To Boost, `ligand` is one multi-valued option, so a list on the
    // command line replaces the config file's rather than adding to it.
    let ligands = if args.ligand.is_empty() {
        config.ligand
    } else {
        args.ligand
    };
    let batch = if args.batch.is_empty() {
        config.batch
    } else {
        args.batch
    };
    let out = args.out.or(config.out);
    let dir = args.dir.or(config.dir);
    if ligands.is_empty() && batch.is_empty() {
        return Err(anyhow::anyhow!(
            "no ligand given: pass at least one --ligand FILE or --batch FILE/DIR"
        ));
    }
    if !batch.is_empty() && dir.is_none() {
        return Err(anyhow::anyhow!("--batch needs an output --dir"));
    }
    if ligands.len() > 1 && out.is_some() {
        // Vina docks several --ligand files together into one --out;
        // rustdock docks each on its own, so one --out cannot hold them.
        return Err(anyhow::anyhow!(
            "--out takes one ligand; rustdock docks each --ligand on its own, \
             so give --dir (or no --out) for {} ligands",
            ligands.len()
        ));
    }

    let score_only = args.score_only || config.score_only.unwrap_or(false);
    let local_only = args.local_only || config.local_only.unwrap_or(false);
    let autobox = args.autobox || config.autobox.unwrap_or(false);
    if score_only && local_only {
        return Err(anyhow::anyhow!(
            "--score_only and --local_only ask for different things; pick one"
        ));
    }
    if autobox && !(score_only || local_only) {
        // Vina computes an autobox only for --score_only and --local_only;
        // for a real search it would quietly use an undefined box.
        return Err(anyhow::anyhow!(
            "--autobox only applies to --score_only and --local_only; give \
             --center and --size to dock"
        ));
    }

    // A search needs a box. Scoring a pose does not, so let it go without
    // one rather than demanding six numbers that nothing reads.
    let (center, size) = match box_vector(
        "center",
        args.center,
        [args.center_x, args.center_y, args.center_z],
        config.center,
    )
    .and_then(|center| {
        let size = box_vector(
            "size",
            args.size,
            [args.size_x, args.size_y, args.size_z],
            config.size,
        )?;
        Ok((center, size))
    }) {
        Ok((center, size)) => (Some(center), Some(size)),
        // Minimizing needs a box, but --autobox will supply it.
        Err(_) if score_only || autobox => (None, None),
        Err(e) => return Err(e),
    };

    Ok(DockSettings {
        receptor,
        ligands,
        batch,
        flex: args.flex.or(config.flex),
        out,
        dir,
        scoring: args
            .scoring
            .or(config.scoring)
            .unwrap_or_else(|| "vina".to_string()),
        center,
        size,
        score_only,
        local_only,
        autobox,
        exhaustiveness: args.exhaustiveness.or(config.exhaustiveness).unwrap_or(8),
        num_modes: args.num_modes.or(config.num_modes).unwrap_or(9),
        min_rmsd: args.min_rmsd.or(config.min_rmsd).unwrap_or(1.0),
        energy_range: args.energy_range.or(config.energy_range).unwrap_or(3.0),
        grid_spacing: args.grid_spacing.or(config.spacing).unwrap_or(0.375),
        no_grid: args.no_grid,
        cpu: args.cpu.or(config.cpu).filter(|&n| n != 0),
        // Vina's seed is an `int`; a negative one is kept as its bits.
        seed: args
            .seed
            .or(config.seed.map(|s| s as u64))
            .filter(|&s| s != 0),
    })
}

/// Vina's default output for `ligand`: its path with a trailing `.pdbqt`
/// removed and `_out.pdbqt` added, so it sits next to the input.
fn default_output(ligand: &Path) -> PathBuf {
    let name = ligand.to_string_lossy();
    let stem = name.strip_suffix(".pdbqt").unwrap_or(&name);
    PathBuf::from(format!("{stem}_out.pdbqt"))
}

/// Parse a `--flex` PDBQT against the rigid receptor it belongs to.
fn load_flex(path: &Path, receptor: &Molecule) -> Result<FlexResidues> {
    info!("Loading flexible residues: {}", path.display());
    parse_flex_pdbqt(path, receptor)
        .with_context(|| format!("Failed to parse flex residue file: {}", path.display()))
}

/// Vina's `--autobox`: the box that holds `ligand` with 4 A to spare.
///
/// As `Vina::grid_dimensions_from_ligand`, the centre is the mean of every
/// atom's coordinates -- hydrogens included -- and each side is the furthest
/// atom from that centre along the axis, plus the padding, doubled and
/// rounded up to a whole Angstrom.
fn autobox_of(ligand: &Molecule) -> (Vector3<f64>, Vector3<f64>) {
    /// Vina's `buffer_size`.
    const PADDING: f64 = 4.0;

    let n = ligand.atoms.len();
    if n == 0 {
        return (Vector3::zeros(), Vector3::repeat((PADDING * 2.0).ceil()));
    }
    let center = ligand
        .atoms
        .iter()
        .map(|atom| atom.coordinates)
        .sum::<Vector3<f64>>()
        / n as f64;
    let mut reach: Vector3<f64> = Vector3::zeros();
    for atom in &ligand.atoms {
        for axis in 0..3 {
            reach[axis] = reach[axis].max((center[axis] - atom.coordinates[axis]).abs());
        }
    }
    let size = Vector3::from_iterator(reach.iter().map(|r| ((r + PADDING) * 2.0).ceil()));
    (center, size)
}

/// A pose's energies in the layout of Vina's `--score_only`: its
/// `show_score`, term for term.
///
/// The affinity is `(1) + (2) + (3) - (4)`. rustdock scores the unbound
/// system as the bound pose's own internal energy, so (2) and (4) cancel
/// and the affinity is the normalized intermolecular energy. Vina agrees on
/// a rigid receptor; with `--flex` it computes (4) from the grids while
/// reporting (2) exactly, so its affinity differs by the gap between the
/// two (see [`score_pose_with_flex`]).
fn score_lines(
    ligand: &Molecule,
    receptor: &Molecule,
    flex: Option<&FlexResidues>,
    forcefield: &dyn ForceField,
) -> Result<Vec<String>> {
    let score = score_pose_with_flex(ligand, receptor, flex, forcefield)?;
    let ligand_internal = ligand_internal_energy(ligand, forcefield)?;

    let intermolecular = score.intermolecular;
    let internal = ligand_internal + score.flex_receptor + score.flex_flex;
    // What the force field's torsional normalization took off.
    let torsional = score.affinity - intermolecular;
    let unbound = internal;

    let line = |label: &str, value: f64| format!("{label:<35}: {value:.3} (kcal/mol)");
    Ok(vec![
        format!(
            "{:<35}: {:.3} (kcal/mol) [=(1)+(2)+(3)-(4)]",
            "Estimated Free Energy of Binding", score.affinity
        ),
        line("(1) Final Intermolecular Energy", intermolecular),
        line("    Ligand - Receptor", intermolecular - score.ligand_flex),
        line("    Ligand - Flex side chains", score.ligand_flex),
        line("(2) Final Total Internal Energy", internal),
        line("    Ligand", ligand_internal),
        line("    Flex   - Receptor", score.flex_receptor),
        line("    Flex   - Flex side chains", score.flex_flex),
        line("(3) Torsional Free Energy", torsional),
        line("(4) Unbound System's Energy", unbound),
    ])
}

/// Print [`score_lines`].
fn show_score(
    ligand: &Molecule,
    receptor: &Molecule,
    flex: Option<&FlexResidues>,
    forcefield: &dyn ForceField,
) -> Result<()> {
    println!();
    for line in score_lines(ligand, receptor, flex, forcefield)? {
        println!("{line}");
    }
    Ok(())
}

/// Warn about flexible side-chain atoms outside the search box.
///
/// As in Vina, every moving atom is charged a steep penalty per Angstrom
/// outside the box, so a side chain that starts or must stay outside it
/// distorts the search.
fn warn_flex_outside_box(flex: &FlexResidues, center: &Vector3<f64>, size: &Vector3<f64>) {
    let half = size * 0.5;
    let outside = flex
        .movable_atoms()
        .filter(|a| !a.atom_type.is_hydrogen())
        .filter(|a| {
            let d = a.coordinates - center;
            (0..3).any(|k| d[k].abs() > half[k])
        })
        .count();
    if outside > 0 {
        warn!(
            "{outside} moving side-chain atoms lie outside the search box; \
             enlarge the box to cover the flexible residues"
        );
    }
}

/// Resolve the `--scoring` flag to a force field.
///
/// An unrecognized name is an error rather than a silent fall back to Vina:
/// the old behaviour handed the user Vina results while they believed they had
/// run something else, and the warning it logged was invisible by default.
fn select_forcefield(name: &str) -> Result<Box<dyn ForceField>> {
    match name.to_lowercase().as_str() {
        "vina" => Ok(Box::new(VinaForceField::default())),
        "vinardo" => Ok(Box::new(VinardoForceField::default())),
        "ad4" => {
            warn!(
                "The AD4 scoring function is experimental and is not validated \
                 against AutoDock4; treat its numbers as indicative only."
            );
            Ok(Box::new(AD4ForceField::default()))
        }
        other => Err(anyhow::anyhow!(
            "unknown scoring function {other:?} (supported: vina, vinardo, ad4)"
        )),
    }
}

/// What every ligand of a `dock` run shares.
struct DockJob<'a> {
    optimizer: &'a MonteCarlo,
    target: &'a Target,
    forcefield: &'a dyn ForceField,
    num_modes: usize,
    energy_range: f64,
}

impl DockJob<'_> {
    /// Dock the ligand in `path`: its modes within `energy_range` of the
    /// best, at most `num_modes`, with affinities torsion-normalized.
    fn dock(&self, path: &Path) -> Result<Vec<DockingResult>> {
        info!("Loading ligand: {}", path.display());
        let ligand = parse_pdbqt(path)
            .with_context(|| format!("Failed to parse ligand file: {}", path.display()))?;
        let results = self.optimizer.dock_into(
            self.target,
            &ligand,
            self.forcefield,
            self.num_modes,
            0, // Vina's step count for this ligand
        )?;

        // Apply the force field's torsional normalization (Vina:
        // affinity = energy / (1 + weight_rot * num_tors))
        let num_tors = ligand.vina_torsion_count();
        let normalized: Vec<_> = results
            .into_iter()
            .map(|mut r| {
                r.energy = self
                    .forcefield
                    .apply_torsional_normalization(r.energy, num_tors);
                r
            })
            .collect();

        let min_energy = normalized.first().map(|r| r.energy).unwrap_or(0.0);
        Ok(normalized
            .into_iter()
            .filter(|r| r.energy <= min_energy + self.energy_range)
            .take(self.num_modes)
            .collect())
    }
}

/// The ligand files named by `--batch`: files as given, directories as
/// their `.pdbqt` files in name order. Two ligands with the same file stem
/// would write the same output, so that is an error.
fn expand_batch(paths: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut ligands = Vec::new();
    for path in paths {
        if path.is_dir() {
            let mut found: Vec<PathBuf> = std::fs::read_dir(path)
                .with_context(|| format!("Failed to read directory: {}", path.display()))?
                .filter_map(|entry| entry.ok().map(|e| e.path()))
                .filter(|p| {
                    p.is_file()
                        && p.extension()
                            .is_some_and(|e| e.eq_ignore_ascii_case("pdbqt"))
                })
                .collect();
            found.sort();
            ligands.extend(found);
        } else {
            ligands.push(path.clone());
        }
    }
    let mut stems = HashSet::new();
    for ligand in &ligands {
        let stem = ligand.file_stem().unwrap_or_default().to_os_string();
        if !stems.insert(stem) {
            return Err(anyhow::anyhow!(
                "two --batch ligands share the name {:?}; their outputs would collide",
                ligand.file_stem().unwrap_or_default()
            ));
        }
    }
    Ok(ligands)
}

/// Dock `ligands` in parallel into `job`'s target, writing each one's poses
/// to `dir/<stem>_out.pdbqt` and every best affinity to `dir/scores.tsv`, in
/// input order. A ligand that fails is logged and listed with its error.
fn run_batch(job: &DockJob, ligands: &[PathBuf], dir: &Path) -> Result<()> {
    std::fs::create_dir_all(dir)
        .with_context(|| format!("Failed to create output directory: {}", dir.display()))?;
    let start = Instant::now();
    let done = AtomicUsize::new(0);
    let total = ligands.len();

    let rows: Vec<String> = ligands
        .par_iter()
        .map(|path| {
            let stem = path.file_stem().unwrap_or_default().to_string_lossy();
            let out = dir.join(format!("{stem}_out.pdbqt"));
            let t0 = Instant::now();
            let outcome = job.dock(path).and_then(|results| {
                write_docking_results(&results, &out)
                    .with_context(|| format!("Failed to write {}", out.display()))?;
                Ok(results)
            });
            let seconds = t0.elapsed().as_secs_f64();
            let n = done.fetch_add(1, Ordering::Relaxed) + 1;
            match outcome {
                Ok(results) => {
                    let best = results.first().map(|r| r.energy);
                    info!(
                        "[{n}/{total}] {stem}: {} kcal/mol, {seconds:.1} s",
                        best.map_or("no pose".to_string(), |e| format!("{e:.3}"))
                    );
                    format!(
                        "{stem}\t{}\t{}\t{seconds:.2}\t{}\t",
                        best.map_or(String::new(), |e| format!("{e:.3}")),
                        results.len(),
                        out.display()
                    )
                }
                Err(e) => {
                    warn!("[{n}/{total}] {stem}: failed: {e:#}");
                    let reason = format!("{e:#}").replace(['\t', '\n'], " ");
                    format!("{stem}\t\t0\t{seconds:.2}\t\t{reason}")
                }
            }
        })
        .collect();

    let summary = dir.join("scores.tsv");
    let mut text = String::from("ligand\taffinity\tmodes\tseconds\toutput\terror\n");
    for row in &rows {
        text.push_str(row);
        text.push('\n');
    }
    std::fs::write(&summary, text)
        .with_context(|| format!("Failed to write {}", summary.display()))?;

    let failed = rows.iter().filter(|r| !r.ends_with('\t')).count();
    let seconds = start.elapsed().as_secs_f64();
    info!(
        "Docked {} of {total} ligands in {seconds:.1} s ({:.1} ligands/min); scores in {}",
        total - failed,
        total as f64 * 60.0 / seconds.max(1e-9),
        summary.display()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustdock_vina::forcefield::vina::VinaForceField;

    /// Parse `dock` arguments as the binary would.
    fn dock_args(argv: &[&str]) -> DockArgs {
        let cli = Cli::try_parse_from(["rustdock", "dock"].iter().chain(argv)).unwrap();
        match cli.command {
            Commands::Dock(args) => args,
            other => panic!("not dock: {other:?}"),
        }
    }

    const VINA_CONFIG: &str = "\
# A config file as written for AutoDock Vina 1.2.7
receptor = rec.pdbqt
ligand = lig.pdbqt   # trailing comment

center_x = 1.5
center_y = -2
center_z = 3
size_x = 20
size_y = 22
size_z = 24
exhaustiveness = 16
num_modes = 20
energy_range = 4
min_rmsd = 2
seed = 7
cpu = 0
spacing = 0.5
scoring = vinardo
out = poses.pdbqt
verbosity = 1
max_evals = 0
weight_gauss1 = -0.035579
score_only = false
";

    #[test]
    fn a_vina_config_runs_unchanged() {
        let settings = resolve(dock_args(&[]), parse_config(VINA_CONFIG).unwrap()).unwrap();
        assert_eq!(
            settings,
            DockSettings {
                receptor: "rec.pdbqt".into(),
                ligands: vec!["lig.pdbqt".into()],
                batch: vec![],
                flex: None,
                out: Some("poses.pdbqt".into()),
                dir: None,
                scoring: "vinardo".into(),
                center: Some(Vector3::new(1.5, -2.0, 3.0)),
                size: Some(Vector3::new(20.0, 22.0, 24.0)),
                score_only: false,
                local_only: false,
                autobox: false,
                exhaustiveness: 16,
                num_modes: 20,
                min_rmsd: 2.0,
                energy_range: 4.0,
                grid_spacing: 0.5,
                no_grid: false,
                cpu: None,
                seed: Some(7),
            }
        );
    }

    #[test]
    fn the_command_line_wins_over_the_config() {
        let args = dock_args(&[
            "--ligand",
            "a.pdbqt",
            "--center_x",
            "-9",
            "--exhaustiveness",
            "32",
            "--num_modes",
            "3",
            "--spacing",
            "0.4",
            "--seed",
            "0",
            "--out",
            "x.pdbqt",
        ]);
        let s = resolve(args, parse_config(VINA_CONFIG).unwrap()).unwrap();
        assert_eq!(s.ligands, vec![PathBuf::from("a.pdbqt")]);
        assert_eq!(s.center, Some(Vector3::new(-9.0, -2.0, 3.0)));
        assert_eq!((s.exhaustiveness, s.num_modes), (32, 3));
        assert_eq!(s.grid_spacing, 0.4);
        assert_eq!(s.seed, None, "seed 0 draws a random seed, as in Vina");
        assert_eq!(s.out, Some(PathBuf::from("x.pdbqt")));
        assert_eq!(s.receptor, PathBuf::from("rec.pdbqt"));
    }

    #[test]
    fn defaults_are_vinas() {
        let args = dock_args(&[
            "--receptor",
            "r.pdbqt",
            "--ligand",
            "l.pdbqt",
            "--center",
            "0,0,0",
            "--size",
            "20,20,20",
        ]);
        let s = resolve(args, VinaConfig::default()).unwrap();
        assert_eq!(s.scoring, "vina");
        assert_eq!((s.exhaustiveness, s.num_modes), (8, 9));
        assert_eq!(
            (s.min_rmsd, s.energy_range, s.grid_spacing),
            (1.0, 3.0, 0.375)
        );
        assert_eq!((s.cpu, s.seed), (None, None));
    }

    #[test]
    fn config_errors_name_the_line() {
        let cases = [
            ("bogus = 1", "unknown option"),
            ("seed = 1\nseed = 2", "set twice"),
            ("center_x", "expected `key = value`"),
            ("size_x = big", "not a valid number"),
            ("randomize_only = true", "not supported"),
            ("weight_gauss1 = -0.1", "only Vina's default"),
            ("maps = rec", "not supported"),
            ("config = other.txt", "unknown option"),
        ];
        for (text, expected) in cases {
            let err = format!("{:#}", parse_config(text).unwrap_err());
            assert!(err.contains("config line"), "{text}: {err}");
            assert!(err.contains(expected), "{text}: {err}");
        }
    }

    #[test]
    fn ligand_lines_accumulate() {
        let config = parse_config("ligand = a.pdbqt\nligand = b.pdbqt").unwrap();
        assert_eq!(
            config.ligand,
            vec![PathBuf::from("a.pdbqt"), "b.pdbqt".into()]
        );
    }

    #[test]
    fn incomplete_box_names_the_missing_axes() {
        let args = dock_args(&["--receptor", "r", "--ligand", "l", "--center", "0,0,0"]);
        let err = resolve(args, parse_config("size_x = 20").unwrap()).unwrap_err();
        assert!(err.to_string().contains("size_y, size_z"), "{err}");
    }

    #[test]
    fn one_out_cannot_hold_several_ligands() {
        let args = dock_args(&[
            "--receptor",
            "r",
            "--ligand",
            "a",
            "b",
            "--center",
            "0,0,0",
            "--size",
            "9,9,9",
            "--out",
            "o.pdbqt",
        ]);
        assert!(resolve(args, VinaConfig::default()).is_err());
    }

    #[test]
    fn center_and_center_x_conflict() {
        let argv = ["rustdock", "dock", "--center", "0,0,0", "--center_x", "1"];
        assert!(Cli::try_parse_from(argv).is_err());
    }

    #[test]
    fn score_only_needs_no_search_box() {
        let args = dock_args(&["--receptor", "r", "--ligand", "l", "--score_only"]);
        let s = resolve(args, VinaConfig::default()).unwrap();
        assert!(s.score_only);
        assert_eq!((s.center, s.size), (None, None));

        // Docking still insists on one.
        let args = dock_args(&["--receptor", "r", "--ligand", "l"]);
        assert!(resolve(args, VinaConfig::default()).is_err());
    }

    #[test]
    fn autobox_only_applies_to_score_only() {
        let args = dock_args(&["--receptor", "r", "--ligand", "l", "--autobox"]);
        let err = resolve(args, VinaConfig::default()).unwrap_err();
        assert!(err.to_string().contains("--score_only"), "{err}");

        let args = dock_args(&[
            "--receptor",
            "r",
            "--ligand",
            "l",
            "--autobox",
            "--score_only",
        ]);
        let s = resolve(args, VinaConfig::default()).unwrap();
        assert!(s.autobox && s.score_only);
    }

    #[test]
    fn a_config_may_ask_for_score_only() {
        let config = parse_config("score_only = yes\nautobox = 1").unwrap();
        let s = resolve(dock_args(&["--receptor", "r", "--ligand", "l"]), config).unwrap();
        assert!(s.score_only && s.autobox);
        // Vina still rejects what rustdock cannot do.
        assert!(parse_config("randomize_only = true").is_err());
    }

    #[test]
    fn local_only_needs_a_box_or_autobox() {
        // Minimizing scores on grids, so it needs a box -- unlike scoring.
        let args = dock_args(&["--receptor", "r", "--ligand", "l", "--local_only"]);
        assert!(resolve(args, VinaConfig::default()).is_err());

        let args = dock_args(&[
            "--receptor",
            "r",
            "--ligand",
            "l",
            "--local_only",
            "--autobox",
        ]);
        let s = resolve(args, VinaConfig::default()).unwrap();
        assert!(s.local_only && s.autobox);
        assert_eq!(s.center, None, "--autobox takes the box from the ligand");

        let args = dock_args(&[
            "--receptor",
            "r",
            "--ligand",
            "l",
            "--local_only",
            "--center",
            "0,0,0",
            "--size",
            "20,20,20",
        ]);
        let s = resolve(args, VinaConfig::default()).unwrap();
        assert!(s.local_only && s.center.is_some());
    }

    #[test]
    fn score_only_and_local_only_are_exclusive() {
        let args = dock_args(&[
            "--receptor",
            "r",
            "--ligand",
            "l",
            "--score_only",
            "--local_only",
        ]);
        let err = resolve(args, VinaConfig::default()).unwrap_err();
        assert!(err.to_string().contains("pick one"), "{err}");
    }

    #[test]
    fn a_config_may_ask_for_local_only() {
        let config = parse_config("local_only = 1\nautobox = true").unwrap();
        let s = resolve(dock_args(&["--receptor", "r", "--ligand", "l"]), config).unwrap();
        assert!(s.local_only && s.autobox);
    }

    /// Vina's `grid_dimensions_from_ligand`: the mean of every atom, and
    /// each side the furthest reach plus 4 A, doubled and rounded up.
    #[test]
    fn autobox_holds_the_ligand_with_padding() {
        let mut molecule = Molecule::new("lig");
        for (x, y, z) in [(0.0, 0.0, 0.0), (4.0, 1.0, 0.0), (2.0, -1.0, 0.5)] {
            molecule.add_atom(rustdock_vina::Atom::new(
                rustdock_vina::AtomType::Carbon,
                Vector3::new(x, y, z),
                "C".to_string(),
                1,
                "LIG".to_string(),
                1,
                'A',
                0.0,
            ));
        }
        let (center, size) = autobox_of(&molecule);
        assert_approx_eq::assert_approx_eq!(center.x, 2.0, 1e-12);
        assert_approx_eq::assert_approx_eq!(center.y, 0.0, 1e-12);
        // Reach 2.0 in x, 1.0 in y, 0.333 in z: (reach + 4) * 2, rounded up.
        assert_eq!(size, Vector3::new(12.0, 10.0, 9.0));

        // Every atom must sit inside the box it produced.
        for atom in &molecule.atoms {
            for axis in 0..3 {
                assert!((atom.coordinates[axis] - center[axis]).abs() <= size[axis] / 2.0);
            }
        }
    }

    /// The breakdown is Vina's `show_score`, line for line, and its own
    /// arithmetic holds: the affinity is (1) + (2) + (3) - (4).
    #[test]
    fn the_breakdown_is_vinas() {
        let data = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("test_data");
        let receptor = parse_pdbqt(data.join("receptor.pdbqt")).unwrap();
        let ligand = parse_pdbqt(data.join("ligand.pdbqt")).unwrap();
        let forcefield = VinaForceField::default();
        let lines = score_lines(&ligand, &receptor, None, &forcefield).unwrap();

        let labels: Vec<&str> = lines
            .iter()
            .map(|l| l.split(':').next().unwrap().trim_end())
            .collect();
        assert_eq!(
            labels,
            [
                "Estimated Free Energy of Binding",
                "(1) Final Intermolecular Energy",
                "    Ligand - Receptor",
                "    Ligand - Flex side chains",
                "(2) Final Total Internal Energy",
                "    Ligand",
                "    Flex   - Receptor",
                "    Flex   - Flex side chains",
                "(3) Torsional Free Energy",
                "(4) Unbound System's Energy",
            ]
        );

        let value = |i: usize| -> f64 {
            lines[i]
                .split(':')
                .nth(1)
                .unwrap()
                .split_whitespace()
                .next()
                .unwrap()
                .parse()
                .unwrap()
        };
        // [=(1)+(2)+(3)-(4)], as the first line claims.
        assert_approx_eq::assert_approx_eq!(
            value(0),
            value(1) + value(4) + value(8) - value(9),
            2e-3
        );
        // Without flexible side chains their terms are zero, and the
        // internal energy is the ligand's alone.
        assert_eq!((value(3), value(6), value(7)), (0.0, 0.0, 0.0));
        assert_approx_eq::assert_approx_eq!(value(4), value(5), 1e-12);
        assert_approx_eq::assert_approx_eq!(value(2), value(1), 1e-12);
    }

    #[test]
    fn default_output_is_beside_the_ligand() {
        assert_eq!(
            default_output(Path::new("dir/lig.pdbqt")),
            PathBuf::from("dir/lig_out.pdbqt")
        );
        assert_eq!(
            default_output(Path::new("lig.mol")),
            PathBuf::from("lig.mol_out.pdbqt")
        );
    }
}
