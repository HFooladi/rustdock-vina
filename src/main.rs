//! Main executable for rustdock-vina

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use log::{error, info, warn};
use nalgebra::Vector3;
use rayon::prelude::*;
use rustdock_vina::forcefield::{
    ad4::AD4ForceField, vina::VinaForceField, vinardo::VinardoForceField, ForceField,
};
use rustdock_vina::io::{parse_pdbqt, write_docking_results};
use rustdock_vina::optimization::monte_carlo::{MonteCarlo, MonteCarloParams, Target};
use rustdock_vina::optimization::DockingResult;
use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

/// Command-line arguments for the application
#[derive(Parser, Debug)]
#[clap(
    name = "rustdock-vina",
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
    Dock {
        /// PDBQT file containing the receptor
        #[clap(long, value_parser)]
        receptor: PathBuf,

        /// PDBQT file containing the ligand to dock
        #[clap(long, value_parser)]
        ligand: Vec<PathBuf>,

        /// PDBQT file containing flexible residues in the receptor
        #[clap(long, value_parser)]
        flex: Option<PathBuf>,

        /// Configuration file for docking
        #[clap(long, short, value_parser)]
        config: Option<PathBuf>,

        /// Center of the search box (x,y,z)
        #[clap(long, value_parser, value_delimiter = ',')]
        center: Option<Vec<f64>>,

        /// Size of the search box (x,y,z)
        #[clap(long, value_parser, value_delimiter = ',')]
        size: Option<Vec<f64>>,

        /// Scoring function to use (vina, ad4, vinardo)
        #[clap(long, default_value = "vina")]
        scoring: String,

        /// Exhaustiveness of the search (higher values increase accuracy but take longer)
        #[clap(long, default_value_t = 8)]
        exhaustiveness: usize,

        /// Number of binding modes to generate
        #[clap(long, default_value_t = 9)]
        num_modes: usize,

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

        /// Number of threads (default: all cores)
        #[clap(long)]
        cpu: Option<usize>,

        /// Energy range for output poses (kcal/mol)
        #[clap(long, default_value_t = 3.0)]
        energy_range: f64,

        /// Spacing of the precomputed affinity grids (Angstroms)
        #[clap(long, default_value_t = 0.375)]
        grid_spacing: f64,

        /// Search on exact pairwise energies instead of precomputed grids (slow)
        #[clap(long)]
        no_grid: bool,

        /// Random seed; the same seed and inputs reproduce the same poses
        #[clap(long)]
        seed: Option<u64>,
    },

    /// Score an existing pose
    Score {
        /// PDBQT file containing the receptor
        #[clap(long, value_parser)]
        receptor: PathBuf,

        /// PDBQT file containing the ligand pose to score
        #[clap(long, value_parser)]
        ligand: PathBuf,

        /// PDBQT file containing flexible residues in the receptor
        #[clap(long, value_parser)]
        flex: Option<PathBuf>,

        /// Scoring function to use (vina, ad4, vinardo)
        #[clap(long, default_value = "vina")]
        scoring: String,
    },
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
        Commands::Dock {
            receptor,
            ligand,
            flex,
            config,
            center,
            size,
            scoring,
            exhaustiveness,
            num_modes,
            out,
            batch,
            dir,
            cpu,
            energy_range,
            grid_spacing,
            no_grid,
            seed,
        } => {
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
            if ligand.is_empty() && batch.is_empty() {
                return Err(anyhow::anyhow!(
                    "no ligand given: pass at least one --ligand FILE or --batch FILE/DIR"
                ));
            }
            if !batch.is_empty() && dir.is_none() {
                return Err(anyhow::anyhow!("--batch needs an output --dir"));
            }

            // Parse configuration
            let (center, box_size) = match (center, size, config) {
                (Some(c), Some(s), _) if c.len() == 3 && s.len() == 3 => {
                    let center = Vector3::new(c[0], c[1], c[2]);
                    let box_size = Vector3::new(s[0], s[1], s[2]);
                    (center, box_size)
                }
                (Some(c), _, _) if c.len() != 3 => {
                    return Err(anyhow::anyhow!(
                        "--center needs exactly 3 comma-separated values, got {}",
                        c.len()
                    ));
                }
                (_, Some(sz), _) if sz.len() != 3 => {
                    return Err(anyhow::anyhow!(
                        "--size needs exactly 3 comma-separated values, got {}",
                        sz.len()
                    ));
                }
                (_, _, Some(config_path)) => {
                    // Parse configuration file
                    let config_str = std::fs::read_to_string(&config_path).with_context(|| {
                        format!("Failed to read config file: {}", config_path.display())
                    })?;

                    parse_config(&config_str)?
                }
                _ => {
                    error!("Either center and size or a config file must be provided");
                    return Err(anyhow::anyhow!(
                        "Either center and size or a config file must be provided"
                    ));
                }
            };

            info!(
                "Search box center: ({}, {}, {})",
                center.x, center.y, center.z
            );
            info!(
                "Search box size: ({}, {}, {})",
                box_size.x, box_size.y, box_size.z
            );

            // Parse receptor
            info!("Loading receptor: {}", receptor.display());
            let receptor_molecule = parse_pdbqt(&receptor).with_context(|| {
                format!("Failed to parse receptor file: {}", receptor.display())
            })?;

            // Parse flexible residues if provided
            let _flex_molecule = if let Some(flex_path) = flex {
                info!("Loading flexible residues: {}", flex_path.display());
                Some(parse_pdbqt(&flex_path).with_context(|| {
                    format!("Failed to parse flex file: {}", flex_path.display())
                })?)
            } else {
                None
            };

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
                seed,
                ..MonteCarloParams::default()
            });
            // Every ligand shares one set of receptor grids.
            let target =
                optimizer.target(receptor_molecule, forcefield.as_ref(), center, box_size)?;
            let job = DockJob {
                optimizer: &optimizer,
                target: &target,
                forcefield: forcefield.as_ref(),
                num_modes,
                energy_range,
            };

            // Process each ligand
            for ligand_path in ligand {
                // Determine output path
                let output_path = if let Some(out_path) = &out {
                    out_path.clone()
                } else if let Some(dir_path) = &dir {
                    let ligand_name = ligand_path.file_stem().unwrap().to_string_lossy();
                    dir_path.join(format!("{}_out.pdbqt", ligand_name))
                } else {
                    let ligand_stem = ligand_path.file_stem().unwrap().to_string_lossy();
                    PathBuf::from(format!("{}_out.pdbqt", ligand_stem))
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
                    let rmsd = r.rmsd.unwrap_or(0.0);
                    println!(
                        "{:4} | {:10.3} | {:8.3} | {:8.3}",
                        i + 1,
                        r.energy,
                        rmsd,
                        rmsd
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

            // Parse flexible residues if provided
            let _flex_molecule = if let Some(flex_path) = flex {
                info!("Loading flexible residues: {}", flex_path.display());
                Some(parse_pdbqt(&flex_path).with_context(|| {
                    format!("Failed to parse flex file: {}", flex_path.display())
                })?)
            } else {
                None
            };

            // Parse ligand
            info!("Loading ligand: {}", ligand.display());
            let ligand_molecule = parse_pdbqt(&ligand)
                .with_context(|| format!("Failed to parse ligand file: {}", ligand.display()))?;

            // Calculate interaction energy between ligand and receptor
            let cutoff = 8.0;
            let mut intermolecular_energy = 0.0;
            let mut num_pairs = 0;
            let mut positive_energy = 0.0;
            let mut negative_energy = 0.0;

            println!("\nLigand atoms: {}", ligand_molecule.atoms.len());
            println!("Ligand bonds: {}", ligand_molecule.bonds.len());
            println!("Ligand torsions: {}", ligand_molecule.torsions.len());
            println!("Receptor atoms: {}", receptor_molecule.atoms.len());

            for lig_atom in &ligand_molecule.atoms {
                for rec_atom in &receptor_molecule.atoms {
                    let distance = lig_atom.distance(rec_atom);
                    if distance < cutoff && distance > 0.01 {
                        if let Ok(pair_energy) =
                            forcefield.atom_pair_energy(lig_atom, rec_atom, distance)
                        {
                            intermolecular_energy += pair_energy;
                            num_pairs += 1;
                            if pair_energy > 0.0 {
                                positive_energy += pair_energy;
                            } else {
                                negative_energy += pair_energy;
                            }
                        }
                    }
                }
            }

            // Vina formula: affinity = intermolecular / (1 + weight_rot * num_tors)
            // The internal energy cancels with unbound energy in the full formula
            let num_tors = ligand_molecule.vina_torsion_count();
            let affinity =
                forcefield.apply_torsional_normalization(intermolecular_energy, num_tors);
            let torsional_effect = intermolecular_energy - affinity;

            println!("\nScoring breakdown:");
            println!("  Atom pairs within {}A: {}", cutoff, num_pairs);
            println!("  Positive energy (repulsion): {:.3}", positive_energy);
            println!("  Negative energy (attraction): {:.3}", negative_energy);
            println!("  Intermolecular energy: {:.3}", intermolecular_energy);
            println!(
                "  Torsional effect: {:.3} ({} torsions, weighted count {})",
                torsional_effect,
                ligand_molecule.torsions.len(),
                num_tors
            );
            println!("\nAffinity: {:.3} (kcal/mol)", affinity);
            info!("Scoring completed successfully");
        }
    }

    Ok(())
}

/// Parse center and box size from a configuration file
fn parse_config(config_str: &str) -> Result<(Vector3<f64>, Vector3<f64>)> {
    let mut center_x = None;
    let mut center_y = None;
    let mut center_z = None;
    let mut size_x = None;
    let mut size_y = None;
    let mut size_z = None;

    for line in config_str.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.splitn(2, '=').collect();
        if parts.len() != 2 {
            continue;
        }

        let key = parts[0].trim();
        let value = parts[1].trim();

        match key {
            "center_x" => center_x = Some(value.parse::<f64>()?),
            "center_y" => center_y = Some(value.parse::<f64>()?),
            "center_z" => center_z = Some(value.parse::<f64>()?),
            "size_x" => size_x = Some(value.parse::<f64>()?),
            "size_y" => size_y = Some(value.parse::<f64>()?),
            "size_z" => size_z = Some(value.parse::<f64>()?),
            // A real Vina config also carries receptor/ligand/out/exhaustiveness/
            // num_modes/cpu/seed. None of those are honoured here, so say so
            // rather than silently running with defaults.
            other => warn!(
                "ignoring unsupported config key {other:?}: only center_x/y/z \
                 and size_x/y/z are recognised"
            ),
        }
    }

    // Check if all values are present
    match (center_x, center_y, center_z, size_x, size_y, size_z) {
        (Some(cx), Some(cy), Some(cz), Some(sx), Some(sy), Some(sz)) => {
            Ok((Vector3::new(cx, cy, cz), Vector3::new(sx, sy, sz)))
        }
        _ => Err(anyhow::anyhow!(
            "Missing configuration values for center or box size"
        )),
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
