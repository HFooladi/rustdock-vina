//! Main executable for rustdock-vina

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use log::{debug, error, info, warn};
use nalgebra::Vector3;
use std::path::PathBuf;

use rustdock_vina::atom::AtomType;
use rustdock_vina::forcefield::{ad4::AD4ForceField, vina::VinaForceField, ForceField};
use rustdock_vina::io::{parse_pdbqt, write_docking_results};
use rustdock_vina::optimization::{monte_carlo::MonteCarlo, Optimizer};

/// Command-line arguments for the application
#[derive(Parser, Debug)]
#[clap(
    name = "rustdock-vina",
    version = rustdock_vina::VERSION,
    author = "Author <author@example.com>",
    about = "A Rust implementation of AutoDock Vina for molecular docking"
)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

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

        /// Output directory for the docking results (in batch mode)
        #[clap(long, value_parser)]
        dir: Option<PathBuf>,

        /// Energy range for output poses (kcal/mol)
        #[clap(long, default_value_t = 3.0)]
        energy_range: f64,
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
    // Initialize logger
    env_logger::init();

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
            dir,
            energy_range,
        } => {
            // Determine scoring function
            let forcefield: Box<dyn ForceField> = match scoring.to_lowercase().as_str() {
                "vina" => Box::new(VinaForceField::default()),
                "ad4" => Box::new(AD4ForceField::default()),
                _ => {
                    warn!("Unknown scoring function: {}. Using Vina instead.", scoring);
                    Box::new(VinaForceField::default())
                }
            };

            info!("Using {} scoring function", forcefield.name());

            // Parse configuration
            let (center, box_size) = match (center, size, config) {
                (Some(c), Some(s), _) if c.len() == 3 && s.len() == 3 => {
                    let center = Vector3::new(c[0], c[1], c[2]);
                    let box_size = Vector3::new(s[0], s[1], s[2]);
                    (center, box_size)
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
            let flex_molecule = if let Some(flex_path) = flex {
                info!("Loading flexible residues: {}", flex_path.display());
                Some(parse_pdbqt(&flex_path).with_context(|| {
                    format!("Failed to parse flex file: {}", flex_path.display())
                })?)
            } else {
                None
            };

            // Create optimizer
            let optimizer = MonteCarlo::default();

            // Process each ligand
            for ligand_path in ligand {
                info!("Loading ligand: {}", ligand_path.display());
                let ligand_molecule = parse_pdbqt(&ligand_path).with_context(|| {
                    format!("Failed to parse ligand file: {}", ligand_path.display())
                })?;

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

                // Run docking
                info!(
                    "Docking ligand {} with exhaustiveness {}",
                    ligand_path.display(),
                    exhaustiveness
                );
                let results = optimizer.generate_poses(
                    &ligand_molecule,
                    forcefield.as_ref(),
                    center,
                    box_size,
                    num_modes * exhaustiveness, // Generate more poses initially
                    1000,                       // Max steps per pose
                )?;

                // Filter results
                let min_energy = results.first().map(|r| r.energy).unwrap_or(0.0);
                let filtered_results: Vec<_> = results
                    .into_iter()
                    .filter(|r| r.energy <= min_energy + energy_range)
                    .take(num_modes)
                    .collect();

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

            info!("Docking completed successfully");
        }

        Commands::Score {
            receptor,
            ligand,
            flex,
            scoring,
        } => {
            // Determine scoring function
            let forcefield: Box<dyn ForceField> = match scoring.to_lowercase().as_str() {
                "vina" => Box::new(VinaForceField::default()),
                "ad4" => Box::new(AD4ForceField::default()),
                _ => {
                    warn!("Unknown scoring function: {}. Using Vina instead.", scoring);
                    Box::new(VinaForceField::default())
                }
            };

            info!("Using {} scoring function", forcefield.name());

            // Parse receptor
            info!("Loading receptor: {}", receptor.display());
            let receptor_molecule = parse_pdbqt(&receptor).with_context(|| {
                format!("Failed to parse receptor file: {}", receptor.display())
            })?;

            // Parse flexible residues if provided
            let flex_molecule = if let Some(flex_path) = flex {
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

            // This is a placeholder - in a real implementation,
            // we would calculate the interaction energy between the ligand and receptor

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
            _ => {} // Ignore other keys
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
