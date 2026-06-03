use std::{env, path::PathBuf, process};

use pharmcat::cli::{CliAction, PHARMCAT_VERSION_PREFIX, help_text, parse_pharmcat_args};
use pharmcat::pipeline::{CliPipelineOptions, PharmcatResourcePaths, run_cli_config};

fn main() {
    match parse_pharmcat_args(env::args().skip(1)) {
        Ok(CliAction::Help) => {
            println!("{}", help_text());
        }
        Ok(CliAction::Version) => {
            println!("{PHARMCAT_VERSION_PREFIX} {}", env!("CARGO_PKG_VERSION"));
        }
        Ok(CliAction::Run(config)) => {
            let resource_root = env::var_os("PHARMCAT_RESOURCE_DIR")
                .map(PathBuf::from)
                .unwrap_or_else(|| {
                    "repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat".into()
                });
            let options = CliPipelineOptions {
                resources: PharmcatResourcePaths::from_resource_root(resource_root),
                mode: pharmcat::pipeline::PipelineMode::Cli,
            };
            if let Err(err) = run_cli_config(&config, &options) {
                eprintln!("{err}");
                process::exit(2);
            }
        }
        Err(err) => {
            eprintln!("{err}");
            process::exit(1);
        }
    }
}
