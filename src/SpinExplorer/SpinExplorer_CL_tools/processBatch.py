import nmrglue as ng # type: ignore
import os
from SpinExplorer.SpinExplorer_CL_tools.bruker_params_cl import ParameterExtractorBruker
from SpinExplorer.SpinExplorer_CL_tools.convert_nmrglue_cl import Convert_nmrglue
from SpinExplorer.SpinExplorer_CL_tools.make_parameter_file_cl import parameter_write_cl
from SpinExplorer.SpinExplorer_CL_tools.pulse_sequence_parsing import PulseSequenceParser
from SpinExplorer.SpinExplorer_CL_tools.config_register import registry
from SpinExplorer.SpinExplorer_CL_tools.processSpec import FindingParameters
from pathlib import Path
import pandas as pd # type: ignore
import argparse
import yaml # type: ignore


def filter_experiments(df: pd.DataFrame, experiment_prefix: str, protein: str) -> pd.DataFrame:
    """
    Filter rows where:
    - Experiment column starts with the given prefix
    - Name column contains the given protein string
    
    Args:
        df: Input dataframe
        experiment_prefix: Prefix to match in Experiment column (e.g. 'zgesgp', 't1rho')
        protein: Protein name to match in Name column (e.g. 'aSyn')
    
    Returns:
        Filtered dataframe
    """
    exp_mask = df["Experiment"].str.startswith(experiment_prefix)
    
    return df[exp_mask]


def group_by_base_title(df: pd.DataFrame, protein: str) -> list[pd.DataFrame]:
    base_titles = df[~df["Title"].str.contains(" ")]["Title"].unique()
    
    groups = []
    titles = []
    for base_title in base_titles:
        base_mask = df["Title"] == base_title
        compound_mask = (
            df["Title"].str.contains(base_title, regex=False) &
            df["Title"].str.contains(protein, case=False)
        )
        group = df[base_mask | compound_mask]
        if len(group) > 1:
            groups.append(group)
        titles.append(base_title)
    
    return groups, titles

def process_from_config(config_path: str, organise_by: str | None = None) -> None:
    with open(config_path, "r") as f:
        config_data = yaml.safe_load(f)

    base_dir = Path(config_data.get("output_dir", Path.cwd()))
    df = pd.read_csv(config_data["csv_path"])

    for exp in config_data["experiments"]:
        if('_icon' in exp):
            exp_name = exp.split('_icon')[0]
            suffix = exp.split('.')[-1]
            register_exp = exp_name + '.' + suffix
        config = registry._registry[register_exp]
        for prot in config_data["proteins"]:
            filtered_exps = filter_experiments(df, exp, prot)
            groups, titles = group_by_base_title(filtered_exps, prot)
            for i, group in enumerate(groups):
                output_dir = _resolve_output_dir(base_dir, organise_by, prot, exp, titles[i])
                write_1d_multi_session(group, exp, prot, config, output_dir=output_dir)


def write_1d_multi_session(df, exp, protein_name, config, outy_name=None, outy_folder=None, output_dir=None):
    # output_dir (from --organise-by) takes precedence over legacy outy_folder
    if output_dir is not None:
        resolved_folder = output_dir
    elif outy_folder is not None:
        resolved_folder = outy_folder
    else:
        resolved_folder = Path('./')

    if outy_name is None:
        titles = df["Title"]
        base = [t for t in titles if " " not in t][0]
        outy_name = base + '_' + protein_name + '_' + exp + '.session'

    with open(resolved_folder / Path(outy_name), 'w') as outy:
        outy.write('1D\n')
        outy.write('MultiplotMode:True\n')

        for i, (_, row) in enumerate(df.iterrows()):
            outy.write(f'file_path:{str(Path.cwd())+'/'+str(int(row['Expno']))+'/test.ft'}\n')
            outy.write(f'title:{str(row['Expno'])}\n')
            outy.write(f'p0_coarse:0.0\n')
            outy.write(f'p1_coarse:0.0\n')
            outy.write(f'p0_fine:0.0\n')
            outy.write(f'p1_fine:0.0\n')
            outy.write(f'colour:{i}\n')
            outy.write(f'linewidth:0.5\n')
            outy.write(f'reference_range:0\n')
            outy.write(f'reference_value:0.0\n')
            outy.write(f'vertical_range:0\n')
            outy.write(f'vertical_value:0.0\n')
            outy.write(f'multiply_range:0\n')
            outy.write(f'multiply_value:1.0\n')
            outy.write(f'pivot_point:0\n')
            outy.write(f'pivot_x:0\n')
            outy.write(f'pivot_visible:False\n')



def main():
    parser = argparse.ArgumentParser(description="SpinExplorer")
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to a YAML config file for sorting batch analyses"
    )
    parser.add_argument(
        "--organise-by",
        dest="organise_by",
        choices=["protein", "experiment", "both", "protein+compound", "compound"],
        default=None,
        help=(
            "Create organised output folders for sessions. "
            "'protein' groups by protein name, "
            "'experiment' groups by pulse sequence, "
            "'both' nests experiment folders inside protein folders."
        )
    )

    args = parser.parse_args()

    parent_folder = Path.cwd()
    child_folders = [f for f in parent_folder.iterdir() if f.is_dir()]

    if not child_folders:
        print("No child folders found to process data")
        return

    print(f"Found {len(child_folders)} folders to process")

    for folder in child_folders:
        print(f"\nProcessing: {folder.name}")
        try:
            os.chdir(folder)

            input_dat = FindingParameters()
            pp_parser = PulseSequenceParser()
            sequence = pp_parser.parse()

            config = registry.get_default_config(sequence)

            nmr_glue_conv = Convert_nmrglue(input_dat.params, input_dat)

            params = parameter_write_cl(nmr_glue_conv, config)
            params.write_out_dict(params.dictionary)

            if(nmr_glue_conv.params.remove_filter_before_processing==True):
                remove_filter=False
            else:
                remove_filter=True


            config.process_data(pseudo_flag=nmr_glue_conv.params.pseudo_flag, filter_removal=remove_filter)

            print(f"Successfully processed: {folder.name}")

        except Exception as e:
            print(f"Processing not possible for: {folder.name} ({e})")
            print(f"Make sure the folder contains NMR data and processing instructions")
            print(f"are in the registry for this pulse sequence")
        finally:
            os.chdir(parent_folder)

    if args.config:
        process_from_config(args.config, organise_by=args.organise_by)


def _resolve_output_dir(
    base_dir: Path,
    organise_by: str | None,
    protein: str,
    experiment: str,
    compound: str,
) -> Path:
    """
    Build and create the output directory for a session based on
    the --organise-by mode.

      protein          → base_dir/sessions/<protein>/
      experiment       → base_dir/sessions/<experiment>/
      compound         → base_dir/sessions/<compound>/
      both             → base_dir/sessions/<protein>/<experiment>/
      protein+compound → base_dir/sessions/<protein>/<compound>/
      None             → base_dir/  (no subfolder created)
    """
    if organise_by is None:
        return base_dir

    sessions_root = base_dir / "sessions"

    if organise_by == "protein":
        output_dir = sessions_root / protein
    elif organise_by == "experiment":
        output_dir = sessions_root / experiment
    elif organise_by == "compound":
        output_dir = sessions_root / compound
    elif organise_by == "both":
        output_dir = sessions_root / protein / experiment
    elif organise_by == "protein+compound":
        output_dir = sessions_root / protein / compound
    else:
        raise ValueError(f"Unknown organise_by value: {organise_by!r}")

    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir



if __name__ == "__main__":
    main()