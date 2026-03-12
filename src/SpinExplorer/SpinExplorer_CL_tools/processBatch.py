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
    for base_title in base_titles:
        base_mask = df["Title"] == base_title
        compound_mask = (
            df["Title"].str.contains(base_title, regex=False) &
            df["Title"].str.contains(protein, case=False)
        )
        group = df[base_mask | compound_mask]
        if len(group) > 1:
            groups.append(group)
    
    return groups

def process_from_config(config_path: str) -> None:
    with open(config_path, "r") as f:
        config_data = yaml.safe_load(f)
    
    df = pd.read_csv(config_data["csv_path"])
    
    for exp in config_data["experiments"]:
        config = registry._registry[exp]
        for prot in config_data["proteins"]:
            filtered_exps = filter_experiments(df, exp, prot)
            groups = group_by_base_title(filtered_exps, prot)
            for group in groups:
                write_1d_multi_session(group, exp, prot, config)

def write_1d_multi_session(df, exp, protein_name, config, outy_name = None, outy_folder = Path('./')):
    if outy_name is None:
        titles = df["Title"]
        base = [t for t in titles if " " not in t][0]
        outy_name = base+'_'+protein_name+'_'+exp+'.session'
    
    with open(outy_folder / Path(outy_name), 'w') as outy:
        outy.write('1D\n')
        outy.write('MultiplotMode:True\n')
        for i, (_, row) in enumerate(df.iterrows()):
            outy.write(f'file_path:{str(Path.cwd())+'/'+str(row['Expno'])+'/test.ft'}\n')
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

    args = parser.parse_args()

    parent_folder = Path.cwd()
    print('hello')
    print(parent_folder)
    child_folders = [f for f in parent_folder.iterdir() if f.is_dir()]
    
    if not child_folders:
        print("No child folders found")
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

            config.process_data()
            
            print(f"Successfully processed: {folder.name}")
            
        except Exception as e:
            print(f"Processing not possible for: {folder.name} ({e})")
            print(f"Make sure the folder contains NMR data and processing instructions")
            print(f"are in the registry for this pulse sequence")
        finally:
            os.chdir(parent_folder)
        
    if args.config:
        process_from_config(args.config)
    

if __name__ == "__main__":
    main()