from SpinExplorer.SpinExplorer_CL_tools.experiment_config import ExperimentConfigStore, DimensionConfig, FTOptions
from SpinExplorer.SpinExplorer_CL_tools.pulse_sequence_parsing import ConfigurationRegistry

proton_nhsqcconfig = DimensionConfig.standard_proton(ph_p0 = 0.0, zf_additional_value=2, ex_flag = True, ex_start_ppm = 6.0, ex_end_ppm = 9.5)
nitrogen_nhsqcconfig = DimensionConfig.standard_nitrogen(ph_p0=-90.0, zf_additional_value=2)
nitrogen_nhsqcconfig2 = DimensionConfig.nitrogen_alt(ph_p0=0.0, zf_additional_value=2)
nitrogen_btrosyconfig = DimensionConfig.standard_nitrogen(ph_p0=-62,ph_p1=-31, zf_additional_value=2)

nitrogen_dim_alt = DimensionConfig.standard_nitrogen(ft_option = FTOptions.ALT)
nitrogen_dim_alt_neg = DimensionConfig.standard_nitrogen(ft_option = FTOptions.ALT_NEG)

carbon_dim_alt = DimensionConfig.standard_carbon(ft_option = FTOptions.ALT)
carbon_dim_alt_neg = DimensionConfig.standard_carbon(ft_option = FTOptions.ALT_NEG)

registry = ConfigurationRegistry()

nhsqc_config = ExperimentConfigStore(['Dimension 0 (1H)', 'Dimension 1 (15N)'], [proton_nhsqcconfig, nitrogen_nhsqcconfig], 'test.fid', 'test.ft2')
nhsqc_config2 = ExperimentConfigStore(['Dimension 0 (1H)', 'Dimension 1 (15N)'], [proton_nhsqcconfig, nitrogen_nhsqcconfig2], 'test.fid', 'test.ft2')
btrosy_config = ExperimentConfigStore(['Dimension 0 (1H)', 'Dimension 1 (15N)'], [proton_nhsqcconfig, nitrogen_btrosyconfig], 'test.fid', 'test.ft2')
standard_1H_1D = ExperimentConfigStore(['Dimension 0 (1H)'], [DimensionConfig.standard_proton(ph_p0 = 0.0)], 'test.fid', 'test.ft')
waterlogsy_icon = ExperimentConfigStore(['Dimension 0 (1H)'], [DimensionConfig.standard_proton(ph_p0 = 90.0)], 'test.fid', 'test.ft')

standard_19F_1D = ExperimentConfigStore(['Dimension 0 (19F)'], [DimensionConfig.standard_fluorine(ph_p0 = 0.0)], 'test.fid', 'test.ft')
triple_res_bruker = ExperimentConfigStore(['Dimension 0 (1H)', 'Dimension 1 (15N)', 'Dimension 2 (13C)'],[proton_nhsqcconfig, nitrogen_dim_alt_neg, carbon_dim_alt], 'test.fid', 'test.ft3')
double_N_assignment_bruker = ExperimentConfigStore(['Dimension 0 (1H)', 'Dimension 1 (15N)', 'Dimension 2 (15N)'],[proton_nhsqcconfig, nitrogen_dim_alt, nitrogen_dim_alt], 'test.fid', 'test.ft3')


# Pulse program names in the auto-process registry (more to be added soon)

# Standard 1D proton experiments
registry.register("zg", standard_1H_1D)
registry.register("zgpr", standard_1H_1D)
registry.register("zgesgp", standard_1H_1D)
registry.register("zgesgp.apk", standard_1H_1D)


registry.register("t1rho.rf", standard_1H_1D)
registry.register("t1rho.apk", standard_1H_1D)
registry.register("wlogsy.rf", standard_1H_1D)
registry.register("PO_waterlogsy.bind", waterlogsy_icon)
registry.register("PO-WaterLOGSY.bind", waterlogsy_icon)


# Diffusion experiments
registry.register("stebpesgp1s", standard_1H_1D)
registry.register("steesgp1s", standard_1H_1D)
registry.register("stebpgp1s", standard_1H_1D)
registry.register("stegp1s", standard_1H_1D)


# Fluorine Experiments
registry.register("19F_R2_cpmg_500_dfh_2", standard_19F_1D)

# 2D 1H-15N experiments
registry.register("hsqcetfpf3gp", nhsqc_config)
registry.register("hsqcfpf3gpphwg.pjs", nhsqc_config2)
registry.register("b_trosyetf3gpsi.3.cw", btrosy_config)

# 3D assignment experiments
registry.register("hncogpwg3d", triple_res_bruker)
registry.register("hncacogpwg3d", triple_res_bruker)
registry.register("cbcaconhgpwg3d", triple_res_bruker)
registry.register("cbcanhgpwg3d", triple_res_bruker)
registry.register("hncagpwg3d", triple_res_bruker)
registry.register("hncocagpwg3d", triple_res_bruker)
registry.register("hncocacbgpwg3d", triple_res_bruker)

registry.register("hncannhgpwg3d", double_N_assignment_bruker)
registry.register("hncocannhgpwg3d", double_N_assignment_bruker)