from SpinExplorer.SpinExplorer_CL_tools.experiment_config import ExperimentConfigStore, DimensionConfig
from SpinExplorer.SpinExplorer_CL_tools.pulse_sequence_parsing import ConfigurationRegistry

proton_nhsqcconfig = DimensionConfig.standard_proton(ph_p0 = 0.0, zf_additional_value=2, ex_flag = True, ex_start_ppm = 6.0, ex_end_ppm = 9.5)
nitrogen_nhsqcconfig = DimensionConfig.standard_nitrogen(ph_p0=-90.0, zf_additional_value=2)
nitrogen_nhsqcconfig2 = DimensionConfig.nitrogen_alt(ph_p0=0.0, zf_additional_value=2)
nitrogen_btrosyconfig = DimensionConfig.standard_nitrogen(ph_p0=-62,ph_p1=-31, zf_additional_value=2)


registry = ConfigurationRegistry()

nhsqc_config = ExperimentConfigStore(['Dimension 0 (1H)', 'Dimension 1 (15N)'], [proton_nhsqcconfig, nitrogen_nhsqcconfig], 'test.fid', 'test.ft2')
nhsqc_config2 = ExperimentConfigStore(['Dimension 0 (1H)', 'Dimension 1 (15N)'], [proton_nhsqcconfig, nitrogen_nhsqcconfig2], 'test.fid', 'test.ft2')
btrosy_config = ExperimentConfigStore(['Dimension 0 (1H)', 'Dimension 1 (15N)'], [proton_nhsqcconfig, nitrogen_btrosyconfig], 'test.fid', 'test.ft2')
standard_1H_1D = ExperimentConfigStore(['Dimension 0 (1H)'], [DimensionConfig.standard_proton(ph_p0 = 0.0)], 'test.fid', 'test.ft')
waterlogsy_icon = ExperimentConfigStore(['Dimension 0 (1H)'], [DimensionConfig.standard_proton(ph_p0 = 90.0)], 'test.fid', 'test.ft')

# Pulse program names in the auto-process registry (more to be added soon)

# Standard 1D proton experiments
registry.register("zg", standard_1H_1D)
registry.register("zgpr", standard_1H_1D)
registry.register("zgesgp", standard_1H_1D)


registry.register("t1rho.rf", standard_1H_1D)
registry.register("wlogsy.rf", standard_1H_1D)
registry.register("PO-WaterLOGSY.bind", waterlogsy_icon)
registry.register("PO_waterlogsy_icon_260304.bind",waterlogsy_icon)
registry.register("zgesgp_icon_bind_260304.apk",standard_1H_1D)
registry.register("t1rho_icon_bind_260304.apk",standard_1H_1D)

# Diffusion experiments
registry.register("stebpesgp1s", standard_1H_1D)
registry.register("steesgp1s", standard_1H_1D)
registry.register("stebpgp1s", standard_1H_1D)
registry.register("stegp1s", standard_1H_1D)

# 2D 1H-15N experiments
registry.register("hsqcetfpf3gp", nhsqc_config)
registry.register("hsqcfpf3gpphwg.pjs", nhsqc_config2)
registry.register("b_trosyetf3gpsi.3.cw", btrosy_config)
