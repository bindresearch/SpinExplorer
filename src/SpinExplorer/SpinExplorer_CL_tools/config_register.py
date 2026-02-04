from experiment_config import ExperimentConfigStore, DimensionConfig
from pulse_sequence_parsing import ConfigurationRegistry

proton_config = DimensionConfig.standard_proton(ph_p0 = 0.0, ex_flag = True, ex_start_ppm = 6.0, ex_end_ppm = 9.5)
nitrogen_config = DimensionConfig.standard_nitrogen(ph_p0=90.0)

registry = ConfigurationRegistry()

hsqc_config = ExperimentConfigStore(['Dimension 0 (1H)', 'Dimension 1 (15N)'], [proton_config, nitrogen_config], 'test.fid', 'test.ft2')
standard_1H_1D = ExperimentConfigStore(['Dimension 0 (1H)'], [DimensionConfig.standard_proton(ph_p0 = 0.0)], 'test.fid', 'test.ft')

registry.register("hsqcetfpf3gp", hsqc_config)
registry.register("zgesgp", standard_1H_1D)
registry.register("t1rho.rf", standard_1H_1D)
registry.register("wlogsy.rf", standard_1H_1D)