from experiment_config import ExperimentConfigStore, DimensionConfig
proton_config = DimensionConfig.standard_proton(ph_p0 = 0.0, ex_flag = True, ex_start_ppm = 6.0, ex_end_ppm = 9.5)
nitrogen_config = DimensionConfig.standard_nitrogen(ph_p0=90.0)

hsqc_config_template = ExperimentConfigStore(['Dimension 0 (1H)', 'Dimension 1 (15N)'], [proton_config, nitrogen_config])