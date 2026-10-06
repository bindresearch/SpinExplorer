import numpy as np
from numpy.typing import NDArray
from typing import Union, Literal 
from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import nmrglue as ng # type: ignore
from scipy.stats import truncnorm # type: ignore

from typing import Literal, Optional
import numpy as np
from datetime import datetime

from SpinExplorer.SpinProcess.Processing.ist import apply_sampling_schedule_to_2D_signal, generate_sampling_schedule_poisson, generate_sampling_schedule_poisson_nd, ist_2d, read_sched
from SpinExplorer.SpinProcess.Processing.ist import apply_sampling_schedule_to_2D_signal_ist, write_sched, inflate_spectra_2D_signal_ist, apply_sampling_schedule_to_nd_signal_ist, inflate_spectra_nd_signal_ist
from SpinExplorer.SpinProcess.Processing.ist import ist_3d, apply_sampling_schedule_nd



# ── Per-dimension parameters ───────────────────────────────────────────────────

@dataclass
class DimParams:
    """Parameters for a single spectral dimension."""
    label: str       # e.g. '1H', '15N', '13C'
    sw: float        # spectral width (Hz)
    obs: float       # observe frequency (MHz)
    car: float       # carrier (ppm)
    td_size: int     # time domain size (complex points)
    orig: float = 0.0     # origin (Hz)
    apod: float = 0.0
    obs_mid: float = 0.0
    units: float = 0.0
    ft_flag: float = 0.0
    aq_sign: float = 0.0      # 0.0 = direct/States, 2.0 = States-TPPI
    quad_flag: float = 0.0
    off_ppm: float = 0.0
    p0: float = 0.0
    p1: float = 0.0
    apod_code: float = 0.0
    apod_q1: float = 0.0
    apod_q2: float = 0.0
    apod_q3: float = 0.0
    lb: float = 0.0
    gb: float = 0.0
    goff: float = 0.0
    c1: float = 0.0
    apod_df: float = 0.0
    zf: float = 0.0
    x1: float = 0.0
    xn: float = 0.0
    ft_size: float = 0.0
    size: float = 0.0         # FDF3SIZE/FDF4SIZE — full size for indirect dims in 3D/4D

    @property
    def center(self) -> float:
        return self.td_size / 2 + 1


# ── NMRPipe dimension index mapping ───────────────────────────────────────────
#
#   dims[0] (direct)       → FDF2
#   dims[1] (1st indirect) → FDF1
#   dims[2] (2nd indirect) → FDF3
#   dims[3] (3rd indirect) → FDF4
#
DIM_INDEX_TO_PIPE = {0: 2, 1: 1, 2: 3, 3: 4}

# Labels for unused/padding dimensions
UNUSED_LABELS = {1: 'Y', 3: 'Z', 4: 'A'}


# ── Header builder ─────────────────────────────────────────────────────────────

class NMRPipeHeader:
    """
    Constructs an NMRPipe header dictionary from supplied dimension parameters
    and spectral data. Supports 1D–4D spectra.

    Dimension ordering convention (NMRPipe):
        dims[0] → direct dimension  → FDF2
        dims[1] → 1st indirect      → FDF1
        dims[2] → 2nd indirect      → FDF3
        dims[3] → 3rd indirect      → FDF4

    Usage (1D):
        dim = DimParams(label='1H', sw=12820.5, obs=800.56,
                        car=4.86, orig=-2521.3, td_size=16384)
        header = NMRPipeHeader(dims=[dim], data=fid)

    Usage (2D):
        dim1 = DimParams(label='1H',  sw=12820.5, obs=800.56,
                         car=4.86,   orig=-2509.5, td_size=1024)
        dim2 = DimParams(label='15N', sw=1785.1,  obs=81.13,
                         car=117.18, orig=8628.3,  td_size=128,
                         aq_sign=2.0)
        header = NMRPipeHeader(dims=[dim1, dim2], data=spectrum)

    Usage (3D):
        dim1 = DimParams(label='1H',  ...)
        dim2 = DimParams(label='15N', ..., aq_sign=2.0)
        dim3 = DimParams(label='13C', ..., aq_sign=2.0, size=256.0)
        header = NMRPipeHeader(dims=[dim1, dim2, dim3], data=spectrum)
    """

    MAGIC      = 0.0
    FLT_FORMAT = 4008636160.0
    FLT_ORDER  = 2.3450000286102295

    def __init__(self,
                 dims: list[DimParams],
                 data: np.ndarray,
                 nus_dim: float = 0.0,
                 dmx_val: float = 0.0,
                 dmx_flag: float = 0.0,
                 title: str = '',
                 comment: str = '',
                 username: str = '',
                 operator: str = '',
                 src_name: str = ''):

        if not 1 <= len(dims) <= 4:
            raise ValueError("dims must have between 1 and 4 entries")
        if data.ndim != len(dims):
            raise ValueError(
                f"data has {data.ndim} dimensions but {len(dims)} DimParams supplied")

        self.dims     = dims
        self.data     = data
        self.n_dims   = len(dims)
        self.now      = datetime.now()
        self.nus_dim  = nus_dim
        self.dmx_val  = dmx_val
        self.dmx_flag = dmx_flag

        self.title    = title
        self.comment  = comment
        self.username = username
        self.operator = operator
        self.src_name = src_name

    # ── Public ────────────────────────────────────────────────────────────────

    def build(self) -> dict:
        h = {}
        h.update(self._global_fields())
        h.update(self._dim_order_fields())
        h.update(self._data_fields())
        h.update(self._datetime_fields())
        h.update(self._text_fields())

        # Active dimensions
        for i, dim in enumerate(self.dims):
            pipe_idx = DIM_INDEX_TO_PIPE[i]
            h.update(self._dim_fields(dim, pipe_idx))

        # Padding dimensions
        active = {DIM_INDEX_TO_PIPE[i] for i in range(self.n_dims)}
        for pipe_idx in [2, 1, 3, 4]:
            if pipe_idx not in active:
                h.update(self._unused_dim_fields(pipe_idx))

        return h

    def update_dim(self, index: int, **kwargs) -> None:
        """Update fields on a dimension. index is 0-based (0 = direct)."""
        dim = self.dims[index]
        for k, v in kwargs.items():
            if hasattr(dim, k):
                setattr(dim, k, v)
            else:
                raise ValueError(f"DimParams has no field '{k}'")

    # ── Private builders ──────────────────────────────────────────────────────

    def _fdspecnum(self) -> float:
        """
        FDSPECNUM = number of FIDs (rows) in the file.
        1D: 1
        2D+: first indirect dimension td_size * 2 (States/hypercomplex factor)
        """
        if self.n_dims == 1:
            return 1.0
        return float(self.dims[1].td_size * 2)

    def _global_fields(self) -> dict:
        is_1d = self.n_dims == 1
        is_2d = self.n_dims == 2
        return {
            'FDMAGIC':       self.MAGIC,
            'FDFLTFORMAT':   self.FLT_FORMAT,
            'FDFLTORDER':    self.FLT_ORDER,
            'FDSIZE':        float(self.data.shape[-1]),
            'FDREALSIZE':    float(self.data.shape[-1]),
            'FDSPECNUM':     self._fdspecnum(),
            'FDQUADFLAG':    0.0,
            'FD2DPHASE':     0.0 if is_1d else 2.0,
            'FDTRANSPOSED':  0.0,
            'FDDIMCOUNT':    float(self.n_dims),
            'FDNUSDIM':      self.nus_dim,
            'FDPIPEFLAG': 0.0 if is_2d else 1.0,
            'FDCUBEFLAG':    0.0,
            'FDPIPECOUNT':   0.0,
            'FDSLICECOUNT':  0.0,
            'FDSLICECOUNT1': 0.0,
            'FDFILECOUNT':   1.0,
            'FDTHREADCOUNT': 0.0,
            'FDTHREADID':    0.0,
            'FDFIRSTPLANE':  0.0,
            'FDLASTPLANE':   0.0,
            'FDPARTITION':   0.0,
            'FDPLANELOC':    0.0,
            'FDMCFLAG':      0.0,
            'FDNOISE':       0.0,
            'FDRANK':        0.0,
            'FDTEMPERATURE': 0.0,
            'FDPRESSURE':    0.0,
            'FD2DVIRGIN':    1.0,
            'FDTAU':         0.0,
            'FDDOMINFO':     0.0,
            'FDMETHINFO':    0.0,
            'FDSCORE':       0.0,
            'FDSCANS':       0.0,
            'FDDMXVAL':      self.dmx_val,
            'FDDMXFLAG':     self.dmx_flag,
            'FDDELTATR':     0.0,
            'FDLASTBLOCK':   0.0,
            'FDCONTBLOCK':   0.0,
            'FDBASEBLOCK':   0.0,
            'FDPEAKBLOCK':   0.0,
            'FDBMAPBLOCK':   0.0,
            'FDHISTBLOCK':   0.0,
            'FD1DBLOCK':     0.0,
            'FDUSER1':       0.0,
            'FDUSER2':       0.0,
            'FDUSER3':       0.0,
            'FDUSER4':       0.0,
            'FDUSER5':       0.0,
            'FDUSER6':       0.0,
        }

    def _dim_order_fields(self) -> dict:
        order = [2.0, 1.0, 3.0, 4.0]
        return {
            'FDDIMORDER':  order,
            'FDDIMORDER1': order[0],
            'FDDIMORDER2': order[1],
            'FDDIMORDER3': order[2],
            'FDDIMORDER4': order[3],
        }

    def _data_fields(self) -> dict:
        real_data = self.data.real if np.iscomplexobj(self.data) else self.data
        has_data  = real_data.size > 0 and real_data.max() != real_data.min()
        return {
            'FDMAX':       float(real_data.max()) if has_data else 0.0,
            'FDMIN':       float(real_data.min()) if has_data else 0.0,
            'FDSCALEFLAG': 1.0 if has_data else 0.0,
            'FDDISPMAX':   float(real_data.max()) if has_data else 0.0,
            'FDDISPMIN':   float(real_data.min()) if has_data else 0.0,
            'FDPTHRESH':   0.0,
            'FDNTHRESH':   0.0,
        }

    def _datetime_fields(self) -> dict:
        return {
            'FDMONTH': float(self.now.month),
            'FDDAY':   float(self.now.day),
            'FDYEAR':  float(self.now.year),
            'FDHOURS': float(self.now.hour),
            'FDMINS':  float(self.now.minute),
            'FDSECS':  float(self.now.second),
        }

    def _text_fields(self) -> dict:
        return {
            'FDSRCNAME':  self.src_name,
            'FDUSERNAME': self.username,
            'FDOPERNAME': self.operator,
            'FDTITLE':    self.title,
            'FDCOMMENT':  self.comment,
        }

    def _dim_fields(self, dim: DimParams, pipe_idx: int) -> dict:
        p = f'FDF{pipe_idx}'
        fields = {
            f'{p}LABEL':    dim.label,
            f'{p}APOD':     float(dim.apod or dim.td_size),
            f'{p}SW':       dim.sw,
            f'{p}OBS':      dim.obs,
            f'{p}OBSMID':   dim.obs_mid,
            f'{p}ORIG':     dim.orig,
            f'{p}UNITS':    dim.units,
            f'{p}FTFLAG':   dim.ft_flag,
            f'{p}AQSIGN':   dim.aq_sign,
            f'{p}QUADFLAG': dim.quad_flag,
            f'{p}CAR':      dim.car,
            f'{p}CENTER':   float(dim.center),
            f'{p}OFFPPM':   dim.off_ppm,
            f'{p}P0':       dim.p0,
            f'{p}P1':       dim.p1,
            f'{p}APODCODE': dim.apod_code,
            f'{p}APODQ1':   dim.apod_q1,
            f'{p}APODQ2':   dim.apod_q2,
            f'{p}APODQ3':   dim.apod_q3,
            f'{p}LB':       dim.lb,
            f'{p}GB':       dim.gb,
            f'{p}GOFF':     dim.goff,
            f'{p}C1':       dim.c1,
            f'{p}APODDF':   dim.apod_df,
            f'{p}ZF':       dim.zf,
            f'{p}X1':       dim.x1,
            f'{p}XN':       dim.xn,
            f'{p}FTSIZE':   dim.ft_size,
            f'{p}TDSIZE':   float(dim.td_size),
        }
        # FDF3SIZE / FDF4SIZE only appear on higher indirect dimensions
        if pipe_idx in (3, 4):
            fields[f'{p}SIZE'] = dim.size if dim.size > 0 else float(dim.td_size * 2)
        return fields

    def _unused_dim_fields(self, pipe_idx: int) -> dict:
        p     = f'FDF{pipe_idx}'
        label = UNUSED_LABELS.get(pipe_idx, 'X')
        return {
            f'{p}LABEL':    label,
            f'{p}APOD':     0.0,
            f'{p}OBS':      0.0,
            f'{p}OBSMID':   0.0,
            f'{p}SW':       0.0,
            f'{p}ORIG':     0.0,
            f'{p}FTFLAG':   0.0,
            f'{p}AQSIGN':   0.0,
            f'{p}SIZE':     1.0,
            f'{p}QUADFLAG': 1.0,
            f'{p}UNITS':    0.0,
            f'{p}P0':       0.0,
            f'{p}P1':       0.0,
            f'{p}CAR':      0.0,
            f'{p}CENTER':   1.0,
            f'{p}OFFPPM':   0.0,
            f'{p}APODCODE': 0.0,
            f'{p}APODQ1':   0.0,
            f'{p}APODQ2':   0.0,
            f'{p}APODQ3':   0.0,
            f'{p}LB':       0.0,
            f'{p}GB':       0.0,
            f'{p}GOFF':     0.0,
            f'{p}C1':       0.0,
            f'{p}ZF':       0.0,
            f'{p}X1':       0.0,
            f'{p}XN':       0.0,
            f'{p}FTSIZE':   0.0,
            f'{p}TDSIZE':   0.0,
        }

@dataclass
class DistributionParam:
    min_val: float
    max_val: float
    dist: Literal["uniform", "normal"] = "uniform"
    dtype: Literal["float", "int"] = "float"
    mean: float | None = None
    std: float | None = None
    seed: int | None = None
    _rng: np.random.Generator = field(init=False, repr=False)

    def __post_init__(self):
        self._rng = np.random.default_rng(self.seed)
        if self.dist == "normal" and (self.mean is None or self.std is None):
            raise ValueError("mean and std are required for normal distribution")

    def draw(self, n: int) -> np.ndarray:
        if self.dist == "uniform":
            samples = self._rng.uniform(low=self.min_val, high=self.max_val, size=n)
        else:
            samples = truncated_normal(self.mean, self.std, self.min_val, self.max_val, n, seed = self.seed)
        return samples.astype(int) if self.dtype == "int" else samples

    def single(self) -> int | float:
        return self.draw(1)[0]

def truncated_normal(mean, std, lower, upper, n, seed=None):
    a = (lower - mean) / std
    b = (upper - mean) / std
    return truncnorm.rvs(a, b, loc=mean, scale=std, size=n, random_state=seed)

def freq_to_rad(val):
    return val*2.0*np.pi

def rad_to_freq(val):
    return val/(2.0*np.pi)


def generate_signal(amp: float,freq: float, r2:float, max_time:float, num_points:int)->NDArray:
    times = np.linspace(0,max_time,num_points)
    return amp*np.exp((1.0j*freq-r2)*times)

# Vectorised signal generation — shape (num_signals, num_points)
def general_signals_vector(amp, freq, r2, times):
    return (amp[:, None] * np.exp(
            (1j * freq[:, None] - r2[:, None]) * times
            )).astype(np.complex64)

def generate_1d_spec_from_ranges(num_signals:int, amp_dist: DistributionParam, freq_dist: DistributionParam, 
                                 r2_dist: DistributionParam, max_time:float, num_points:int):
    

    amp_vals,freq_vals,r2_vals = [d.draw(num_signals) for d in (amp_dist,freq_dist,r2_dist)]
    return np.sum([generate_signal(amp,freq,r2, max_time, num_points) for (amp,freq,r2) in zip(amp_vals,freq_vals,r2_vals)], axis = 0)

def generate_2d_spec_from_ranges(num_signals:int, amp_dist1: DistributionParam, freq_dist1: DistributionParam, r2_dist1: DistributionParam,
                                 sw1:float, num_points1: int,
                                 amp_dist2: DistributionParam, freq_dist2: DistributionParam, r2_dist2: DistributionParam,
                                 sw2:float, num_points2: int):
    
    amp_vals1,freq_vals1,r2_vals1 = [d.draw(num_signals) for d in (amp_dist1,freq_dist1,r2_dist1)]
    amp_vals2,freq_vals2,r2_vals2 = [d.draw(num_signals) for d in (amp_dist2,freq_dist2,r2_dist2)]
    
    max_time1 = (num_points1 - 1) / sw1
    max_time2 = (num_points2 - 1) / sw2
    
    times1 = np.linspace(0, max_time1, num_points1, dtype=np.float32)
    times2 = np.linspace(0, max_time2, num_points2, dtype=np.float32)

    sigs_direct    = general_signals_vector(amp_vals1, freq_vals1, r2_vals1, times1)  # (num_signals, num_points1)
    sigs_indirect  = general_signals_vector(amp_vals2, freq_vals2, r2_vals2, times2)  # (num_signals, num_points2)

    cos_indirect = sigs_indirect.real  # (num_signals, num_points2)
    sin_indirect = sigs_indirect.imag  # (num_signals, num_points2)

    # Vectorised outer product summed over signals — shape (num_points2, num_points1)
    fid_2d_r = np.einsum('si,sj->ij', cos_indirect, sigs_direct, optimize=True)
    fid_2d_i = np.einsum('si,sj->ij', sin_indirect, sigs_direct, optimize=True)

    n_indirect, n_direct = fid_2d_r.shape
    fid_interleaved = np.empty((n_indirect * 2, n_direct), dtype=np.complex64)
    fid_interleaved[0::2] = fid_2d_r  # cos-modulated rows
    fid_interleaved[1::2] = fid_2d_i  # sin-modulated rows

    # fid_2d_r = np.sum([
    #     np.outer(np.real(generate_signal(amp2, freq2, r2_2, max_time2, int(num_points2))),
    #              generate_signal(amp1, freq1, r2_1, max_time1, num_points1))
    #     for (amp1, freq1, r2_1, amp2, freq2, r2_2) 
    #     in zip(amp_vals1, freq_vals1, r2_vals1, amp_vals2, freq_vals2, r2_vals2)
    # ], axis=0)  # shape: (num_points2, num_points1)

    # fid_2d_i = np.sum([
    #     np.outer(np.imag(generate_signal(amp2, freq2, r2_2, max_time2, int(num_points2))),
    #              generate_signal(amp1, freq1, r2_1, max_time1, num_points1))
    #     for (amp1, freq1, r2_1, amp2, freq2, r2_2) 
    #     in zip(amp_vals1, freq_vals1, r2_vals1, amp_vals2, freq_vals2, r2_vals2)
    # ], axis=0)  

    # # Hypercomplex combination along indirect dimension (axis=0)
    # # Split real/imag of indirect dimension before FT

    # n_indirect, n_direct = fid_2d_r.shape

    # # NMRPipe expects each row as interleaved [re0, im0, re1, im1, ...]
    # # so the direct dimension doubles in size on disk

    # fid_interleaved = np.empty((n_indirect * 2, n_direct), dtype=np.complex64)
    # fid_interleaved[0::2] = fid_2d_r.astype(np.complex64)  # cos-modulated rows
    # fid_interleaved[1::2] = fid_2d_i.astype(np.complex64)  # sin-modulated rows

    return fid_interleaved


def generate_3d_spec_from_ranges(num_signals: int,
                                 amp_dist1: DistributionParam, freq_dist1: DistributionParam, r2_dist1: DistributionParam,
                                 sw1: float, num_points1: int,
                                 amp_dist2: DistributionParam, freq_dist2: DistributionParam, r2_dist2: DistributionParam,
                                 sw2: float, num_points2: int,
                                 amp_dist3: DistributionParam, freq_dist3: DistributionParam, r2_dist3: DistributionParam,
                                 sw3: float, num_points3: int):

    amp_vals1, freq_vals1, r2_vals1 = [d.draw(num_signals) for d in (amp_dist1, freq_dist1, r2_dist1)]
    amp_vals2, freq_vals2, r2_vals2 = [d.draw(num_signals) for d in (amp_dist2, freq_dist2, r2_dist2)]
    amp_vals3, freq_vals3, r2_vals3 = [d.draw(num_signals) for d in (amp_dist3, freq_dist3, r2_dist3)]

    max_time1 = (num_points1 - 1) / sw1
    max_time2 = (num_points2 - 1) / sw2
    max_time3 = (num_points3 - 1) / sw3

    # Four hypercomplex components for 3D States:
    # indirect dim 2 (slowest) × indirect dim 1 × direct
    fid_rr = np.zeros((num_points3, num_points2, num_points1), dtype=complex)  # cos2 × cos1 × direct
    fid_ri = np.zeros((num_points3, num_points2, num_points1), dtype=complex)  # cos2 × sin1 × direct
    fid_ir = np.zeros((num_points3, num_points2, num_points1), dtype=complex)  # sin2 × cos1 × direct
    fid_ii = np.zeros((num_points3, num_points2, num_points1), dtype=complex)  # sin2 × sin1 × direct


    times1 = np.linspace(0, max_time1, num_points1, dtype=np.float32)
    times2 = np.linspace(0, max_time2, num_points2, dtype=np.float32)
    times3 = np.linspace(0, max_time3, num_points3, dtype=np.float32)

    sig_direct = general_signals_vector(amp_vals1, freq_vals1, r2_vals1, times1)
    sig_indirect1 = general_signals_vector(amp_vals2, freq_vals2, r2_vals2, times2)
    sig_indirect2 = general_signals_vector(amp_vals3, freq_vals3, r2_vals3, times3)

    cos1 = np.real(sig_indirect1)  
    sin1 = np.imag(sig_indirect1)  
    cos2 = np.real(sig_indirect2)  
    sin2 = np.imag(sig_indirect2)  

    fid_rr = np.einsum('si,sj,sk->ijk', cos2, cos1, sig_direct, optimize = True)
    fid_ri = np.einsum('si,sj,sk->ijk', cos2, sin1, sig_direct, optimize = True)
    fid_ir = np.einsum('si,sj,sk->ijk', sin2, cos1, sig_direct, optimize = True)
    fid_ii = np.einsum('si,sj,sk->ijk', sin2, sin1, sig_direct, optimize = True)

    # NMRPipe hypercomplex interleaving for 3D:
    # The two indirect dimensions are interleaved as pairs of rows.
    # For each indirect2 point, there are 4 rows:
    #   [cos2_cos1, cos2_sin1, sin2_cos1, sin2_sin1]
    # shape: (num_points3 * 4, num_points2, num_points1) — no wait,
    # NMRPipe streams 3D as a series of 2D planes, so we interleave
    # along the slowest indirect dim (axis=0) with factor 2,
    # and along the faster indirect dim (axis=1) with factor 2.

    n3, n2, n1 = fid_rr.shape

    # Interleave faster indirect (dim2/axis=1): cos1/sin1 pairs per row
    # giving shape (n3, n2*2, n1)
    fid_indirect1_interleaved = np.empty((n3, n2 * 2, n1), dtype=np.complex64)
    fid_indirect1_interleaved[:, 0::2, :] = fid_rr.astype(np.complex64)  # cos2_cos1
    fid_indirect1_interleaved[:, 1::2, :] = fid_ri.astype(np.complex64)  # cos2_sin1

    fid_indirect1_interleaved_sin = np.empty((n3, n2 * 2, n1), dtype=np.complex64)
    fid_indirect1_interleaved_sin[:, 0::2, :] = fid_ir.astype(np.complex64)  # sin2_cos1
    fid_indirect1_interleaved_sin[:, 1::2, :] = fid_ii.astype(np.complex64)  # sin2_sin1

    # Interleave slowest indirect (dim3/axis=0): cos2/sin2 pairs
    # giving final shape (n3*2, n2*2, n1)
    fid_interleaved = np.empty((n3 * 2, n2 * 2, n1), dtype=np.complex64)
    fid_interleaved[0::2, :, :] = fid_indirect1_interleaved        # cos2 rows
    fid_interleaved[1::2, :, :] = fid_indirect1_interleaved_sin    # sin2 rows

    return fid_interleaved


def write_as_nmrpipe(array: NDArray, spec_dic: dict, outfile: str):
    #print(array.shape)
    ng.pipe.write(outfile, spec_dic, array, overwrite=True)

def direct_dimension_process(infile: str):
    dic, data = ng.pipe.read(infile)
    dic, data = ng.pipe_proc.em(dic,data, lb = 5.0, c = 0.5)
    dic, data = ng.pipe_proc.zf(dic, data, zf = 1)
    dic, data = ng.pipe_proc.ft(dic, data)
    dic, data = ng.pipe_proc.ps(dic, data, p0=0.0, p1=0.0)
    dic, data = ng.pipe_proc.di(dic, data)
    dic, data = ng.pipe_proc.tp(dic, data)

    return dic, data

def direct_dimension_process_3d(infile: str):
    dic, data = ng.pipe.read(infile)
    dic, data = ng.pipe_proc.em(dic,data, lb = 5.0, c = 0.5)
    dic, data = ng.pipe_proc.zf(dic, data, zf = 1)
    dic, data = ng.pipe_proc.ft(dic, data)
    dic, data = ng.pipe_proc.ps(dic, data, p0=0.0, p1=0.0)
    dic, data = ng.pipe_proc.di(dic, data)

    return dic, data


# import time 
# start = time.time()

sig_dist = DistributionParam(50,100,'uniform','int')
amp_dist = DistributionParam(0.2,1.0,'normal', 'float', 0.5, 0.2)
freq_dist = DistributionParam(-5000.0*2.0*np.pi,5000.0*2.0*np.pi,'uniform','float')
freq_dist2 = DistributionParam(-900.0*2.0*np.pi,900.0*2.0*np.pi,'uniform','float')
freq_dist3 = DistributionParam(-1000.0*2.0*np.pi, 1000.0*2.0*np.pi, 'uniform','float')
r2_dist = DistributionParam(10.0,20.0,'normal','float',12.0,4.0)
sw1 = 12500.
sw2 = 2000.0
sw3 = 3000.0
num_points1 = 512
num_points2 = 20
num_points3 = 10
num_sigs = int(sig_dist.single())



num_points1 = 64
num_points2 = 40
num_points3 = 20
trial_3d = generate_3d_spec_from_ranges(num_sigs, amp_dist,freq_dist,r2_dist,sw1,num_points1,
                                         amp_dist,freq_dist2,r2_dist,sw2,num_points2,
                                         amp_dist, freq_dist3, r2_dist,sw3, num_points3)

dim1 = DimParams(label='1H',  sw=sw1, obs=800.56,
                    car=4.86, td_size=num_points1)
dim2 = DimParams(label='15N', sw=sw2,  obs=81.13,
                    car=117.18,  td_size=num_points2,
                    aq_sign=0.0)
dim3 = DimParams(label='13C', sw=sw3,  obs=200.14,
                    car=35.0,  td_size=num_points3,
                    aq_sign=0.0)

sampling_sched = generate_sampling_schedule_poisson_nd(0.25, [num_points3,num_points2], 2)
write_sched(sampling_sched,'sched')
header = NMRPipeHeader(dims=[dim1, dim2, dim3], data=trial_3d)
dic = header.build()
sampling_sched = read_sched('sched')
write_as_nmrpipe(trial_3d, dic, 'test_3d.fid')

#nus_data, nus_dic = apply_sampling_schedule_to_nd_signal_ist(trial_3d,dic, sampling_sched)
#write_as_nmrpipe(nus_data, nus_dic, 'hms_3d.fid')
#print(nus_data.shape)
#nus_data, nus_dic = inflate_spectra_nd_signal_ist(nus_data, nus_dic, sampling_sched, [num_points3, num_points2])


nus_dic,nus_data = direct_dimension_process_3d('test_3d.fid')
nus_data = apply_sampling_schedule_nd(nus_data, sampling_sched,indirect_axes=(0,1))

write_as_nmrpipe(nus_data, nus_dic, 'nus_3d.ft1')

nus_data = np.transpose(nus_data, axes=(2,0,1))

recon_data = ist_3d(nus_data, sampling_sched,mode=1,max_iter = 1000)
recon_data = np.transpose(recon_data, axes=(1,2,0))

ng.pipe.write('nus_3d_recon.ft1', nus_dic, recon_data, True)



def check_3d_proc(array_3d:NDArray, indirect_axes):
    signal_ft = np.fft.fftn(array_3d, axes=indirect_axes)
    return signal_ft

def plot_2d_nmr_contour(spectrum: NDArray,
                        min_contour_fraction: float = 0.1,
                        n_levels: int = 10,
                        fac: float = 1.2,
                        pos_color: str = 'blue',
                        neg_color: str = 'red',
                        ax: Optional[plt.Axes] = None,
                        figsize: tuple = (8, 8)) -> plt.Axes:
    """
    Contour plot of 2D NMR data with geometrically spaced contour levels.

    Parameters
    ----------
    spectrum              : 2D numpy array of spectral intensities
    min_contour_fraction  : minimum contour level as fraction of max intensity (default 0.1)
    n_levels              : number of contour levels for positive and negative (default 10)
    fac                   : multiplicative spacing between contour levels (default 1.2)
    pos_color             : colour for positive contours (default 'blue')
    neg_color             : colour for negative contours (default 'red')
    ax                    : existing matplotlib axes to plot on (default None — creates new)
    figsize               : figure size if creating new axes (default (8, 8))

    Returns
    -------
    ax : matplotlib Axes object
    """
    max_val = np.max(np.abs(spectrum))
    cl_start = min_contour_fraction * max_val

    # Geometrically spaced contour levels: cl_start * fac^0, fac^1, ..., fac^(n-1)
    levels = cl_start * fac ** np.arange(n_levels)

    pos_levels = levels          #  positive contours
    neg_levels = -levels[::-1]   #  negative contours (mirrored, highest first)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    ax.contour(spectrum,  levels=pos_levels, colors=pos_color, linewidths=0.8)
    ax.contour(spectrum,  levels=neg_levels, colors=neg_color, linewidths=0.8)

    ax.set_xlabel('Direct dimension')
    ax.set_ylabel('Indirect dimension')

    return ax

#dic, data = ng.pipe.read('test_3d.ft1')
#dic, data = ng.pipe_proc.tp(dic, data)

# print(data.shape)
# signal_ft = check_3d_proc(data, indirect_axes=(2,))
# signal_ft = np.transpose(signal_ft, axes = (0,2,1))
# print(signal_ft.shape)
# signal_ft = check_3d_proc(signal_ft, indirect_axes=(2,))


"""
dic, data = ng.pipe.read('test_3d.fid')
#dic,data = ng.pipe_proc.tp(dic,data)

lb = 10.0

apod = np.exp(-np.pi * np.arange(data.shape[-1]) * lb/sw1).astype(data.dtype)

data[...,:] = data[...,:]*apod
data = np.fft.fftshift(np.fft.fft(data), axes = -1)
data = np.real(data)
plt.plot(data[0,0,:])
plt.show()
data = np.transpose(data, axes=(2,0,1))
data = data[...,0::2] + 1.0j*data[...,1::2]

apod = np.exp(-np.pi * np.arange(data.shape[-1]) * lb/sw2).astype(data.dtype)
data[...,:] = data[...,:]*apod

n_zeros = 128-data.shape[-1]
zeros = np.zeros((*data.shape[:-1], n_zeros), dtype=data.dtype)
data = np.concatenate([data, zeros], axis=-1)

data = np.fft.fftshift(np.fft.fft(data), axes = -1)
data = np.real(data)
print(data.shape)

data = np.transpose(data, axes=(0,2,1))

apod = np.exp(-np.pi * np.arange(data.shape[-1]) * lb/sw3).astype(data.dtype)
data[...,:] = data[...,:]*apod

n_zeros = 128-data.shape[-1]
zeros = np.zeros((*data.shape[:-1], n_zeros), dtype=data.dtype)
data = np.concatenate([data, zeros], axis=-1)

data = data[...,0::2] + 1.0j*data[...,1::2]
data = np.fft.fftshift(np.fft.fft(data), axes = -1)
#data = np.real(data)
print(data.shape)

for i in range(data.shape[0]):
    ax = plot_2d_nmr_contour(np.real(data[i,:,:]), min_contour_fraction=0.1, n_levels=12, fac = 1.2)
    plt.show()
"""