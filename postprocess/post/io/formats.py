# post/io/formats.py
from __future__ import annotations
import numpy as np

# Ajuste aqui quando rodar REAL=double no solver:
REAL_DTYPE = np.dtype("<f4")  # "<f8" se double

# ---- TKE ----
TKE_DTYPE = np.dtype(
    [("t_star", "<i8"), ("ke", "<f8")]
)  # isso é agregado; pode ficar f8

# ---- Centerline ----
CLN_MAGIC = b"CLN1"

# ---- Annular profile ----
APR_MAGIC = b"APR1"
APP_MAGIC = b"APP1"

# ---- Channel / Couette / Jet ----
CHP_MAGIC = b"CHP1"
CUP_MAGIC = b"CUP1"
JCL_MAGIC = b"JCL1"
JSC_MAGIC = b"JSC1"
