# post/io/formats.py
from __future__ import annotations

import numpy as np

# ---- TKE ----
TKE_DTYPE = np.dtype([("t_star", "<i8"), ("ke", "<f8")])  # little-endian

# ---- Centerline ----
CLN_MAGIC = b"CLN1"

# ---- Annular profile ----
APR_MAGIC = b"APR1"

APP_MAGIC = b"APP1"
