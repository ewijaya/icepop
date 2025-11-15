#!/usr/bin/env python
"""
Constants used across the icepop package.
Centralizes magic numbers and configuration values for better maintainability.
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"

# ============================================================================
# Cell Population Scoring Constants
# ============================================================================

# Decimal precision for rounding cell population scores
CELLPOP_SCORE_DECIMALS = 5

# Decimal precision for display/output
CELLPOP_DISPLAY_DECIMALS = 3

# Replacement value for infinity in calculations
INFINITY_REPLACEMENT = 999

# Default threshold calculation method
DEFAULT_THRESHOLD_METHOD = 'median'

# Available threshold methods
THRESHOLD_METHODS = ['median', 'q1', 'q1std', 'uniform']

# Percentile for quartile-based threshold
FIRST_QUARTILE_PERCENTILE = 25

# ============================================================================
# Deconvolution Constants
# ============================================================================

# Tikhonov regularization lambda parameter
TIKHONOV_LAMBDA = 0.5

# Default normalization method for deconvolution
DEFAULT_NORMALIZATION_METHOD = 'v3'

# Available normalization methods
NORMALIZATION_METHODS = ['v1', 'v2', 'v3', 'v4']

# ============================================================================
# Specificity Score Constants
# ============================================================================

# Default specificity score limit
DEFAULT_SPECIFICITY_LIMIT = 0.8

# Default top-k genes to select per cell type
DEFAULT_TOP_K_GENES = 1

# Number of threshold steps for condition number enumeration
SPECIFICITY_THRESHOLD_STEPS = 10

# ============================================================================
# Correlation Constants
# ============================================================================

# Default number of top correlated phenotypes to select (when all negative)
DEFAULT_TOP_PHENOTYPES = 10

# Default number of top phenotypes (version 3)
DEFAULT_TOP_PHENOTYPES_V3 = 100

# ============================================================================
# Enumerate Output Constants
# ============================================================================

# Timeout for network requests (in seconds)
NETWORK_REQUEST_TIMEOUT = 30

# ============================================================================
# General Constants
# ============================================================================

# CSV file opening mode (use 'r' for Python 3 instead of 'rb')
CSV_READ_MODE = 'r'

# Default number of rows to skip when reading CSV
DEFAULT_SKIPROWS = [2]

# Column names for standard DEG input
DEG_COLUMN_NAMES = ["probe", "Genes"]
