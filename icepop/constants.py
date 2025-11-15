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

# ============================================================================
# Circos Visualization Constants
# ============================================================================

# Default circle radius for Circos plots
DEFAULT_CIRCOS_RADIUS = 10000

# Cell type RGB color mapping for Circos visualizations
CELLTYPE_RGB_COLORS = {
    "Bcells": "(23,190,207)",
    "DendriticCells": "(188,189,34)",
    "Macrophages": "(127,127,127)",
    "Macrophage": "(127,127,127)",
    "Monocytes": "(207,236,249)",
    "Monocyte": "(207,236,249)",
    "NKCells": "(140,86,75)",
    "Neutrophil": "(148,103,189)",
    "Neutrophils": "(148,103,189)",
    "StemCells": "(214,39,40)",
    "TCells": "(214,39,40)",
    "StromalCells": "(44,160,44)",
    "abTcells": "(255,127,14)",
    "gdTCells": "(31,119,180)"
}

# Cell type abbreviations for Circos labels
CELLTYPE_ABBREVIATIONS = {
    "Bcells": "B",
    "DendriticCells": "DC",
    "Macrophages": "Mac",
    "Macrophage": "Mac",
    "Monocytes": "Mo",
    "Monocyte": "Mo",
    "NKCells": "NK",
    "Neutrophil": "Neu",
    "Neutrophils": "Neu",
    "StemCells": "Stem",
    "StromalCells": "Stro",
    "TCells": "T",
    "abTcells": "abT",
    "gdTCells": "gdT"
}

# ============================================================================
# Fold Change Analysis Constants
# ============================================================================

# Default fold change range for enumeration
DEFAULT_FOLD_CHANGE_RANGE = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
