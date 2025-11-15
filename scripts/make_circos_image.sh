#!/bin/bash
#--------------------------------------------------
# Script to make circos image.
#
# Usage:
#   $ make_circos_image.sh WORKING_DIR
#
# Environment variables (optional):
#   CIRCOS_DIR  - Path to Circos installation (default: auto-detect from PATH)
#   PERL        - Perl executable (default: perl)
#   BASENAME    - Base name for output (default: image)
#--------------------------------------------------
export LC_ALL=C

# Parse command line arguments
WORKING_DIR=$1

if [ -z "$WORKING_DIR" ]; then
    echo "Error: WORKING_DIR argument is required"
    echo "Usage: $0 WORKING_DIR"
    exit 1
fi

# Set defaults from environment or use fallback values
PERL=${PERL:-perl}
BASENAME=${BASENAME:-image}

# Auto-detect CIRCOS_DIR if not set
if [ -z "$CIRCOS_DIR" ]; then
    # Try to find circos in PATH
    CIRCOS_BIN=$(which circos 2>/dev/null)
    if [ -n "$CIRCOS_BIN" ]; then
        CIRCOS_DIR=$(dirname $(dirname "$CIRCOS_BIN"))
    else
        echo "Error: CIRCOS_DIR not set and circos not found in PATH"
        echo "Please set CIRCOS_DIR environment variable or install circos in PATH"
        exit 1
    fi
fi

# Verify CIRCOS_DIR exists
if [ ! -d "$CIRCOS_DIR" ]; then
    echo "Error: CIRCOS_DIR=$CIRCOS_DIR does not exist"
    exit 1
fi

# Verify circos binary exists
if [ ! -f "$CIRCOS_DIR/bin/circos" ]; then
    echo "Error: circos binary not found at $CIRCOS_DIR/bin/circos"
    exit 1
fi

cd "$WORKING_DIR" || exit 1
$PERL "$CIRCOS_DIR/bin/circos" -param random_string="$BASENAME" -conf ./etc/circos-medium.conf
# Uncomment for cluster-based circos:
# $PERL "$CIRCOS_DIR/bin/circos" -param random_string="$BASENAME" -conf ./etc/circos-clusterbased-medium.conf





