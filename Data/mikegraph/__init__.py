"""
MikeGraph - Python toolkit for hydraulic network analysis with MIKE+ databases

A comprehensive toolkit for hydraulic engineers and water management professionals
working with MIKE+ databases. Provides tools for network analysis, catchment modeling,
pipe network extraction, and hydraulic calculations.

Key Features:
    - Network graph creation from MIKE+ databases
    - Time-area catchment analysis
    - Pipe network extraction and analysis  
    - Full-flow capacity calculations using Colebrook-White equation
    - Integration with ArcGIS Pro toolboxes
    - Support for both .mdb and .sqlite database formats
"""

# Import main classes
from .graph import MikeNetwork
from .timearea import TimeAreaAnalyzer
from .network import PipeNetwork
from .utils import calculate_full_flow

# Package metadata
__version__ = "2.0.0"
__author__ = "Emil Nielsen"
__email__ = "enielsen93@hotmail.com"
__description__ = "Python toolkit for hydraulic network analysis with MIKE+ databases"
__url__ = "https://github.com/enielsen93/mikegraph"

# Define what gets imported with "from mikegraph import *"
__all__ = [
    "Graph",
    "TimeAreaAnalyzer",
    "PipeNetwork",
    "calculate_full_flow"
]

# Version info tuple (major, minor, patch)
__version_info__ = tuple(map(int, __version__.split('.')))

def get_version():
    """Return the version string."""
    return __version__

# Optional: Check for required dependencies on import
def _check_dependencies():
    """Check if critical dependencies are available."""
    missing_deps = []
    
    try:
        import pandas
    except ImportError:
        missing_deps.append("pandas")
        
    try:
        import numpy
    except ImportError:
        missing_deps.append("numpy")
        
    # Optional dependencies
    optional_missing = []
    try:
        import networkx
    except ImportError:
        optional_missing.append("networkx (required for network analysis)")
        
    try:
        import mikeio
    except ImportError:
        optional_missing.append("mikeio (required for DFS0 files)")
    
    if missing_deps:
        raise ImportError(
            f"Missing required dependencies: {', '.join(missing_deps)}\n"
            f"Install with: pip install {' '.join(missing_deps)}"
        )
    
    if optional_missing:
        import warnings
        warnings.warn(
            f"Optional dependencies not found: {', '.join(optional_missing)}",
            ImportWarning
        )

# Run dependency check on import
_check_dependencies()