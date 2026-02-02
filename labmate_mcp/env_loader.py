"""
Environment variable loader for labmate-mcp.

Automatically loads API keys from ~/.labmate-mcp.env on import.
Call load_env() at the top of server.py before any API clients are initialized.
"""

import os
from pathlib import Path

ENV_FILE = Path.home() / ".labmate-mcp.env"


def load_env():
    """
    Load environment variables from ~/.labmate-mcp.env if it exists.
    
    This function should be called at the very start of server.py,
    before any API clients are initialized.
    
    The .env file format is simple:
        KEY=value
        # comments are ignored
        
    Existing environment variables are NOT overwritten.
    """
    if not ENV_FILE.exists():
        return
    
    loaded = []
    with open(ENV_FILE) as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith("#"):
                continue
            
            # Parse KEY=value
            if "=" not in line:
                continue
            
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()
            
            # Don't overwrite existing env vars
            if key not in os.environ and value:
                os.environ[key] = value
                loaded.append(key)
    
    return loaded


# Auto-load on import
_loaded = load_env()
