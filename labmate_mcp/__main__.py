"""
labmate-mcp entry point.

Handles:
    labmate-mcp          → Run MCP server (normal operation)
    labmate-mcp --setup  → Interactive API key setup wizard
    labmate-mcp --help   → Show help
"""

import sys


def main():
    """Main entry point for labmate-mcp."""
    
    # Handle --setup flag
    if "--setup" in sys.argv or "-s" in sys.argv:
        from labmate_mcp.setup_wizard import run_setup
        run_setup()
        return
        
    # Handle --help flag
    if "--help" in sys.argv or "-h" in sys.argv:
        print("""
🧪 labmate-mcp — Your AI lab companion

Usage:
    labmate-mcp           Run the MCP server
    labmate-mcp --setup   Configure API keys interactively
    labmate-mcp --help    Show this help message
    
API keys are stored in ~/.labmate-mcp.env and loaded automatically.
Run --setup to add keys for IBM RXN, Rowan Science, and more.

Documentation: https://github.com/JonasRackl/labmate-mcp
""")
        return
        
    # Load env vars from ~/.labmate-mcp.env before starting server
    from labmate_mcp.env_loader import load_env
    load_env()
    
    # Run the MCP server
    from labmate_mcp.server import main as server_main
    server_main()


if __name__ == "__main__":
    main()
