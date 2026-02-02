# Contributing to labmate-mcp

Thanks for considering a contribution! Whether it's a new tool, a bug fix, a new named reaction, or better documentation — it all helps.

## Quick Start

```bash
git clone https://github.com/JonasRackl/labmate-mcp.git
cd labmate-mcp
pip install -e ".[rowan]"
labmate-mcp  # verify it runs
```

## Project Structure

```
labmate_mcp/
├── server.py       Tool definitions (the @mcp.tool() functions)
├── apis.py         External API clients
├── bench.py        Reference databases + calculators
├── writing.py      Writing/publication support
├── peptide.py      Peptide chemistry integrations
├── chemistry.py    Chemistry utilities
└── __init__.py     Version
```

## Adding a New Tool

Every tool follows the same three-part pattern:

### 1. Write the backend function

Add it to the appropriate module (`apis.py` for API calls, `bench.py` for reference data, `writing.py` for writing support, etc.):

```python
# In apis.py or the relevant module
async def my_new_lookup(query: str) -> dict:
    """Fetch data from SomeAPI."""
    async with httpx.AsyncClient() as client:
        resp = await client.get(f"https://api.example.com/{query}")
        resp.raise_for_status()
        return resp.json()
```

### 2. Create the Pydantic input model + tool in `server.py`

```python
class MyNewInput(BaseModel):
    """Description of the tool."""
    query: str = Field(description="What to search for")

@mcp.tool()
async def my_new_tool(params: MyNewInput) -> str:
    """One-line description shown in Claude's tool list. Be specific and useful."""
    result = await my_new_lookup(params.query)
    if not result:
        return "No results found."
    # Format into clean Markdown
    lines = [f"**{result['name']}**", f"Value: {result['value']}"]
    return "\n".join(lines)
```

### 3. Import in `server.py`

Add the import at the top of `server.py` alongside the other imports from that module.

## Adding Named Reactions

Named reactions live in `bench.py` in the `NAMED_REACTIONS` dict. Each entry follows this structure:

```python
_nr("Reaction Name",
    aliases=["alternate name", "abbreviation"],
    category="coupling",  # or reduction, oxidation, rearrangement, etc.
    type_="transition-metal-catalyzed",
    summary="One-sentence description.",
    mechanism="General mechanism description.",
    conditions="Typical reagents, solvents, temperatures.",
    scope="What substrates work well.",
    limitations="What doesn't work.",
    examples=["Example transformation description"],
    references=["Key reference DOI or citation"])
```

## Adding Experimental Templates

Templates live in `writing.py` in the `EXPERIMENTAL_TEMPLATES` dict. Follow the existing format — each template should include placeholder fields in `[BRACKETS]` that the user fills in.

## Code Standards

- **Formatting**: Use clear, readable Python. No strict linter enforced yet.
- **Docstrings**: Every tool needs a clear one-line docstring (this is what Claude sees).
- **Error handling**: Return helpful error messages, never raw stack traces. If an API is unreachable, suggest alternatives.
- **No API keys required** for core tools. If a tool requires a key, gate it with a clear setup message.

## Pull Requests

1. Fork the repo and create a branch (`feature/my-new-tool`)
2. Make your changes
3. Test that the server starts: `labmate-mcp`
4. Update the tool count in `README.md` and `__init__.py` if adding tools
5. Add a CHANGELOG entry
6. Open a PR with a clear description

## Reporting Issues

Open an issue on GitHub. Include:
- What you tried (the question you asked Claude)
- What happened (error message or unexpected result)
- Your Python version and OS

## Questions?

Open an issue or reach out. No question is too small.
