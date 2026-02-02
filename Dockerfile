FROM python:3.12-slim

LABEL maintainer="Jonas Rackl"
LABEL description="labmate-mcp — 78-tool MCP server for chemistry research"
LABEL version="7.0.0"

WORKDIR /app

# Install dependencies first (cacheable layer)
COPY pyproject.toml README.md ./
RUN pip install --no-cache-dir .

# Copy source code
COPY labmate_mcp/ labmate_mcp/

# Re-install to register the package properly
RUN pip install --no-cache-dir .

ENV PYTHONUNBUFFERED=1

ENTRYPOINT ["labmate-mcp"]
