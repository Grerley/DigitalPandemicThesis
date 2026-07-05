#!/bin/bash
# SessionStart hook — prepare the Digital Pandemic Simulation Lab for web sessions.
# The app is pure client-side ES modules with NO runtime dependencies, so there
# is nothing to install. This hook instead (1) runs the headless acceptance
# self-test so correctness is verified before the session begins, and (2) serves
# the app on a local static server so it is immediately reachable in the browser.
# Idempotent and non-interactive.
set -euo pipefail

# Only act in the Claude Code on the web (remote) environment.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

PROJECT_DIR="${CLAUDE_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
APP_DIR="$PROJECT_DIR/app"
PORT="${DPL_PORT:-8000}"

echo "[digital-pandemic-lab] preparing simulation lab…"

# 1) Headless acceptance self-test (asserts the thesis targets under Node).
if command -v node >/dev/null 2>&1; then
  if node "$APP_DIR/test/selftest.mjs"; then
    echo "[digital-pandemic-lab] acceptance self-test: PASSED"
  else
    echo "[digital-pandemic-lab] acceptance self-test: FAILED (see output above)" >&2
  fi
else
  echo "[digital-pandemic-lab] node not found — skipping self-test" >&2
fi

# 2) Serve the app on a static server (background), if the port is free.
if command -v python3 >/dev/null 2>&1; then
  if ! curl -s -o /dev/null "http://localhost:${PORT}/" 2>/dev/null; then
    nohup python3 -m http.server "$PORT" --directory "$APP_DIR" \
      >/tmp/dpl-server.log 2>&1 &
    disown || true
    sleep 1
    echo "[digital-pandemic-lab] serving app at http://localhost:${PORT}/ (logs: /tmp/dpl-server.log)"
  else
    echo "[digital-pandemic-lab] a server is already listening on port ${PORT}"
  fi
else
  echo "[digital-pandemic-lab] python3 not found — start a static server manually (see app/README.md)" >&2
fi

# Expose the app dir + URL to the rest of the session.
if [ -n "${CLAUDE_ENV_FILE:-}" ]; then
  {
    echo "export DPL_APP_DIR=\"$APP_DIR\""
    echo "export DPL_URL=\"http://localhost:${PORT}/\""
  } >> "$CLAUDE_ENV_FILE"
fi

echo "[digital-pandemic-lab] ready."
