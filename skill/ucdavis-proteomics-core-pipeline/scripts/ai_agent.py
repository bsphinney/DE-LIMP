#!/usr/bin/env python3
"""
ai_agent.py -- the ONE definition of how the skill names the AI assistant that ran it (2.15.1).

Before 2.15.1 the report, the reproducibility log and the podcast paragraph each said "Claude,
Anthropic's AI model" whatever actually ran the skill, and docs/how-reports-are-made.html named a
model version. A skill can be driven by Claude Code, Claude Desktop or another agent, so the name is
read from what the run records and never assumed:

  1. an `ai_agent` string in the session's session.json (what a session or a staff member recorded);
  2. the environment: Claude Code sets CLAUDECODE=1 (and CLAUDE_CODE_ENTRYPOINT); the model, only
     when ANTHROPIC_MODEL names it;
  3. otherwise "an AI assistant" -- no product and no model.

Stdlib only.
"""
import json
import os

GENERIC = "An AI assistant"


def _session_agent(session):
    if not session:
        return None
    try:
        with open(os.path.join(os.path.expanduser(session), "session.json"),
                  encoding="utf-8") as fh:
            obj = json.load(fh)
    except (OSError, ValueError):
        return None
    name = obj.get("ai_agent") if isinstance(obj, dict) else None
    return name.strip()[:120] if isinstance(name, str) and name.strip() else None


def detected(env=None):
    """The assistant the environment names (step 2 above), or None -- never a guess. What
    session.py init records as session.json `ai_agent`, so a report written later in a job on a
    compute node (no CLAUDECODE there) still names what ran the analysis."""
    env = os.environ if env is None else env
    if env.get("CLAUDECODE") or env.get("CLAUDE_CODE_ENTRYPOINT"):
        model = (env.get("ANTHROPIC_MODEL") or "").strip()
        return (f"Claude Code (Anthropic's AI coding assistant, model {model})" if model
                else "Claude Code (Anthropic's AI coding assistant)")
    return None


def agent_label(session=None, env=None):
    """The assistant's name for a sentence's subject: "Claude Code (Anthropic's AI coding
    assistant)", "Claude Code (model ...)", the session's recorded name, or "An AI assistant"."""
    return _session_agent(session) or detected(env) or GENERIC
