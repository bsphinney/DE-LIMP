# Recording skill problems (`report_issue.sh`)

The skill is fixed from these reports (golden rule 8). SKILL.md keeps the command and the
quoting rule; this file holds the rest: where a report lands, how each argument is passed, and
what to record and what not.

## When to record (golden rule 8)

8. **Record every skill problem the moment it happens — `scripts/report_issue.sh`.**
   A script that errors or answers wrongly, a documented step that did not match reality, a
   workaround you had to invent, something you had to find or install that the skill should
   have, the user correcting how you used the skill: record it **then**, not at the end —
   sessions end without warning, and an unrecorded problem is one nobody fixes. Tell the user
   in one line when you do. → "The command, and what to record" below.

## The command, and what to record

The skill is fixed from these reports. For a Core member they land in the Core's shared
folder on HIVE, `/quobyte/proteomics-grp/skill_issues/`, one file per entry, where the
maintainers read them. A report holds only what you write in it. The Slack post at the end of
the search or analysis counts them (`references/notifications.md`), so the Core sees them
without anyone forwarding a file.

```
bash scripts/report_issue.sh --title '<short name>' --what-file - \
    --impact '<time, compute or storage lost, or what the result would have got wrong>' \
    [--workaround '<what you did instead>'] [--fix '<what the skill should do>'] \
    [--kind bug|docs|missing|agent_mistake] [--severity high|medium|low] \
    [--step '<SKILL.md step>'] --session '<session name>' <<'EOF'
<what happened: the exact command, the exact error text, the path -- verbatim>
EOF
```

**The exact text goes in on standard input** (`--what-file -` and a heredoc quoted `<<'EOF'`,
in which nothing is expanded), **never inside `--what "..."`**: in double quotes bash runs
anything in backticks or `$(...)` — an error text full of them — on this computer, and again
on HIVE when the command goes through `hive_exec.sh`. Every other argument in single quotes;
one that holds a single quote, or a second long text, goes in a file (`--impact-file F`,
`--workaround-file F`, `--fix-file F`, `--title-file F`; only one field reads standard input).
It runs **locally** and needs only bash (Git Bash on Windows, no local Python): in
`hive_remote` it writes to HIVE through `hive_exec.sh` by itself; on HIVE it writes directly;
anywhere else — or if the account cannot write the Core folder — it keeps the file in
`~/.proteomics-pipeline/issues/` and prints how to send it to the Core. It never fails the
session. Use one `--session` for the whole session (the step-3b name, or a short tag before
there is one).

**Record:** a skill script that errors, hangs, or returns something wrong or misleading; a
step in SKILL.md or a reference that did not match what you found; a tool, path or setting
the skill should have found and you had to find by hand; a workaround you invented; anything
that cost the user significant time, cluster time or storage; the user correcting the skill's
behaviour. **Also record your own mistakes** that cost the user something (`--kind
agent_mistake`) — telling the user setup installed a tool it had not, searching the wrong
filesystem before asking where a drive points. **Do not record** problems with the user's
data (a failed injection, a bad run — those belong in the analysis report) or a network blip
that a retry fixed.

Write it so someone who was not there can fix it: the exact command, the error text verbatim,
the paths, the numbers, and what you expected instead. **Never** a password, private key or
token (the script refuses text that looks like one).
