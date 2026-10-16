"""
client_text.py -- wording rules for text a client reads: the report, the podcast transcript.

DRAFT HISTORY. A client deliverable states what the analysis found and what was checked, in the
present tense. It never narrates drafts the client never saw. A staff report (PROT_0000, skill
2.11.2): after the step-8d expert review, the report's Expert Review Notes read as a changelog
("checked the first draft", "two protein names were wrong ... corrected"), the audit section said
an earlier normalisation reason was "corrected", and the podcast said "the first draft of this
report fell for it". It confuses a reader who knows no mass spectrometry and makes a sound result
look shaky. The review -> fix history belongs in logs/decisions.md (staff only).

The ONE list of those phrases. make_podcast.py check FAILS on them in the transcript;
make_analysis_html.py WARNS on them in AI_Analysis_Report.md. Each pattern names draft or
revision history, never ordinary science: "p-values were corrected for multiple testing",
"batch-corrected", "fixed modification", "fixed at 15 ppm", "fell for most proteins", "draft
genome", "previously reported interactors", "an earlier version of UniProt" all pass
(tests/test_client_text.py holds both lists). Stdlib only.

  draft_history(text) -> [snippet, ...]   one per place the text narrates a draft
"""
import re

_DOC = (r"(?:report|analysis|write-?up|summary|script|podcast|episode|text|document|results|"
        r"review|section|notes)")
# a draft that is not ours: a genome, the client's manuscript
_NOT_OURS = (r"(?!\s+(?:genome|assembl\w*|sequence\w*|transcriptome|annotation\w*|manuscript|paper|"
             r"grant|thesis|proposal|of\s+(?:your|the|a)\s+(?:manuscript|paper|grant|thesis|"
             r"proposal)))")
# what "corrected" is about in statistics and methods: for multiple testing, using BH, (BH) ...
_STATS_CORRECTION = (r"(?!\s*(?:\(|for\b|using\b|with\b|by\b|via\b|according\b|to\s+account\b|"
                     r"across\b|at\b|against\b|on\b|in\s+the\s+model\b))")
_FIXABLE = (r"(?:it|this|that|them|these|those|the\s+(?:\w+\s+){0,2}?(?:mistakes?|errors?|typos?|"
            r"names?|labels?|legends?|captions?|text|wording|draft|report|sections?|categor\w+|"
            r"reasons?|claims?|statements?|figures?))")

DRAFT_HISTORY = tuple(re.compile(p, re.I) for p in (
    # drafts and versions of THIS report
    r"(?<!your\s)\b(?:first|earlier|previous|prior|original|initial|old|rough)\s+drafts?\b" + _NOT_OURS,
    r"\b(?:the|this|our|my|that)\s+draft\b" + _NOT_OURS,
    r"\bdrafts?\s+of\s+(?:this|the|our|my)\s+" + _DOC + r"\b",
    r"\b(?:earlier|previous|prior|older|first|original|initial)\s+(?:versions?|revisions?|"
    r"iterations?)\s+of\s+(?:this|the|our|my)\s+" + _DOC + r"\b",
    r"\b(?:earlier|previous|prior|older)\s+(?:versions?|revisions?)\b(?!\s+of\b)",
    r"\bthis\s+(?:version|revision)\s+of\s+(?:the|this|our)\s+" + _DOC + r"\b",
    r"\b(?:the|this)\s+" + _DOC + r"\s+now\s+(?:says|reads|states|shows|lists|names|uses|gives|"
    r"correctly)\b",
    r"\bnow\s+correctly\b",
    # corrections and mistakes of our own
    r"\b(?:was|were|has\s+been|have\s+been|had\s+been|is|are)\s+(?:now\s+|since\s+|later\s+|"
    r"then\s+|subsequently\s+)?corrected\b" + _STATS_CORRECTION,
    r"\b(?:now|since|later|subsequently)\s+(?:been\s+)?corrected\b" + _STATS_CORRECTION,
    r"\b(?:we|i)(?:'ve|\s+have|\s+had)?\s+(?:corrected|fixed|amended|revised)\s+" + _FIXABLE + r"\b",
    r"(?<!mass\s)(?<!ppm\s)\b(?:mistakes?|errors?|typos?)\s+(?:was|were|has\s+been|have\s+been)"
    r"\s+(?:now\s+)?(?:corrected|fixed|removed)\b" + _STATS_CORRECTION,
    r"\b(?:was|were)\s+wrong\b[^.!?;]{0,60}\b(?:corrected|fixed|amended)\b",
    r"\b(?:we|i|the\s+" + _DOC + r"|this\s+" + _DOC + r"|our\s+" + _DOC + r")\s+(?:was|were|had\s+"
    r"been|got\s+(?:it|this|that|these|those|them)?)\s*wrong\b",
    r"\bgot\s+(?:it|this|that|these|those|them)\s+wrong\b",
    r"\bhonest\s+correction\b",
    r"\bcorrection\s+(?:to|of)\s+(?:the|an|our|my)\s+(?:earlier|previous|first|original|prior)\b",
    r"\bfell\s+(?:for|into)\s+(?:it|this|that|the\s+trap|that\s+trap|this\s+trap)\b",
    r"\b(?:earlier|previous|original|initial)\s+(?:reason|verdict|wording)\b",
))


def draft_history(text):
    """Each place `text` narrates a draft or a correction of its own: a snippet around it, one
    per sentence. [] when there is none."""
    text = text or ""
    out, seen = [], set()
    for rx in DRAFT_HISTORY:
        for m in rx.finditer(text):
            a = max(text.rfind(c, 0, m.start()) for c in ".!?\n") + 1
            if a in seen:
                continue
            seen.add(a)
            out.append(" ".join(text[max(a, m.start() - 60):m.end() + 40].split()))
    return out
