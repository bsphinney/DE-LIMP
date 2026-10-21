"""
client_text.py -- wording rules for text a client reads: the report, the podcast transcript.

DRAFT HISTORY. A client deliverable states what the analysis found and what was checked, in the
present tense. It never narrates drafts the client never saw. A staff report (PROT_0000, skill
2.11.2): after the step-8d expert review, the report's Expert Review Notes read as a changelog
("checked the first draft", "two protein names were wrong ... corrected"), the audit section said
an earlier normalisation reason was "corrected", and the podcast said "the first draft of this
report fell for it". It confuses a reader who knows no mass spectrometry and makes a sound result
look shaky. The review -> fix history belongs in logs/decisions.md (staff only).

The ONE list of those phrases, in two tiers:
  draft_history(text)        STRICT: phrases that only narrate a draft or a fix of this text --
                             "first draft", "an earlier version of this report", "fell for it",
                             "we got this wrong", "the report originally said", "we've updated
                             the wording", "... was incorrect, and it has been changed".
                             make_podcast.py check FAILS on them.
  draft_history_loose(text)  LOOSE: words that are draft history only in context -- "corrected",
                             "fixed", "updated" in a sentence that also speaks of a review, a
                             mistake, a draft or the report; "the draft"; "an earlier pass";
                             "initially reported as". check WARNS on them.
make_analysis_html.py WARNS on both in AI_Analysis_Report.md. "corrected" alone is science ("the
p-values are corrected", "batch effects were corrected"), so it is never a FAIL and never a WARN
without that context. tests/test_client_text.py holds the science sentences that must pass both
tiers, and what each tier must catch. Draft history in words neither tier has is the meaning
check's (f) (make_podcast.MEANING_CHECKS). Stdlib only.

PEAK-CHECK OVERCLAIMS (peak_check_overclaims, at the end): Skyline's peak check is a visual check
of DIA-NN's peaks, never "confirmed by Skyline" or "validated" (skyline_peaks.py).
"""
import re

_DOC = (r"(?:report|analysis|write-?up|summary|script|podcast|episode|text|document|results|"
        r"review|section|notes)")
# a draft that is not ours: a genome, a proteome map, the client's manuscript
_NOT_OURS = (r"(?!\s+(?:genome|proteome|map|assembl\w*|sequence\w*|transcriptome|annotation\w*|"
             r"structure|model|manuscript|paper|grant|thesis|proposal|of\s+(?:your|the|a)\s+"
             r"(?:manuscript|paper|grant|thesis|proposal)))")
# what only a text has (never "the gene names", which a method also corrects)
_OURS = r"(?:wording|labels?|report|text|sections?|draft|legends?|captions?|summary)"
# a fix done FOR the data, not to a text: "corrected for by median normalisation", "with ComBat"
_BY_METHOD = r"(?![^.!?;]{0,40}\b(?:for|using|by|with)\b)"
# a fix to the sample records, not to a text: "corrected in the sample sheet"
_RECORDS = re.compile(r"\s+(?:in|on)\s+the\s+(?:sample\s+sheet|submission(?:\s+form)?|form|"
                      r"manifest|metadata|sheet|tube|vial|label)\b", re.I)
_IN_RECORDS = r"(?!" + _RECORDS.pattern + r")"


def _rx(*patterns):
    return tuple(re.compile(p, re.I) for p in patterns)


DRAFT_HISTORY = _rx(
    # drafts and versions of THIS text
    r"(?<!your\s)\b(?:first|earlier|previous|prior|original|initial|old|rough)\s+drafts?\b" + _NOT_OURS,
    r"\bdrafts?\s+of\s+(?:this|the|our|my)\s+" + _DOC + r"\b",
    r"\b(?:earlier|previous|prior|older|first|original|initial)\s+(?:versions?|revisions?|"
    r"iterations?)\s+of\s+(?:this|the|our|my)\s+" + _DOC + r"\b",
    r"\b(?:an?|the|our|my)\s+(?:earlier|previous|prior|older)\s+version\b(?!\s+of\b)(?!\s*\()"
    r"(?!\s*,?\s*v?\d)",                                    # "the previous version, 2.6.1," is software
    r"\bthis\s+(?:version|revision)\s+of\s+(?:the|this|our)\s+" + _DOC + r"\b",
    r"\b(?:the|this)\s+" + _DOC + r"\s+now\s+(?:says|reads|states|shows|lists|names|uses|gives|"
    r"correctly)\b",
    r"\b(?:the|this|our)\s+" + _DOC + r"\s+(?:originally|initially|first)\s+(?:said|stated|read|"
    r"called|claimed|listed|named|reported)\b",
    r"\boriginally,?\s+(?:the|this|our)\s+" + _DOC + r"\s+(?:said|stated|read|called|claimed)\b",
    r"\b(?:our|the|an?)\s+(?:first|earlier|previous)\s+pass\s+(?:of|at|on)\s+(?:this|the|our)\s+"
    + _DOC + r"\b",
    r"\b(?:earlier|previous|original|initial)\s+(?:\w+\s+)?(?:reason|verdict|wording)\b",
    # a fix of our own, told as one
    r"\b(?:wrong|incorrect|mislabell?ed)\b[^.!?;]{0,80}\b(?:has|have)\s+been\s+(?:changed|"
    r"corrected|fixed|updated|amended)\b" + _IN_RECORDS,
    r"\b(?:caught|flagged|spotted)\b[^.!?;]{0,80}\b(?:which|that)?\s*(?:we|then)\s+(?:then\s+)?"
    r"(?:fixed|changed|updated|corrected|amended)\b" + _BY_METHOD,
    r"\b(?:we|i)(?:[’']ve|\s+have)?\s+(?:updated|changed|revised|amended|reworded|fixed|corrected)"
    r"\s+the\s+(?:\w+\s+){0,2}?" + _OURS + r"\b" + _BY_METHOD,
    r"\b(?:they|it)\s+(?:are|is)\s+(?:now\s+)?fixed\s+now\b|\b(?:is|are)\s+fixed\s+now\b",
    r"\b(?:we|i)\s+got\s+(?:it|this|that|these|those|them)\s+wrong\b",
    r"\bhonest\s+correction\b",
    r"\bcorrection\s+(?:to|of)\s+(?:the|an|our|my)\s+(?:earlier|previous|first|original|prior)\b",
    r"\bfell\s+for\s+it\b|\bfell\s+into\s+(?:the|that|this)\s+trap\b",
)
# Draft history only in context: a sentence that also speaks of a review, a mistake, a draft or
# the report. "mass error" and "error rate" are science, not a mistake.
_CONTEXT = re.compile(r"\b(?:review(?:er|ed)?|draft|earlier|originally|initially|mistakes?|"
                      r"(?<!mass\s)(?<!ppm\s)errors?(?!\s+rates?)|wrong(?:ly)?|incorrect(?:ly)?|"
                      r"mislabell?ed|typos?|trap|backwards|(?:the|this)\s+report)\b", re.I)
_FIX_VERB = re.compile(r"\b(?:corrected|fixed|updated|changed|amended|revised)\b", re.I)
DRAFT_HISTORY_LOOSE = _rx(
    r"\b(?:the|this|our|my|that)\s+draft\b" + _NOT_OURS,
    r"\b(?:an?|the|our)\s+(?:earlier|previous)\s+pass\b",
    r"\b(?:initially|originally)\s+(?:reported|called|listed|named|described)\s+as\b",
    r"\bhad\s+(?:it|this|that)\s+backwards\b",
    r"\bfell\s+(?:for|into)\s+it\b",
    r"\bnow\s+correctly\b",
)


def _sentence_start(text, at):
    return max(text.rfind(c, 0, at) for c in ".!?\n") + 1


def _hits(text, patterns, seen):
    out = []
    for rx in patterns:
        for m in rx.finditer(text):
            a = _sentence_start(text, m.start())
            if a in seen:
                continue
            seen.add(a)
            out.append(" ".join(text[max(a, m.start() - 60):m.end() + 40].split()))
    return out


def draft_history(text):
    """STRICT: each place `text` narrates a draft or a fix of its own -- a snippet, one per
    sentence. [] when there is none."""
    return _hits(text or "", DRAFT_HISTORY, set())


def draft_history_loose(text):
    """LOOSE: places that may narrate a draft, by context, and that draft_history did not
    already name -- a snippet, one per sentence."""
    text = text or ""
    seen = set()
    _hits(text, DRAFT_HISTORY, seen)                  # a sentence the strict tier names: once
    out = _hits(text, DRAFT_HISTORY_LOOSE, seen)
    for m in _FIX_VERB.finditer(text):
        a = _sentence_start(text, m.start())
        ends = [i for i in (text.find(c, m.end()) for c in ".!?\n") if i >= 0]
        if (a in seen or not _CONTEXT.search(text[a:min(ends or [len(text)])])
                or _RECORDS.match(text, m.end())):
            continue
        seen.add(a)
        out.append(" ".join(text[max(a, m.start() - 60):m.end() + 40].split()))
    return out


# PEAK-CHECK OVERCLAIMS. Skyline's peak check (skyline_peaks.py) draws DIA-NN's peaks again from
# the raw files, but Skyline integrates inside DIA-NN's peak edges by default, so a peak that
# matches is not independent confirmation. A client text never says the results were confirmed,
# validated, verified or corroborated by Skyline (or by "the peak check"). A sentence that SAYS it
# is not ("... is not independent confirmation") is fine: a negation earlier in the same clause
# clears it. Two tiers, as for draft history:
#   peak_check_overclaims        STRICT: the claim is about the peak check -- its verb in the
#                                same clause as Skyline / the peak check, or in a clause that
#                                names no other actor ("Skyline agreed, which confirms ...").
#                                make_podcast.py check FAILS on them.
#   peak_check_overclaims_loose  LOOSE: the verb belongs to another check or actor in the
#                                sentence ("..., and the contaminant check confirmed that ...";
#                                "..., we confirmed ..."), or the reader is asked to do it ("open
#                                the Skyline file to verify the peaks yourself"). check WARNS
#                                (review of 2.14, LOW-2).
# skyline_peaks.report_section() is tested to have neither; make_analysis_html.py warns on both
# in AI_Analysis_Report.md.
_PEAK_SUBJECT = re.compile(r"\b(?:skyline|peak[\s-]+checks?|second\s+program)\b", re.I)
_OVERCLAIM = re.compile(r"\b(?:confirm\w*|validat\w*|verif\w*|corroborat\w*|prov(?:e|es|ed|en)\b|"
                        r"proof|independent(?:ly)?\s+(?:check\w*|evidence|measure\w*))", re.I)
_NEGATION = re.compile(r"\b(?:not|no|never|nor|without|cannot|isn['’]t|aren['’]t|doesn['’]t|"
                       r"don['’]t|wasn['’]t|weren['’]t|neither)\b|n['’]t\b", re.I)
_SENTENCE = re.compile(r"[^.!?\n]+[.!?]?")
_CLAUSE_BREAKS = ",;:—("
# another check, method or person that a clause's verb can belong to
_OTHER_ACTOR = re.compile(r"\b(?!peak\b)[\w-]+\s+(?:checks?|tests?|audit|analys[ie]s|blots?|"
                          r"assays?)\b|\b(?:we|they|DIA-NN|the\s+(?:audit|DE|data|statistics|"
                          r"review(?:ers?)?))\b", re.I)
# the reader is invited to do it: not a claim of what was done
_INVITATION = re.compile(r"\byourself\b|\byou\s+(?:can|could|may|might|will)\b", re.I)


def _peak_check_tiers(text):
    """-> (strict, loose): one snippet per sentence, in its stricter tier."""
    strict, loose = [], []
    for m in _SENTENCE.finditer(text or ""):
        sent = m.group(0)
        if not _PEAK_SUBJECT.search(sent):
            continue
        tier = None
        for c in _OVERCLAIM.finditer(sent):
            start = max(sent.rfind(ch, 0, c.start()) for ch in _CLAUSE_BREAKS) + 1
            if _NEGATION.search(sent[start:c.start()]):
                continue
            ends = [i for i in (sent.find(ch, c.end()) for ch in _CLAUSE_BREAKS) if i >= 0]
            clause = sent[start:min(ends) if ends else len(sent)]
            mine = _PEAK_SUBJECT.search(clause) or not _OTHER_ACTOR.search(clause)
            if mine and not _INVITATION.search(clause):
                tier = "strict"
                break
            tier = "loose"
        if tier:
            (strict if tier == "strict" else loose).append(" ".join(sent.split()))
    return strict, loose


def peak_check_overclaims(text):
    """STRICT: each sentence that claims Skyline (or the peak check) confirmed, validated,
    verified, corroborated or proved a result -- a snippet per sentence; [] when there is none."""
    return _peak_check_tiers(text)[0]


def peak_check_overclaims_loose(text):
    """LOOSE: sentences with such a verb near Skyline or the peak check that belongs to another
    check or actor, or that asks the reader to do it -- for a person to read, never a failure."""
    return _peak_check_tiers(text)[1]
