# ==============================================================================
#  HELPER FUNCTIONS — AI Integration (Google Gemini + OpenAI-compatible)
# ==============================================================================
#
#  Two providers are supported:
#    gemini        — Google's generativeLanguage API (File API for data upload)
#    openai_compat — any OpenAI-compatible /v1/chat/completions endpoint
#                    (vLLM, llama.cpp, LiteLLM, the UC Davis metabolomics
#                    gateway, ...). No File API: the data table is inlined.
#
#  Call sites MUST use the ask_ai_* dispatchers, never ask_gemini_* directly,
#  so that adding a third provider requires no edits outside this file
#  (CLAUDE.md architectural rule #3: concepts have one definition).
#
#  SECURITY (v4.1.0 review). These requests run inside the single R process that
#  serves every Shiny session, and on the public Hugging Face Space the person
#  typing the endpoint URL and API key is an anonymous visitor. Three rules follow
#  and are enforced here rather than at call sites:
#    1. A typed endpoint URL is validated before any request (ai_validate_endpoint):
#       no credentials/query in the URL; on a public deployment https only and the
#       host must resolve to public addresses only; on a local install plain http
#       is also accepted to this computer / the local network / the container host
#       (with a warning). Link-local and metadata addresses are refused everywhere.
#       The validated addresses are pinned for the request and redirects are
#       refused, so DNS cannot be swapped between the check and the connection.
#    2. Upstream error bodies are never echoed. A failure becomes one short,
#       sanitised status line (ai_error_from_condition).
#    3. Every request has a timeout, capped on a public deployment because a
#       blocked R process blocks every visitor (ai_request_timeout).
# ==============================================================================

# --- PROVIDER REGISTRY (the single definition) ---
ai_providers <- function() {
  list(
    gemini = list(
      label             = "Google Gemini (cloud)",
      destination       = "Google's Gemini API",
      key_label         = "Gemini API Key",
      key_placeholder   = "AIzaSy...",
      # Google's documented key format. Used ONLY to refuse sending a key that is
      # recognisably a Gemini key to any other provider's endpoint.
      key_pattern       = "^AIza[0-9A-Za-z_-]{35}$",
      # Hosts a key of this provider may legitimately be sent to under ANOTHER
      # provider setting, e.g. Google's own OpenAI-compatible endpoint
      # https://generativelanguage.googleapis.com/v1beta/openai
      key_hosts         = c("generativelanguage.googleapis.com"),
      default_model     = "gemini-2.5-flash",
      base_url          = "https://generativelanguage.googleapis.com/v1beta",
      supports_file_api = TRUE,
      # Google responds quickly; no need for a long ceiling
      timeout_s         = 300,
      # ~1M token window: the full 800-protein table fits comfortably
      max_payload_chars = 400000
    ),
    openai_compat = list(
      label             = "OpenAI-compatible endpoint",
      destination       = "the OpenAI-compatible endpoint you configured",
      key_label         = "API Key",
      key_placeholder   = "sk-...",
      key_pattern       = NULL,   # gateway keys have no fixed format
      key_hosts         = NULL,
      # Measured on the UC Davis gateway (2026-09-10): flash-next answered a
      # trivial prompt in 3s vs 35s / 67s / >90s for the other three, and
      # produced AI Summary output of equal quality. Fastest wins here because
      # every one of these is a thinking model and latency is queue-dominated.
      default_model     = "qwen3.8-flash-next",
      base_url          = "https://llm.metabolomics.us/v1",
      supports_file_api = FALSE,
      # Measured 2026-09-10: qwen3.8-27b exceeded 300s on a Data Chat request
      # and 90s on a two-token prompt. Latency here is queue-dominated, so the
      # ceiling is deliberately generous while the gateway is being tuned.
      # NOTE: this is the LOCAL ceiling. On a public deployment every request is
      # capped at AI_PUBLIC_MAX_TIMEOUT_S — see ai_deployment_policy().
      timeout_s         = 900,
      # The smallest gateway models expose 65,536 tokens. DE tables tokenise at
      # ~1.01 chars/token (measured — prose is ~6.0), so chars ~= tokens here.
      # 40,000 leaves ~25k for the QC table, system prompt, reasoning and answer.
      max_payload_chars = 40000
    )
  )
}

# Measured chars-per-token for the numeric DE tables this app sends.
# Calibrated against the gateway's own prompt_tokens accounting on 2026-09-10:
# 50/200/800-protein tables gave 1.00 / 1.01 / 1.01. Do NOT substitute the
# usual ~4 chars/token rule of thumb — that is a PROSE figure (measured 6.0
# here) and under-counts a numeric table by roughly 6x.
ai_chars_per_token <- function() 1.01

# Read one field for a provider, falling back to gemini for unknown/NULL input
# so a stale saved session can never blank out a user-facing label.
ai_provider_field <- function(provider, field) {
  p <- ai_providers()
  key <- if (is.null(provider) || length(provider) != 1 || is.na(provider) ||
             !nzchar(provider) || !(provider %in% names(p))) "gemini" else provider
  p[[key]][[field]]
}

# Human-readable record of WHO produced a piece of AI text, for exports and the
# chat transcript (CLAUDE.md rule #1: an export must say what actually ran, not
# a hardcoded "Gemini"). The endpoint host is included for openai_compat
# because "an OpenAI-compatible endpoint" alone does not identify the service.
ai_source_label <- function(provider, model = NULL, base_url = NULL) {
  prov <- if (is.null(provider) || length(provider) != 1 || is.na(provider) ||
              !(provider %in% names(ai_providers()))) "gemini" else provider
  label <- ai_provider_field(prov, "label")
  if (identical(prov, "openai_compat")) {
    u <- if (is.null(base_url) || length(base_url) != 1 || is.na(base_url) ||
             !nzchar(trimws(base_url))) ai_provider_field(prov, "base_url") else trimws(base_url)
    host <- tryCatch(httr2::url_parse(u)$hostname, error = function(e) NULL)
    if (length(host) == 1 && !is.na(host) && nzchar(host)) label <- paste0(label, " at ", host)
  }
  m <- if (is.null(model) || length(model) != 1 || is.na(model) || !nzchar(trimws(model)))
    "model not recorded" else paste0("model ", trimws(model))
  paste0(label, ", ", m)
}

# ==============================================================================
#  DEPLOYMENT POLICY — the single definition of what a deployment may do
# ==============================================================================

# Longest a single AI request may block the R process on a PUBLIC deployment.
# Every Shiny session shares one single-threaded R process, so a request that
# waits N seconds freezes the app for every visitor for N seconds. The previous
# ceiling (900 s) meant one slow gateway model could take the public Space down
# for 15 minutes. 180 s is the trade-off: the default gateway model answered a
# trivial prompt in 3 s (measured 2026-09-10) and Gemini Flash is typically far
# faster than this, while the one model measured beyond 300 s (qwen3.8-27b)
# could not have been served by the old Gemini-sized 300 s ceiling either.
# Local installs keep the provider's own ceiling — there the only session a slow
# request blocks is the user's own. Asynchronous requests would remove the
# trade-off, but `future` is not in the Docker base image and new packages need
# a base-image rebuild (CLAUDE.md), so that is deferred.
AI_PUBLIC_MAX_TIMEOUT_S <- 180

# Names a container uses to reach the machine it runs on. In Docker mode (the
# recommended Windows setup) "localhost" is the container itself, so a model
# server on the user's PC is reached as host.docker.internal.
AI_CONTAINER_HOST_ALIASES <- c("host.docker.internal", "gateway.docker.internal",
                               "host.containers.internal", "host.lima.internal")

# `public = TRUE` is the SAFE default: a caller that forgets to say what kind of
# deployment it is gets the strict policy, never the permissive one.
ai_deployment_policy <- function(public = TRUE) {
  public <- !isFALSE(public)
  list(
    public                  = public,
    # Loopback / RFC1918 / CGNAT / ULA endpoints (a vLLM server on this machine
    # or the lab LAN) are legitimate for a local install and are exactly the
    # SSRF targets on a public one. Link-local and cloud-metadata ranges are
    # refused everywhere.
    allow_private_endpoints = !public,
    # Plain http — the default for vLLM, llama.cpp and Ollama — only to a server
    # on this computer, the local network, or the container host, and only on a
    # local install (the request carries the API key unencrypted, so it is
    # warned about). A public deployment is https-only.
    allow_http_private      = !public,
    max_timeout_s           = if (public) AI_PUBLIC_MAX_TIMEOUT_S else Inf,
    # Activity-log notes are attributed by output folder + OS user name. On a
    # public deployment every visitor runs as the same OS user in the same HOME,
    # so no row can be attributed to the visitor in front of the screen: send none.
    send_activity_notes     = !public
  )
}

# Is this a PUBLIC deployment? DELIMP_PUBLIC_DEPLOYMENT overrides detection:
# 1/true/yes forces the public policy, 0/false/no forces the local policy. Unset
# or empty falls back to `detected` (app.R passes is_hf_space). Any other value
# is treated as public — a typo must fail closed, not open.
ai_resolve_public_deployment <- function(env_value = Sys.getenv("DELIMP_PUBLIC_DEPLOYMENT", ""),
                                         detected = FALSE) {
  v <- tolower(trimws(as.character(env_value %||% "")))
  if (length(v) != 1 || is.na(v) || !nzchar(v)) return(isTRUE(detected))
  if (v %in% c("1", "true", "yes", "on")) return(TRUE)
  if (v %in% c("0", "false", "no", "off")) return(FALSE)
  message("[DE-LIMP] DELIMP_PUBLIC_DEPLOYMENT='", env_value,
          "' not recognised (use 1 or 0); applying the public AI policy")
  TRUE
}

# The timeout actually applied to a request: the explicit override or the
# provider's registry value, never above the deployment cap.
ai_request_timeout <- function(provider, override = NULL, policy = ai_deployment_policy()) {
  t <- suppressWarnings(as.numeric(override %||% ai_provider_field(provider, "timeout_s")))
  if (length(t) != 1 || is.na(t) || t <= 0) t <- ai_provider_field(provider, "timeout_s")
  min(t, policy$max_timeout_s)
}

# Which other provider a key recognisably belongs to (NULL if none). Guards the
# credential leak where a Gemini key typed earlier would otherwise go out as a
# Bearer token to a different provider's endpoint. `host` is the destination:
# a key sent to one of its own provider's hosts (Google's OpenAI-compatible
# endpoint for a Gemini key) is not a leak.
ai_key_owner <- function(provider, api_key, host = NULL) {
  if (is.null(api_key) || length(api_key) != 1 || is.na(api_key) || !nzchar(api_key)) return(NULL)
  h <- if (is.null(host) || length(host) != 1 || is.na(host)) "" else tolower(host)
  for (o in setdiff(names(ai_providers()), provider)) {
    pat <- ai_provider_field(o, "key_pattern")
    if (!is.null(pat) && grepl(pat, trimws(api_key)) &&
        !(nzchar(h) && h %in% tolower(ai_provider_field(o, "key_hosts") %||% character(0))))
      return(o)
  }
  NULL
}

ai_key_belongs_to_other_provider <- function(provider, api_key, host = NULL) {
  !is.null(ai_key_owner(provider, api_key, host))
}

ai_foreign_key_message <- function(owner, host) {
  owner_label <- ai_provider_field(owner, "label")
  paste0("the API key looks like a ", owner_label, " key, so it was not sent to '", host,
         "': a key sent to a server that does not belong to its provider can be read and reused ",
         "by that server. To continue, enter the key issued for '", host, "' (a proxy such as ",
         "LiteLLM normally has its own key), or choose ", owner_label, " as the AI provider",
         if (identical(owner, "gemini"))
           ", or use Google's OpenAI-compatible endpoint https://generativelanguage.googleapis.com/v1beta/openai" else "",
         ".")
}

# The host a request for `provider` would actually go to: the provider's fixed
# host for Gemini, the typed endpoint (or its registry default when empty) for
# openai_compat. No DNS. Used to bind a typed API key to where it may be sent.
ai_effective_host <- function(provider, base_url = NULL) {
  prov <- if (is.null(provider) || length(provider) != 1 || is.na(provider) ||
              !(provider %in% names(ai_providers()))) "gemini" else provider
  u <- if (identical(prov, "openai_compat") && !is.null(base_url) && length(base_url) == 1 &&
           !is.na(base_url) && nzchar(trimws(base_url))) trimws(base_url) else ai_provider_field(prov, "base_url")
  h <- tryCatch(httr2::url_parse(u)$hostname, error = function(e) NULL)
  if (length(h) != 1 || is.na(h)) "" else tolower(gsub("^\\[|\\]$", "", h))
}

# A typed key is bound to the provider + effective host current when it was
# entered. A request may use the key only while both are unchanged; this is the
# SERVER-side guard (the browser-side clearing can lose a race with a click).
ai_key_binding <- function(provider, base_url, api_key) {
  list(provider = provider %||% "gemini", host = ai_effective_host(provider, base_url),
       has_key = !is.null(api_key) && length(api_key) == 1 && !is.na(api_key) && nzchar(api_key))
}

ai_key_binding_ok <- function(binding, provider, base_url) {
  !is.null(binding) && identical(binding$provider, provider %||% "gemini") &&
    identical(binding$host, ai_effective_host(provider, base_url))
}

# NULL when the key may be used for this request; otherwise the message to show.
ai_key_binding_message <- function(binding, provider, base_url) {
  if (ai_key_binding_ok(binding, provider, base_url)) return(NULL)
  now <- ai_effective_host(provider, base_url)
  paste0("Not sent: the API key was entered for ",
         if (is.null(binding)) "an unknown endpoint" else
           paste0(ai_provider_field(binding$provider, "label"), " (", binding$host, ")"),
         ", but this request would go to ", ai_provider_field(provider, "label"), " (", now, "). ",
         "Re-enter the API key for that endpoint.")
}

# ==============================================================================
#  ENDPOINT VALIDATION (SSRF guard)
# ==============================================================================

# Canonical dotted-decimal only: no leading zeros (012 is octal 10 to curl and
# the C resolver), no short forms (127.1), no hex (0x7f).
ai_parse_ipv4 <- function(ip) {
  if (length(ip) != 1 || is.na(ip) ||
      !grepl("^(0|[1-9][0-9]{0,2})(\\.(0|[1-9][0-9]{0,2})){3}$", ip)) return(NULL)
  o <- as.integer(strsplit(ip, ".", fixed = TRUE)[[1]])
  if (any(o > 255L)) return(NULL)
  o
}

# Returns the 8 groups of an IPv6 address as integers, or NULL.
ai_parse_ipv6 <- function(ip) {
  if (length(ip) != 1 || is.na(ip)) return(NULL)
  ip <- tolower(sub("%.*$", "", ip))              # drop a zone id (fe80::1%en0)
  if (!grepl(":", ip, fixed = TRUE) || !grepl("^[0-9a-f:.]+$", ip)) return(NULL)
  if (grepl(".", ip, fixed = TRUE)) {              # embedded IPv4 tail (::ffff:1.2.3.4)
    v4 <- ai_parse_ipv4(sub("^.*:", "", ip))
    if (is.null(v4)) return(NULL)
    ip <- paste0(sub("[^:]*$", "", ip),
                 sprintf("%x:%x", v4[1] * 256L + v4[2], v4[3] * 256L + v4[4]))
  }
  n_dc <- lengths(regmatches(ip, gregexpr("::", ip, fixed = TRUE)))
  if (n_dc > 1) return(NULL)
  groups <- function(s) if (!nzchar(s)) character(0) else strsplit(s, ":", fixed = TRUE)[[1]]
  if (n_dc == 1) {
    pos <- regexpr("::", ip, fixed = TRUE)
    left <- groups(substr(ip, 1, pos - 1))
    right <- groups(substr(ip, pos + 2, nchar(ip)))
    fill <- 8L - length(left) - length(right)
    if (fill < 1) return(NULL)
    g <- c(left, rep("0", fill), right)
  } else {
    g <- groups(ip)
  }
  if (length(g) != 8 || !all(grepl("^[0-9a-f]{1,4}$", g))) return(NULL)
  strtoi(g, 16L)
}

ai_ipv4_class <- function(o) {
  a <- o[1]; b <- o[2]; c3 <- o[3]
  if (a == 0) return("blocked")                                   # 0.0.0.0/8 "this host"
  if (a == 169 && b == 254) return("blocked")                     # link-local, incl. 169.254.169.254 metadata
  if (a >= 224) return("blocked")                                 # multicast, reserved, broadcast
  if (a == 192 && b == 0 && c3 %in% c(0, 2)) return("blocked")    # 192.0.0/24 (Oracle metadata 192.0.0.192), TEST-NET-1
  if (a == 198 && b == 51 && c3 == 100) return("blocked")         # TEST-NET-2
  if (a == 203 && b == 0 && c3 == 113) return("blocked")          # TEST-NET-3
  if (a == 127) return("private")                                 # loopback
  if (a == 10) return("private")
  if (a == 172 && b >= 16 && b <= 31) return("private")
  if (a == 192 && b == 168) return("private")
  if (a == 100 && b >= 64 && b <= 127) return("private")          # CGNAT, incl. 100.100.100.200 metadata
  if (a == 198 && b %in% c(18, 19)) return("private")             # benchmarking
  "public"
}

# "public", "private" (loopback / RFC1918 / CGNAT / IPv6 ULA) or "blocked"
# (unparseable, unspecified, link-local, cloud metadata, multicast, reserved).
ai_ip_class <- function(ip) {
  v4 <- ai_parse_ipv4(ip)
  if (!is.null(v4)) return(ai_ipv4_class(v4))
  v6 <- ai_parse_ipv6(ip)
  if (is.null(v6)) return("blocked")
  embedded_v4 <- function(hi, lo) c(hi %/% 256L, hi %% 256L, lo %/% 256L, lo %% 256L)
  if (all(v6 == 0)) return("blocked")                                        # ::
  if (all(v6[1:7] == 0) && v6[8] == 1) return("private")                     # ::1
  if (all(v6[1:5] == 0) && v6[6] == 0xffff) return(ai_ipv4_class(embedded_v4(v6[7], v6[8])))  # ::ffff:a.b.c.d
  if (all(v6[1:6] == 0)) return("blocked")                                   # deprecated ::a.b.c.d
  if (v6[1] == 0x64 && v6[2] == 0xff9b && all(v6[3:6] == 0))
    return(ai_ipv4_class(embedded_v4(v6[7], v6[8])))                          # NAT64
  if (v6[1] == 0x2002) return(ai_ipv4_class(embedded_v4(v6[2], v6[3])))      # 6to4
  if (v6[1] == 0x2001 && v6[2] == 0) return("blocked")                       # Teredo
  if (v6[1] == 0x2001 && v6[2] == 0x0db8) return("blocked")                  # documentation
  if (bitwAnd(v6[1], 0xffc0) == 0xfe80) return("blocked")                    # link-local
  if (bitwAnd(v6[1], 0xffc0) == 0xfec0) return("blocked")                    # site-local (deprecated)
  if (bitwAnd(v6[1], 0xff00) == 0xff00) return("blocked")                    # multicast
  if (bitwAnd(v6[1], 0xfe00) == 0xfc00) return("private")                    # ULA fc00::/7
  "public"
}

ai_ip_allowed <- function(ip, allow_private = FALSE) {
  cls <- ai_ip_class(ip)
  identical(cls, "public") || (isTRUE(allow_private) && identical(cls, "private"))
}

ai_ip_is_loopback <- function(ip) {
  v4 <- ai_parse_ipv4(ip)
  if (!is.null(v4)) return(v4[1] == 127)
  v6 <- ai_parse_ipv6(ip)
  if (is.null(v6)) return(FALSE)
  if (all(v6[1:7] == 0) && v6[8] == 1) return(TRUE)
  all(v6[1:5] == 0) && v6[6] == 0xffff && (v6[7] %/% 256L) == 127
}

# TRUE for a host written as numbers (every label decimal or 0x-hex, or a
# trailing dot) — i.e. something a URL parser or resolver may read as an IPv4
# address in a non-obvious way.
ai_host_is_numeric_form <- function(host) {
  h <- tolower(host)
  labels <- strsplit(h, ".", fixed = TRUE)[[1]]
  length(labels) > 0 && all(grepl("^(0x[0-9a-f]*|[0-9]+)$", labels))
}

# The host exactly as typed in the URL (before any parser normalisation).
ai_raw_url_host <- function(url) {
  auth <- sub("^[A-Za-z][A-Za-z0-9+.-]*://", "", url)
  if (identical(auth, url)) return("")
  auth <- sub("[/?#].*$", "", auth)
  auth <- sub("^.*@", "", auth)
  if (startsWith(auth, "[")) return(sub("\\].*$", "]", auth))
  sub(":[^:]*$", "", auth)
}

ai_default_resolver <- function(host) {
  curl::nslookup(host, ipv4_only = FALSE, multiple = TRUE, error = FALSE)
}

# Validate a user-typed endpoint base URL. Empty input means the provider's
# registry default. Returns list(ok = TRUE, url, scheme, host, port, ips, literal)
# or list(ok = FALSE, error). `resolver` is injectable so tests never touch DNS.
ai_validate_endpoint <- function(base_url, provider = "openai_compat",
                                 policy = ai_deployment_policy(),
                                 resolver = ai_default_resolver) {
  fail <- function(msg) list(ok = FALSE, error = msg)
  url <- if (is.null(base_url) || length(base_url) != 1 || is.na(base_url)) "" else
    trimws(as.character(base_url))
  if (!nzchar(url)) url <- ai_provider_field(provider, "base_url")
  url <- sub("/+$", "", url)
  if (nchar(url) > 2048) return(fail("the endpoint URL is too long."))
  if (grepl("[[:space:][:cntrl:]]", url)) return(fail("the endpoint URL contains spaces or control characters."))

  parsed <- tryCatch(httr2::url_parse(url), error = function(e) NULL)
  host <- if (is.null(parsed)) NULL else parsed$hostname
  if (is.null(host) || length(host) != 1 || is.na(host) || !nzchar(host))
    return(fail("the endpoint URL could not be parsed. Use a form like https://host.example.org/v1"))
  scheme <- tolower(parsed$scheme %||% "")
  if (!scheme %in% c("https", "http"))
    return(fail("the endpoint URL must start with https://"))
  if (!is.null(parsed$username) || !is.null(parsed$password))
    return(fail("the endpoint URL must not contain a user name or password. Put the key in the API Key box."))
  if (length(parsed$query) > 0 || !is.null(parsed$fragment))
    return(fail("the endpoint URL must not contain a query string or #fragment."))

  # The host as typed must be the host the parser returned, and a numeric host
  # must be canonical dotted-decimal. Older httr2 parsed hosts itself while curl
  # reads 012.0.0.1 as 10.0.0.1 and 0x7f.1 as 127.0.0.1; refusing anything
  # ambiguous makes the check independent of the httr2 version.
  raw_host <- tolower(ai_raw_url_host(url))
  if (!identical(raw_host, tolower(host)))
    return(fail("the endpoint host could not be read unambiguously. Write it as a name or as a plain IP address like 192.168.1.50."))
  host <- tolower(gsub("^\\[|\\]$", "", host))
  if (grepl("\\.$", host) || (ai_host_is_numeric_form(host) && is.null(ai_parse_ipv4(host))))
    return(fail("the endpoint host must be a name or a plain dotted IP address like 192.168.1.50 (no leading zeros, hex or short forms)."))
  port <- if (!is.null(parsed$port)) suppressWarnings(as.integer(parsed$port)) else
    if (scheme == "https") 443L else 80L
  if (length(port) != 1 || is.na(port) || port < 1L || port > 65535L)
    return(fail("the endpoint port is not valid."))

  literal <- !is.null(ai_parse_ipv4(host)) || !is.null(ai_parse_ipv6(host))
  ips <- if (literal) host else tryCatch(resolver(host), error = function(e) NULL)
  ips <- unique(as.character(ips))
  ips <- ips[!is.na(ips) & nzchar(ips)]
  allowed <- length(ips) > 0 &&
    all(vapply(ips, ai_ip_allowed, logical(1), allow_private = policy$allow_private_endpoints))
  if (!allowed) {
    # One message for "unresolvable" and "internal" on a public deployment, so
    # the error cannot be used to map which internal names exist.
    if (policy$public)
      return(fail(sprintf("the endpoint '%s' is not allowed here: it must be a public internet address reachable over https.", host)))
    if (length(ips) == 0)
      return(fail(sprintf("could not resolve the endpoint host '%s'.", host)))
    return(fail(sprintf("the endpoint '%s' resolves to a link-local, cloud-metadata or reserved address, which is never allowed.", host)))
  }

  warning <- NULL
  if (scheme == "http") {
    local_net <- all(vapply(ips, function(ip) identical(ai_ip_class(ip), "private"), logical(1)))
    container_host <- host %in% AI_CONTAINER_HOST_ALIASES
    if (!(isTRUE(policy$allow_http_private) && (local_net || container_host))) {
      return(fail(if (policy$public)
        "the endpoint URL must use https:// on this public site."
      else paste("plain http:// is only accepted for a model server on this computer, the local network",
                 "or the container host (e.g. http://localhost:8000/v1, http://192.168.1.50:8000/v1,",
                 "http://host.docker.internal:11434/v1). Use https:// for anything else.")))
    }
    warning <- ai_http_key_warning(host)
  }

  list(ok = TRUE, url = url, scheme = scheme, host = host, port = port,
       ips = ips, literal = literal, warning = warning)
}

# One-line warning for a plain-http endpoint. The single wording, used by the
# validator and by the sidebar notification.
ai_http_key_warning <- function(host) {
  paste0("The endpoint '", host, "' uses plain http: your API key and data are sent unencrypted ",
         "and can be read by anyone on that network path.")
}

# Validation + credential check before any openai_compat request. Returns
# list(ok, endpoint) or list(ok = FALSE, error). The endpoint is resolved FIRST
# so the key check knows the destination host (item: Google's own
# OpenAI-compatible endpoint accepts a Gemini key).
ai_openai_preflight <- function(api_key, base_url, policy = ai_deployment_policy(),
                                resolver = ai_default_resolver) {
  endpoint <- ai_validate_endpoint(base_url, "openai_compat", policy, resolver)
  if (!isTRUE(endpoint$ok)) return(list(ok = FALSE, error = endpoint$error))
  owner <- ai_key_owner("openai_compat", api_key, endpoint$host)
  if (!is.null(owner)) return(list(ok = FALSE, error = ai_foreign_key_message(owner, endpoint$host)))
  list(ok = TRUE, endpoint = endpoint)
}

# Apply the guard to a request: no redirects (a validated host could otherwise
# bounce the request to an internal one) and, for a hostname, pin the addresses
# that were validated so a second DNS answer cannot differ (DNS rebinding).
ai_harden_request <- function(req, endpoint = NULL) {
  req <- httr2::req_options(req, followlocation = 0L)
  if (!is.null(endpoint) && isTRUE(endpoint$ok) && !isTRUE(endpoint$literal) && length(endpoint$ips) > 0) {
    addrs <- vapply(endpoint$ips, function(ip)
      if (is.null(ai_parse_ipv4(ip))) paste0("[", ip, "]") else ip, character(1))
    req <- httr2::req_options(req, resolve = sprintf("%s:%d:%s", endpoint$host, endpoint$port,
                                                     paste(addrs, collapse = ",")))
  }
  req
}

# ==============================================================================
#  ERRORS — a failure is never a result
# ==============================================================================

# Every provider failure is returned as a string starting with this prefix, and
# every consumer checks is_ai_error() before showing, storing or exporting a
# reply. An error must never be rendered as "Analysis Complete", saved as a
# narrative, or fed back to the model as a prior assistant turn.
AI_ERROR_PREFIX <- "API Error:"

ai_error <- function(msg) paste(AI_ERROR_PREFIX, msg)

is_ai_error <- function(x) {
  if (is.null(x) || !is.character(x) || length(x) != 1 || is.na(x)) return(TRUE)
  t <- trimws(x)
  !nzchar(t) || startsWith(t, AI_ERROR_PREFIX)
}

# Upstream text reduced to something safe to show: secrets redacted, no markup,
# no control characters, one line, bounded length.
ai_sanitize_text <- function(x, max_chars = 200, secrets = character(0)) {
  if (is.null(x) || length(x) == 0) return("")
  x <- paste(as.character(x[!is.na(x)]), collapse = " ")
  secrets <- as.character(secrets)
  for (s in secrets[!is.na(secrets) & nchar(secrets) >= 6]) x <- gsub(s, "[redacted]", x, fixed = TRUE)
  x <- gsub("<[^>]*>", " ", x)
  x <- gsub("[[:cntrl:]]+", " ", x)
  x <- gsub("[[:space:]]+", " ", trimws(x))
  if (nchar(x) > max_chars) x <- paste0(substr(x, 1, max_chars), "...")
  x
}

# Condition from req_perform() -> one sanitised "API Error: ..." line.
# HTTP errors: status + reason, plus the provider's own short JSON error message
# when there is one (e.g. "API key not valid"). Never the raw body: on a public
# deployment that body could come from anywhere the URL pointed.
ai_error_from_condition <- function(e, secrets = character(0), timeout_s = NULL) {
  resp <- tryCatch(e$resp, error = function(z) NULL)
  if (inherits(resp, "httr2_response")) {
    status <- tryCatch(httr2::resp_status(resp), error = function(z) NA)
    desc <- tryCatch(httr2::resp_status_desc(resp), error = function(z) NA_character_)
    detail <- tryCatch({
      ct <- httr2::resp_content_type(resp)
      if (length(ct) == 1 && !is.na(ct) && grepl("json", ct, ignore.case = TRUE)) {
        b <- httr2::resp_body_json(resp, simplifyVector = FALSE)
        m <- if (is.list(b$error)) b$error$message else if (is.character(b$error)) b$error else b$message
        if (is.character(m) && length(m) >= 1) m[1] else ""
      } else ""
    }, error = function(z) "")
    line <- paste0("HTTP ", status, if (!is.na(desc) && nzchar(desc)) paste0(" ", desc) else "")
    detail <- ai_sanitize_text(detail, max_chars = 200, secrets = secrets)
    return(ai_error(if (nzchar(detail)) paste0(line, " - ", detail) else line))
  }
  msgs <- c(tryCatch(conditionMessage(e), error = function(z) ""),
            tryCatch(conditionMessage(e$parent), error = function(z) ""))
  full <- paste(msgs[nzchar(msgs)], collapse = " ")
  if (grepl("Timeout was reached|timed out", full, ignore.case = TRUE)) {
    return(ai_error(paste0(
      "the request timed out",
      if (!is.null(timeout_s) && is.finite(timeout_s)) paste0(" after ", round(timeout_s), " s") else "",
      ". Try a faster model, a shorter question, or run DE-LIMP locally for long requests.")))
  }
  ai_error(ai_sanitize_text(if (nzchar(full)) full else "request failed.", max_chars = 200, secrets = secrets))
}

# ==============================================================================
#  PAYLOAD — protein IDs, contaminant flags, tables
# ==============================================================================

# The ONE shortening applied to protein IDs sent to a model: sp|P12345|ALBU_HUMAN
# -> P12345 (first entry of a group; anything else unchanged). The table, the
# selection list and the reverse mapping of [[SELECT: ...]] all use this, so what
# the model sees and what it can send back cannot drift apart.
ai_short_protein_id <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  dec <- grepl("^[a-z]{2}\\|[^|]+\\|", x)
  x[dec] <- sub("^[a-z]{2}\\|([^|]+)\\|.*$", "\\1", x[dec])
  x
}

# Contaminant flag from the single definition (is_contaminant_accession, helpers.R)
# rather than a pattern described to the model in prose (CLAUDE.md rule #3). The
# column is added only when at least one row is a contaminant, so its absence
# means "none of these rows is a flagged contaminant".
ai_flag_contaminants <- function(df, id_col = "Protein") {
  if (is.null(df) || nrow(df) == 0 || !id_col %in% names(df)) return(df)
  flag <- is_contaminant_accession(df[[id_col]])
  if (any(flag)) df$Contaminant <- ifelse(flag, "yes", "")
  df
}

# --- SHARED PAYLOAD FORMATTER ---
# Formats a DE result table for an LLM prompt. Used by BOTH providers so the
# model sees identical data regardless of where the request is sent.
#
# Deliberate reductions (measured ~55% token saving vs. the raw CSV):
#   1. 3 significant figures  — 6dp of a log2 intensity is precision the
#                               measurement does not have, and costs tokens.
#   2. bare accessions        — "sp|P12345|ALBU_HUMAN" -> "P12345"
#   3. an explicit allowlist  — DE statistics, Genes, a Contaminant flag, and the
#                               pipeline's evidence columns passed in `extra_cols`.
#                               Per-sample intensities, Protein.Names (redundant
#                               with Genes) and per-group Mean_/SD_ are NOT sent.
#
# `max_chars` caps the result. Rows are dropped from the BOTTOM, which is the
# least-significant end because callers pass a topTable() ordered by p-value —
# so a trim costs the weakest evidence, never the strongest.
#
# Rows whose Protein is in `pin` (the user's selection) are moved to the top and
# fill the budget FIRST: the prompt tells the model to focus on them. They are
# still bounded by `max_chars`, and each one is also charged for its entry in the
# prompt's selection list (ai_selection_context()), so table + list fit together.
# When a selection is given, the result carries attributes:
#   pinned_sent    — full IDs of selected rows that are in the table
#   pinned_dropped — full IDs of selected rows left out to fit the budget
# The caller passes pinned_sent (and the dropped count) to ask_ai_data() so the
# selection list names exactly what was sent, and tells the user what was not.
format_ai_table <- function(df, sig = 3, max_chars = Inf, pin = NULL, extra_cols = character(0)) {
  core_cols <- c("Protein", "Gene", "Genes",
                 "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  if (is.null(df) || nrow(df) == 0) {
    cols <- if (is.null(df)) character(0) else names(df)
    return(paste(cols[cols %in% c(core_cols, extra_cols, "Contaminant")], collapse = "\t"))
  }

  df <- as.data.frame(df, stringsAsFactors = FALSE)

  # Promote rownames to a Protein column if the caller did not supply one
  if (!"Protein" %in% names(df) && !is.null(rownames(df))) {
    df <- cbind(Protein = rownames(df), df, stringsAsFactors = FALSE)
  }

  # Selected rows first (stable order otherwise)
  pinned <- if (length(pin) > 0 && "Protein" %in% names(df))
    as.character(df$Protein) %in% as.character(pin) else rep(FALSE, nrow(df))
  if (any(pinned)) {
    df <- df[c(which(pinned), which(!pinned)), , drop = FALSE]
    pinned <- c(rep(TRUE, sum(pinned)), rep(FALSE, nrow(df) - sum(pinned)))
  }

  allow <- c(core_cols, as.character(extra_cols), "Contaminant")
  keep <- names(df)[names(df) %in% allow]
  # Preserve a readable order rather than whatever the caller cbind()ed
  df <- df[, intersect(allow, keep), drop = FALSE]

  full_ids <- if ("Protein" %in% names(df)) as.character(df$Protein) else rep("", nrow(df))
  if ("Protein" %in% names(df)) df$Protein <- ai_short_protein_id(df$Protein)

  # Round numerics; render NA as empty so the model does not read "NA" as data
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      # "g" picks the shorter of fixed/scientific per value: logFC stays
      # "-2.12", a p-value stays "1.2e-09". Forcing either wastes tokens.
      out <- formatC(signif(col, sig), digits = sig, format = "g")
      out <- trimws(out)
      out[is.na(col)] <- ""
      out
    } else {
      out <- as.character(col)
      out[is.na(out)] <- ""
      # A tab or newline inside a cell would corrupt the TSV
      gsub("[\t\r\n]+", " ", out)
    }
  })

  header <- paste(names(df), collapse = "\t")
  rows   <- apply(df, 1, function(r) paste(r, collapse = "\t"))

  n_pin <- sum(pinned)
  pin_idx <- seq_len(n_pin)
  rest <- setdiff(seq_along(rows), pin_idx)
  n_pin_keep <- n_pin
  n_keep <- length(rest)
  if (is.finite(max_chars)) {
    # Selected rows first, each charged for its row AND its selection-list entry
    list_overhead <- if (n_pin > 0) nchar(ai_selection_context("X", n_omitted = 99999L)) else 0L
    pin_cost <- nchar(rows[pin_idx]) + 1L + nchar(ai_short_protein_id(full_ids[pin_idx])) + 2L
    running_pin <- nchar(header) + list_overhead + cumsum(pin_cost)
    n_pin_keep <- sum(running_pin <= max_chars)
    used <- nchar(header) + list_overhead + if (n_pin_keep > 0) running_pin[n_pin_keep] - nchar(header) - list_overhead else 0L
    running <- used + cumsum(nchar(rows[rest]) + 1L)
    n_keep  <- sum(running <= max_chars)
    if (n_pin_keep < n_pin || n_keep < length(rest)) {
      message(sprintf("[DE-LIMP] AI payload trimmed to fit %d chars: %d of %d selected rows and %d of %d other rows kept",
                      max_chars, n_pin_keep, n_pin, n_keep, length(rest)))
    }
  }
  out <- paste(c(header, rows[c(pin_idx[seq_len(n_pin_keep)], rest[seq_len(n_keep)])]), collapse = "\n")
  if (length(pin) > 0) {
    attr(out, "pinned_sent") <- full_ids[pin_idx[seq_len(n_pin_keep)]]
    attr(out, "pinned_dropped") <- full_ids[pin_idx[setdiff(seq_len(n_pin), seq_len(n_pin_keep))]]
  }
  out
}

# The selection block of the Data Chat prompt. `sent_ids` must be the selected
# IDs whose rows are actually in the table (format_ai_table's pinned_sent).
ai_selection_context <- function(sent_ids, n_omitted = 0L) {
  shown <- unique(ai_short_protein_id(sent_ids))
  shown <- shown[nzchar(shown)]
  n_omitted <- suppressWarnings(as.integer(n_omitted %||% 0L))
  if (length(n_omitted) != 1 || is.na(n_omitted)) n_omitted <- 0L
  if (length(shown) == 0 && n_omitted == 0) return("")
  paste0(
    "\n!!! USER SELECTION ACTIVE !!!\n",
    if (length(shown) > 0) paste0(
      "Focus analysis on these specific proteins (IDs as in the Protein column; their rows ",
      "are listed first in the DE table):\n", paste(shown, collapse = ", "), "\n") else "",
    if (n_omitted > 0) paste0(
      n_omitted, " further selected protein", if (n_omitted == 1) " was" else "s were",
      " not sent because of the size limit; say so if the user asks about the whole selection.\n") else ""
  )
}

# Column names of a table produced by format_ai_table()
ai_table_columns <- function(data_table) {
  if (is.null(data_table) || !nzchar(data_table)) return(character(0))
  strsplit(strsplit(data_table, "\n", fixed = TRUE)[[1]][1], "\t", fixed = TRUE)[[1]]
}

# How the data prompts describe evidence columns. `evidence` is the named vector
# from pipeline_evidence_columns() (helpers.R); nothing is said about a column
# the pipeline did not supply.
ai_evidence_description <- function(evidence = NULL, guidance = "") {
  if (length(evidence) == 0) return("")
  paste0(paste(sprintf("%s is the %s", names(evidence), evidence), collapse = "; "), ". ",
         if (!is.null(guidance) && nzchar(guidance)) paste0(guidance, " ") else "")
}

# Per-contrast summary for the AI Summary / Claude export prompts.
# `tt` is a topTable-like frame with Protein.Group, logFC, adj.P.Val and optional
# Gene + evidence columns. Sends BOTH views the rules ask the model to keep apart:
# the most significant hits and the largest effect sizes in each direction.
ai_contrast_summary_text <- function(tt, contrast, top_n = 30, top_fc_n = ceiling(top_n / 2),
                                     evidence_cols = character(0), sig_threshold = 0.05) {
  tt <- as.data.frame(tt, stringsAsFactors = FALSE)
  sig <- tt[!is.na(tt$adj.P.Val) & tt$adj.P.Val < sig_threshold & !is.na(tt$logFC), , drop = FALSE]
  n_up <- sum(sig$logFC > 0)
  n_down <- sum(sig$logFC < 0)
  cols <- c("Protein.Group", intersect("Gene", names(tt)), "logFC", "adj.P.Val",
            intersect(evidence_cols, names(tt)))
  tbl <- function(d) {
    if (nrow(d) == 0) return("(none)")
    d <- d[, cols, drop = FALSE]
    names(d)[1] <- "Protein"
    d <- ai_flag_contaminants(d, "Protein")
    format_ai_table(d, extra_cols = evidence_cols)
  }
  by_sig <- utils::head(sig[order(sig$adj.P.Val), , drop = FALSE], top_n)
  up <- sig[sig$logFC > 0, , drop = FALSE]
  up <- utils::head(up[order(-up$logFC), , drop = FALSE], top_fc_n)
  down <- sig[sig$logFC < 0, , drop = FALSE]
  down <- utils::head(down[order(down$logFC), , drop = FALSE], top_fc_n)

  paste0(
    "### ", contrast, "\n",
    "Significant proteins (adj.P.Val < ", sig_threshold, "): ", nrow(sig),
    " (", n_up, " up, ", n_down, " down)\n\n",
    "Top ", nrow(by_sig), " by SIGNIFICANCE (smallest adj.P.Val):\n", tbl(by_sig), "\n\n",
    "Largest INCREASES by EFFECT SIZE (highest logFC among significant proteins; top ",
    nrow(up), "):\n", tbl(up), "\n\n",
    "Largest DECREASES by EFFECT SIZE (lowest logFC among significant proteins; top ",
    nrow(down), "):\n", tbl(down)
  )
}

# ==============================================================================
#  MEMORY — what the client resends, because the model remembers nothing
# ==============================================================================
#
#  No LLM has memory between API calls: every request is stateless and the model
#  sees only the `messages` array. What looks like memory in any chat product is
#  the client resending history. DE-LIMP kept `values$chat_history` for the UI
#  and the download button but never sent it, so a follow-up like "why is that?"
#  reached the model with no referent.
#
#  ONLY two sources feed memory here, both of which are incapable of fabricating:
#    - prior turns of THIS conversation (what was actually said)
#    - the user's own project name and notes from the activity log, for the
#      loaded dataset only (activity_project_context(), helpers_search.R)
#
#  Model conclusions are deliberately NOT persisted across sessions. A confident
#  fabrication ("P30838 is a hemoglobin" - measured, and wrong) would become
#  durable project knowledge and compound. That is CLAUDE.md rule #2 at scale.
# ==============================================================================

# Convert DE-LIMP's chat_history into OpenAI messages, oldest first.
#
# Strict chat templates (vLLM Mistral / Gemma / Llama) reject a messages array
# that does not alternate user/assistant after the system prompt, so the result:
#   - never contains error replies or app notices (they are not conversation,
#     and a failed turn would otherwise be "remembered" as the model's answer);
#   - starts with a user turn and ends with an assistant turn (the caller appends
#     the current question as the next user turn);
#   - never has two consecutive turns from the same role — of two user turns in a
#     row (the first one's answer failed), only the later is kept.
ai_history_messages <- function(chat_history, max_turns = 6, max_chars = 6000) {
  if (is.null(chat_history) || length(chat_history) == 0) return(list())

  clean <- Filter(Negate(is.null), lapply(chat_history, function(m) {
    if (!is.list(m)) return(NULL)
    txt <- m$content %||% ""
    if (length(txt) != 1 || is.na(txt) || !nzchar(trimws(txt))) return(NULL)
    if (isTRUE(m$error) || isTRUE(m$notice)) return(NULL)
    is_user <- identical(m$role, "user")
    if (!is_user && is_ai_error(txt)) return(NULL)
    list(role = if (is_user) "user" else "assistant", content = as.character(txt))
  }))

  normalize <- function(turns) {
    out <- list()
    for (t in turns) {
      if (length(out) == 0) {
        if (t$role == "user") out <- list(t)
      } else if (identical(out[[length(out)]]$role, t$role)) {
        out[[length(out)]] <- t            # keep the later of two same-role turns
      } else {
        out[[length(out) + 1]] <- t
      }
    }
    if (length(out) > 0 && identical(out[[length(out)]]$role, "user")) out <- out[-length(out)]
    out
  }
  clean <- normalize(clean)
  if (length(clean) == 0) return(list())

  # Keep the most recent turns: drop from the front, whole exchanges at a time
  if (length(clean) > max_turns) {
    clean <- normalize(clean[(length(clean) - max_turns + 1):length(clean)])
  }

  # Then trim to the character budget, again dropping the oldest exchange first
  total <- function(x) sum(vapply(x, function(t) nchar(t$content), numeric(1)))
  while (length(clean) > 2 && total(clean) > max_chars) {
    clean <- normalize(clean[-(1:2)])
  }
  # A single surviving exchange may still exceed the budget on its own: shorten
  # the longer turn(s) rather than drop the whole exchange
  if (length(clean) > 0 && total(clean) > max_chars) {
    per_turn <- floor(max_chars / length(clean))
    for (i in seq_along(clean)) {
      if (nchar(clean[[i]]$content) > per_turn) clean[[i]]$content <- substr(clean[[i]]$content, 1, per_turn)
    }
  }
  clean
}

# The Gemini path takes one text blob, not a messages array.
ai_history_text <- function(chat_history, max_turns = 6, max_chars = 6000) {
  m <- ai_history_messages(chat_history, max_turns, max_chars)
  if (length(m) == 0) return("")
  paste0("\n--- EARLIER IN THIS CONVERSATION ---\n",
         paste(vapply(m, function(x)
           paste0(if (x$role == "user") "User: " else "Assistant: ", x$content),
           character(1)), collapse = "\n"),
         "\n--- END ---\n")
}

# Plain-text transcript for the chat download. Each AI turn names the provider
# and model that produced it (rule #1); errors and app notices are labelled as
# such rather than attributed to a model.
ai_chat_transcript <- function(chat_history) {
  if (is.null(chat_history) || length(chat_history) == 0) return(character(0))
  vapply(chat_history, function(msg) {
    who <- if (identical(msg$role, "user")) {
      "YOU"
    } else if (isTRUE(msg$notice)) {
      "DE-LIMP"
    } else {
      src <- msg$source %||% "provider and model not recorded"
      if (isTRUE(msg$error) || is_ai_error(msg$content)) paste0("ERROR (", src, ")") else paste0("AI (", src, ")")
    }
    paste0(who, ": ", msg$content %||% "", "\n---\n")
  }, character(1))
}

# Project context from the user's OWN notes. Emits nothing when there is nothing
# to say - never a placeholder, never an invented default (rule #2).
format_project_context <- function(project = NULL, notes = NULL, max_chars = 1500) {
  ok <- function(x) !is.null(x) && length(x) == 1 && !is.na(x) && nzchar(trimws(as.character(x)))
  parts <- character(0)
  if (ok(project)) parts <- c(parts, paste0("Project: ", trimws(as.character(project))))
  if (ok(notes)) {
    n <- trimws(as.character(notes))
    if (nchar(n) > max_chars) n <- paste0(substr(n, 1, max_chars), " [truncated]")
    parts <- c(parts, paste0("Notes: ", n))
  }
  if (length(parts) == 0) return("")
  paste0("\n--- PROJECT CONTEXT (recorded by the user, not generated by AI) ---\n",
         paste(parts, collapse = "\n"), "\n--- END ---\n")
}

# ==============================================================================
#  [[SELECT: ...]] — the model asking the app to highlight proteins
# ==============================================================================

# Split a reply into the text to show and the IDs the model asked to select.
ai_parse_select_directive <- function(reply) {
  if (is.null(reply) || length(reply) != 1 || is.na(reply)) return(list(text = reply %||% "", ids = character(0)))
  m <- regmatches(reply, gregexpr("\\[\\[SELECT:.*?\\]\\]", reply, perl = TRUE))[[1]]
  if (length(m) == 0) return(list(text = reply, ids = character(0)))
  raw <- gsub("^\\[\\[SELECT:|\\]\\]$", "", m)
  ids <- trimws(unlist(strsplit(raw, "[,;]")))
  ids <- unique(ids[nzchar(ids)])
  list(text = trimws(gsub("\\[\\[SELECT:.*?\\]\\]", "", reply, perl = TRUE)), ids = ids)
}

# Map IDs the model returned (it saw ai_short_protein_id() forms) back to the
# app's real row names. Tries the exact row name, then the shortened form, then
# any member of a ";"-separated protein group. Unmatched IDs are reported, not
# silently dropped, so the chat never claims to have highlighted nothing.
ai_match_protein_ids <- function(ids, full_ids) {
  ids <- unique(trimws(as.character(ids)))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  full_ids <- as.character(full_ids)
  short <- ai_short_protein_id(full_ids)
  members <- lapply(strsplit(full_ids, ";", fixed = TRUE), function(p) ai_short_protein_id(trimws(p)))
  matched <- character(0)
  unmatched <- character(0)
  for (id in ids) {
    hit <- full_ids[full_ids == id]
    if (length(hit) == 0) hit <- full_ids[short == id]
    if (length(hit) == 0) hit <- full_ids[vapply(members, function(p) id %in% p, logical(1))]
    if (length(hit) > 0) matched <- c(matched, hit) else unmatched <- c(unmatched, id)
  }
  list(matched = unique(matched), unmatched = unmatched)
}

# ==============================================================================
#  AI SUMMARY PROMPT
# ==============================================================================

# `ctx` is build_ai_data_context()'s list. `evidence` / `evidence_guidance` come
# from the pipeline descriptor (pipeline_evidence_columns() /
# pipeline_evidence_guidance(), helpers.R) — CLAUDE.md rule #1: the prompt only
# describes measurement-depth columns the pipeline that ran actually produced.
# The EVIDENCE STRENGTH rule, in the one wording every DE prompt uses (AI Summary
# and the Claude export). `evidence` / `guidance` come from
# pipeline_evidence_columns() / pipeline_evidence_guidance() (CLAUDE.md rule #3).
ai_evidence_rule <- function(number, evidence = NULL, guidance = "", pipeline = NULL,
                             where = "the tables") {
  pipe_txt <- if (!is.null(pipeline) && length(pipeline) == 1 && !is.na(pipeline) && nzchar(pipeline))
    paste0(" (quantification pipeline: ", pipeline, ")") else ""
  if (length(evidence) > 0) {
    paste0(number, ". EVIDENCE STRENGTH. ", where, " include ", paste(names(evidence), collapse = " and "),
           ": ", ai_evidence_description(evidence, guidance),
           "Cite them when calling a result reliable or unreliable.\n")
  } else {
    paste0(number, ". EVIDENCE STRENGTH. No per-protein measurement-depth statistics (such as precursor ",
           "counts) are supplied for this analysis", pipe_txt, ". Do not describe any hit as ",
           "well or poorly measured on grounds the data does not show.\n")
  }
}

# What a stable-biomarker assessment may weigh: measurement depth only when the
# pipeline supplied it.
ai_biomarker_criteria <- function(evidence = NULL) {
  if (length(evidence) > 0) "low CV, significant p-value, meaningful fold-change and measurement depth" else
    "low CV, significant p-value and meaningful fold-change"
}

build_ai_summary_prompt <- function(ctx, evidence = NULL, evidence_guidance = "",
                                    pipeline = NULL) {
  has_ev <- length(evidence) > 0
  rule5 <- ai_evidence_rule(5, evidence, evidence_guidance, pipeline, where = "The tables")
  evidence_section <- if (has_ev) {
    paste0("For the headline hits, assess the strength of the underlying measurement (rule 5). Name any ",
           "hit whose statistics look strong but whose measurement is thin.")
  } else {
    paste0("Measurement-depth statistics were not supplied (rule 5), so do not grade the headline hits ",
           "on measurement depth. Comment only on what the statistics shown support.")
  }

  system_prompt <- paste0(
    "You are a senior proteomics and systems biology consultant. Write a comprehensive ",
    "analysis of the differential expression results across ALL comparisons below.\n\n",

    "## RULES - these override any stylistic instruction below\n\n",
    "1. IDENTITY. Name a protein ONLY using the gene name supplied in the data (the Gene column). Never ",
    "from memory. Where no gene name is supplied, use the accession alone and say the ",
    "identity was not provided.\n",
    "2. BIOLOGY. Discuss function, pathway or disease association ONLY for proteins whose ",
    "identity was supplied. If you cannot support a claim from the data given, omit it. ",
    "Do not write biology you merely recognise.\n",
    "3. NUMBERS. Do not state any numeric fact - fold-change, p-value, residue count, ",
    "molecular weight - unless it appears in the data below.\n",
    "4. EFFECT SIZE vs SIGNIFICANCE. 'Most increased' and 'most decreased' mean the ",
    "largest and smallest logFC. That is a question about effect size, not significance; ",
    "the protein with the best p-value is often a different one. Each comparison below lists ",
    "the top hits by significance AND, separately, the largest increases and decreases by logFC. ",
    "Take effect-size claims from the logFC lists. If they differ, say both.\n",
    rule5,
    "6. TECHNICAL vs BIOLOGICAL. If a comparison is strongly one-sided, consider whether ",
    "it reflects a global normalisation or loading difference rather than biology, and ",
    "say which you think it is. If the contrast is a method or sample-preparation ",
    "comparison rather than a biological one, say so plainly instead of constructing a ",
    "biological narrative.\n",
    "7. CONTAMINANTS. A row with Contaminant = yes was flagged by DE-LIMP's contaminant list. ",
    "A table without a Contaminant column contains no flagged contaminants. Do not decide from ",
    "an accession or name yourself whether a protein is a contaminant.\n\n",

    "Structure your response with these markdown sections:\n\n",
    "## Overview\n",
    "Number of comparisons analyzed, total significant proteins per comparison (up/down split). ",
    "Overall assessment of the experiment's quality and scope, including whether any comparison ",
    "looks technical rather than biological (rule 6).\n\n",
    "## Key Findings Per Comparison\n",
    "For each comparison: highlight the top upregulated and downregulated proteins by fold-change ",
    "(from the logFC lists; use the supplied gene names). Note any comparison with unusually few or many significant hits.\n\n",
    "## Evidence Quality\n",
    evidence_section, " Also name any flagged contaminant entries (rule 7) that reached significance.\n\n",
    "## Cross-Comparison Biomarkers\n",
    "Proteins significant in multiple comparisons are highest-confidence candidates. ",
    "Discuss consistency of direction (always up, always down, or mixed across comparisons).\n\n",
    "## High-Confidence Biomarker Insights\n",
    "For the most stable proteins (lowest coefficient of variation): assess their potential as ",
    "reliable biomarkers based on the combination of ", ai_biomarker_criteria(evidence),
    ". Discuss biology only within rule 2.\n\n",
    "## Biological Interpretation\n",
    "Suggest what biological processes or pathways may be affected, within rule 2. ",
    "If the data does not support a biological narrative, say so rather than constructing one.\n\n",
    "Use markdown formatting with headers. Be scientific but accessible."
  )

  paste0(
    system_prompt,
    "\n\n--- DATA FOR ANALYSIS ---\n\n",
    "Number of comparisons: ", ctx$n_contrasts, "\n\n",
    ctx$contrast_text, "\n\n",
    "--- CROSS-COMPARISON PROTEINS (significant in >= 2 comparisons) ---\n",
    ctx$cross_text, "\n\n",
    "--- MOST STABLE SIGNIFICANT PROTEINS (lowest CV across replicates) ---\n",
    ctx$stable_prots_text
  )
}

# ==============================================================================
#  OPENAI-COMPATIBLE PROVIDER (/v1/chat/completions)
# ==============================================================================

# Strip reasoning traces. Qwen3 and other thinking models served through vLLM
# emit <think>...</think> inline; some gateways return a separate
# reasoning_content field. Neither belongs in a user-facing summary.
#
# Some chat templates inject the OPENING <think> into the prompt, so the reply
# contains only the trace and a lone closing </think> — everything up to that
# closing tag is reasoning too.
strip_reasoning <- function(txt) {
  if (is.null(txt) || length(txt) != 1 || is.na(txt) || !nzchar(txt)) return("")
  txt <- gsub("(?s)<think>.*?</think>", "", txt, perl = TRUE)
  # A closing tag with no opening one: the template supplied <think>
  if (grepl("</think>", txt, fixed = TRUE)) txt <- sub("(?s)^.*</think>", "", txt, perl = TRUE)
  # An unterminated <think> means the model ran out of budget mid-thought
  txt <- gsub("(?s)<think>.*$", "", txt, perl = TRUE)
  trimws(txt)
}

# Turn one /chat/completions message object into the text a user should see.
#
# Reasoning models expose their chain-of-thought either inline as <think>...
# </think> or, on this gateway, as a separate `reasoning_content` field. Both
# are internal deliberation. Neither may reach a scientific summary — a reader
# would take "let me check the logFC column" as a finding. If the model spent
# its whole budget thinking and produced no answer, say so plainly rather than
# dressing the trace up as analysis (CLAUDE.md rule #2: never silently
# substitute a fabricated value for a missing one).
extract_openai_reply <- function(msg) {
  if (is.null(msg)) return(ai_error("endpoint returned an empty response."))
  content <- msg$content
  content <- if (is.character(content) && length(content) == 1 && !is.na(content)) strip_reasoning(content) else ""
  if (nzchar(content)) return(content)

  reasoning <- msg$reasoning_content %||% ""
  if (is.character(reasoning) && length(reasoning) == 1 && !is.na(reasoning) && nzchar(reasoning)) {
    return(ai_error(paste("the model used its entire output budget on internal",
                          "reasoning and returned no answer. Raise max_tokens, or choose a",
                          "model with a smaller reasoning overhead.")))
  }
  ai_error("endpoint returned an empty response.")
}

# max_tokens covers reasoning AND the answer on this gateway. Measured: 4096
# was entirely consumed by deepseek-v4-flash's internal reasoning on a Data Chat
# prompt, leaving an empty answer. 8192 leaves room for both.
ask_openai_compat <- function(system_prompt, user_prompt, api_key, model_name,
                              base_url = NULL, max_tokens = 8192,
                              timeout_s = NULL, history = list(),
                              policy = ai_deployment_policy()) {
  pre <- ai_openai_preflight(api_key, base_url, policy)
  if (!isTRUE(pre$ok)) return(ai_error(pre$error))
  endpoint <- pre$endpoint
  timeout_s <- ai_request_timeout("openai_compat", timeout_s, policy)

  body <- list(
    model = model_name,
    messages = c(
      list(list(role = "system", content = system_prompt)),
      history,                                    # prior turns, oldest first
      list(list(role = "user", content = user_prompt))
    ),
    temperature = 0.3,
    max_tokens  = max_tokens
  )

  req <- request(paste0(endpoint$url, "/chat/completions")) %>%
    req_headers(
      "Authorization" = paste("Bearer", api_key),
      "Content-Type"  = "application/json"
    ) %>%
    req_body_json(body) %>%
    req_timeout(timeout_s) %>%
    ai_harden_request(endpoint)

  tryCatch({
    resp <- req_perform(req)
    payload <- resp_body_json(resp)
    if (length(payload$choices) == 0) return(ai_error("endpoint returned no choices."))
    extract_openai_reply(payload$choices[[1]]$message)
  }, error = function(e) ai_error_from_condition(e, secrets = api_key, timeout_s = timeout_s))
}

list_openai_compat_models <- function(api_key, base_url = NULL, policy = ai_deployment_policy()) {
  pre <- ai_openai_preflight(api_key, base_url, policy)
  if (!isTRUE(pre$ok)) return(ai_error(pre$error))
  endpoint <- pre$endpoint
  timeout_s <- min(60, policy$max_timeout_s)
  req <- request(paste0(endpoint$url, "/models")) %>%
    req_headers("Authorization" = paste("Bearer", api_key)) %>%
    req_timeout(timeout_s) %>%
    ai_harden_request(endpoint)
  tryCatch({
    data <- resp_body_json(req_perform(req))
    ids <- vapply(data$data %||% list(), function(x) as.character(x$id %||% ""), character(1))
    ids[nzchar(ids)]
  }, error = function(e) ai_error_from_condition(e, secrets = api_key, timeout_s = timeout_s))
}

# ==============================================================================
#  DISPATCHERS — the ONLY entry points call sites may use
# ==============================================================================

# Plain text prompt, no data table (AI Summary, Comparator hypothesis engine)
ask_ai_text <- function(user_query, api_key, model_name, provider = "gemini",
                        base_url = NULL, timeout_s = NULL, chat_history = NULL,
                        policy = ai_deployment_policy()) {
  if (identical(provider, "openai_compat")) {
    ask_openai_compat(
      system_prompt = "You are a PhD-level expert in proteomics and systems biology.",
      user_prompt   = user_query,
      api_key = api_key, model_name = model_name, base_url = base_url,
      timeout_s = timeout_s, history = ai_history_messages(chat_history),
      policy = policy
    )
  } else {
    ask_gemini_text_chat(paste0(ai_history_text(chat_history), user_query),
                         api_key, model_name, timeout_s = timeout_s, policy = policy)
  }
}

# Data-grounded chat. `data_table` is pre-formatted text from format_ai_table().
# Both providers inline it: at ~10k tokens there is nothing for Gemini's File
# API to solve, and one code path is easier to reason about than two.
#
# `evidence` / `evidence_guidance`: the pipeline's evidence columns and how to
# read them (pipeline_evidence_columns() / pipeline_evidence_guidance()). The
# prompt describes only columns that are both supplied by the pipeline AND
# present in `data_table`.
ask_ai_data <- function(user_query, data_table, qc_df, api_key, model_name,
                        provider = "gemini", selected_ids = NULL, base_url = NULL,
                        timeout_s = NULL, chat_history = NULL,
                        project = NULL, notes = NULL,
                        evidence = NULL, evidence_guidance = "",
                        n_selected_omitted = 0L,
                        policy = ai_deployment_policy()) {

  qc_text <- if (is.null(qc_df)) "No QC data available." else
    paste(capture.output(write.table(qc_df, sep = "\t", row.names = FALSE, quote = FALSE)), collapse = "\n")

  cols <- ai_table_columns(data_table)
  if (length(evidence) > 0) {
    n_declared <- length(evidence)
    evidence <- evidence[names(evidence) %in% cols]
    # The guidance sentence names every evidence column; drop it if any is absent
    if (length(evidence) < n_declared) evidence_guidance <- ""
  }
  if (length(evidence) == 0) {
    evidence <- NULL
    evidence_guidance <- ""
  }

  # `selected_ids` = the selected proteins actually present in data_table
  selection_context <- ai_selection_context(selected_ids, n_selected_omitted)

  column_notes <- paste0(
    "logFC is log2 fold-change; adj.P.Val is BH-adjusted; B is the log-odds of differential ",
    "expression (higher = stronger evidence, NOT a batch effect). ",
    if ("Genes" %in% cols) "Genes is the gene name from the search FASTA. " else
      "No gene names are supplied in this table. ",
    if (length(evidence) > 0) paste0(
      ai_evidence_description(evidence, evidence_guidance),
      "Report results that rest on thin measurement with that caveat. ") else
      "No per-protein measurement-depth statistics are supplied; do not claim a hit is well or poorly measured. ",
    if ("Contaminant" %in% cols)
      "Contaminant = yes marks an entry flagged by DE-LIMP's contaminant list; do not judge contamination from names yourself. " else
      "None of the rows is flagged as a contaminant by DE-LIMP's contaminant list. "
  )

  system_prompt <- paste0(
    "You are a PhD-level expert in proteomics. You have two data sources.\n\n",
    "SOURCE 1 - QC METRICS (tab-separated):\n",
    "Use the 'Group' column to compare technical quality (Precursors, MS1) between groups.\n",
    "--- START QC DATA ---\n", qc_text, "\n--- END QC DATA ---\n\n",
    "SOURCE 2 - DIFFERENTIAL EXPRESSION RESULTS (tab-separated):\n",
    "These are the output of limpa/limma. The statistics are already computed - ",
    "interpret them, do not recompute. ", column_notes, "\n\n",
    "RULE ON EFFECT SIZE vs SIGNIFICANCE: 'most increased' and 'most decreased' ",
    "mean the LARGEST and SMALLEST logFC. That is a question about effect size, ",
    "NOT about significance - the protein with the best p-value or the highest B ",
    "is often a different one. If they differ, report both and label which is ",
    "which.\n",
    "RULE ON PROTEIN IDENTITY: name a protein ONLY using the Genes column supplied ",
    "above. If that column is absent or blank for a row, refer to it by accession alone ",
    "and say the identity was not supplied. ",
    "Never supply a protein name, family, or function from memory - a wrong ",
    "identity invalidates the interpretation built on it.\n",
    "--- START DE DATA ---\n", data_table, "\n--- END DE DATA ---\n\n",
    "BI-DIRECTIONAL CONTROL:\n",
    "1. If the user asks about 'selected proteins', see the USER SELECTION section.\n",
    "2. If you find interesting proteins, output their IDs exactly as written in the Protein ",
    "column at the end like this:\n",
    "   [[SELECT: P12345; P67890]]\n"
  )

  proj_context <- format_project_context(project, notes)

  if (identical(provider, "openai_compat")) {
    ask_openai_compat(
      system_prompt = paste0(system_prompt, proj_context, selection_context),
      user_prompt   = user_query,
      api_key = api_key, model_name = model_name, base_url = base_url,
      timeout_s = timeout_s, history = ai_history_messages(chat_history),
      policy = policy
    )
  } else {
    ask_gemini_text_chat(
      paste0(system_prompt, proj_context, selection_context,
             ai_history_text(chat_history), "\n\nUser Question: ", user_query),
      api_key, model_name, timeout_s = timeout_s, policy = policy
    )
  }
}

# Model listing, dispatched by provider
list_ai_models <- function(api_key, provider = "gemini", base_url = NULL,
                           policy = ai_deployment_policy()) {
  if (identical(provider, "openai_compat")) {
    list_openai_compat_models(api_key, base_url, policy)
  } else {
    list_google_models(api_key, policy)
  }
}

# --- CHECK AVAILABLE MODELS ---
list_google_models <- function(api_key, policy = ai_deployment_policy()) {
  timeout_s <- min(60, policy$max_timeout_s)
  req <- request(paste0(ai_provider_field("gemini", "base_url"), "/models")) %>%
    req_url_query(key = api_key) %>%
    req_timeout(timeout_s) %>%
    ai_harden_request()
  tryCatch({
    resp <- req_perform(req)
    data <- resp_body_json(resp)
    models <- vapply(data$models %||% list(), function(x) as.character(x$name %||% ""), character(1))
    models <- gsub("^models/", "", models)
    models[nzchar(models)]
  }, error = function(e) ai_error_from_condition(e, secrets = api_key, timeout_s = timeout_s))
}

# ==============================================================================
#  GEMINI
# ==============================================================================

# Gemini puts the model name in the URL path. Anything but a plain model id is
# refused so a typed "model name" cannot rewrite the request path.
ai_gemini_model_ok <- function(model_name) {
  m <- gsub("^models/", "", as.character(model_name %||% ""))
  length(m) == 1 && !is.na(m) && grepl("^[A-Za-z0-9._-]+$", m)
}

# A narrative restored from a saved session. Sessions saved before 4.1.0 could
# store a failed request ("API Error: <upstream body>") as the narrative; it is
# not shown or exported as analysis, only reported as unavailable.
AI_SAVED_ERROR_LABEL <- "not available (saved error)"
ai_restored_narrative <- function(narrative, source = NULL) {
  if (is.null(narrative)) return(list(narrative = NULL, source = NULL))
  if (is_ai_error(narrative)) return(list(narrative = NULL, source = AI_SAVED_ERROR_LABEL))
  list(narrative = narrative, source = source)
}

# The text of a generateContent response. Thought parts (thinking models) are
# excluded; an empty or blocked answer is an error, never a NULL "result".
extract_gemini_reply <- function(payload) {
  cands <- payload$candidates
  if (length(cands) == 0) {
    reason <- payload$promptFeedback$blockReason %||% ""
    return(ai_error(if (is.character(reason) && nzchar(reason))
      paste0("Gemini returned no answer (prompt blocked: ", ai_sanitize_text(reason, 60), ").") else
      "Gemini returned no answer."))
  }
  parts <- cands[[1]]$content$parts %||% list()
  txt <- paste(vapply(parts, function(p) {
    if (isTRUE(p$thought) || !is.character(p$text)) "" else paste(p$text, collapse = "")
  }, character(1)), collapse = "")
  txt <- strip_reasoning(txt)
  if (!nzchar(txt)) {
    fr <- cands[[1]]$finishReason %||% ""
    return(ai_error(if (is.character(fr) && nzchar(fr))
      paste0("Gemini returned an empty answer (finish reason: ", ai_sanitize_text(fr, 60), ").") else
      "Gemini returned an empty answer."))
  }
  txt
}

# --- FILE API UPLOADER (LEGACY — no longer on the Data Chat path) ---
# Kept for reference and any out-of-band use. Since v4.1.0 the DE table is
# trimmed by format_ai_table() to ~10k tokens and inlined for both providers,
# so there is nothing left for the File API to solve. Note the File API never
# reduced context cost — uploaded files are tokenised into the request exactly
# like inline text; it only avoided re-uploading bytes.
upload_csv_to_gemini <- function(df, api_key, policy = ai_deployment_policy()) {
  temp_file <- tempfile(fileext = ".csv")
  write.csv(df, temp_file, row.names = FALSE)
  file_size <- file.size(temp_file)

  req <- request("https://generativelanguage.googleapis.com/upload/v1beta/files") %>%
    req_url_query(key = api_key) %>%
    req_headers(
      "X-Goog-Upload-Protocol" = "raw",
      "X-Goog-Upload-Command" = "start, upload, finalize",
      "X-Goog-Upload-Header-Content-Length" = as.character(file_size),
      "X-Goog-Upload-Header-Content-Type" = "text/csv",
      "Content-Type" = "text/csv"
    ) %>%
    req_body_file(temp_file) %>%
    req_timeout(ai_request_timeout("gemini", NULL, policy)) %>%
    ai_harden_request()

  resp <- req_perform(req)
  file_info <- resp_body_json(resp)
  return(file_info$file$uri)
}

# --- AI CHAT FUNCTION (LEGACY — File API path, not used by the app) ---
ask_gemini_file_chat <- function(user_query, file_uri, qc_df, api_key, model_name,
                                 selected_ids = NULL, timeout_s = NULL,
                                 policy = ai_deployment_policy()) {
  if (!ai_gemini_model_ok(model_name)) return(ai_error("the model name is not a valid Gemini model id."))
  timeout_s <- ai_request_timeout("gemini", timeout_s, policy)

  qc_text <- "No QC Data Available"
  if(!is.null(qc_df)) {
    qc_text <- paste(capture.output(write.csv(qc_df, row.names=FALSE)), collapse="\n")
  }

  selection_context <- ""
  if (!is.null(selected_ids) && length(selected_ids) > 0) {
    selection_context <- paste0(
      "\n!!! URGENT: USER SELECTION ACTIVE !!!\n",
      "Focus analysis on these specific proteins:\n",
      paste(selected_ids, collapse=", "), "\n"
    )
  }

  system_instruction <- paste0(
    "You are a PhD-level expert in proteomics. ",
    "You have access to two data sources:\n",
    "SOURCE 1: QC METRICS (In Text Below)\n",
    "This table includes 'Group' columns. Use it to compare technical quality (Precursors, MS1) between experimental groups.\n",
    "--- START QC DATA ---\n",
    qc_text,
    "\n--- END QC DATA ---\n\n",
    "SOURCE 2: EXPRESSION DATA (In Uploaded File)\n",
    "Use this file to answer biological questions. It contains the Top 800 proteins.\n\n",
    "IMPORTANT: BI-DIRECTIONAL CONTROL.\n",
    "1. If the user asks about 'selected proteins', refer to the 'URGENT' section below.\n",
    "2. If you find interesting proteins, OUTPUT their IDs at the end like this:\n",
    "   [[SELECT: P12345; P67890]]\n"
  )

  full_url <- paste0(ai_provider_field("gemini", "base_url"), "/models/",
                     gsub("^models/", "", model_name), ":generateContent")

  body <- list(contents = list(list(parts = list(
    list(text = paste0(system_instruction, selection_context, "\n\nUser Question: ", user_query)),
    list(file_data = list(file_uri = file_uri, mime_type = "text/csv"))
  ))))

  req <- request(full_url) %>%
    req_url_query(key = api_key) %>%
    req_headers("Content-Type" = "application/json") %>%
    req_body_json(body) %>%
    req_timeout(timeout_s) %>%
    ai_harden_request()

  tryCatch({
    extract_gemini_reply(resp_body_json(req_perform(req)))
  }, error = function(e) ai_error_from_condition(e, secrets = api_key, timeout_s = timeout_s))
}

# --- AI TEXT CHAT FUNCTION ---
ask_gemini_text_chat <- function(user_query, api_key, model_name, timeout_s = NULL,
                                 policy = ai_deployment_policy()) {
  if (!ai_gemini_model_ok(model_name)) return(ai_error("the model name is not a valid Gemini model id."))
  timeout_s <- ai_request_timeout("gemini", timeout_s, policy)
  full_url <- paste0(ai_provider_field("gemini", "base_url"), "/models/",
                     gsub("^models/", "", model_name), ":generateContent")

  body <- list(contents = list(list(parts = list(list(text = user_query)))))

  req <- request(full_url) %>%
    req_url_query(key = api_key) %>%
    req_headers("Content-Type" = "application/json") %>%
    req_body_json(body) %>%
    req_timeout(timeout_s) %>%
    ai_harden_request()

  tryCatch({
    extract_gemini_reply(resp_body_json(req_perform(req)))
  }, error = function(e) ai_error_from_condition(e, secrets = api_key, timeout_s = timeout_s))
}
