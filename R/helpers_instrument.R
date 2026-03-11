# ==============================================================================
#  helpers_instrument.R — Instrument metadata extraction from raw MS files
#  Pure functions, no Shiny reactivity. Auto-sourced by app.R.
# ==============================================================================

#' Parse timsTOF metadata from analysis.tdf SQLite file
#' @param tdf_path Path to the analysis.tdf file directly
#' @return Named list of instrument metadata, or NULL on failure
parse_timstof_from_tdf <- function(tdf_path) {
  if (!file.exists(tdf_path)) return(NULL)
  if (!requireNamespace("DBI", quietly = TRUE) ||
      !requireNamespace("RSQLite", quietly = TRUE)) {
    return(list(instrument_type = "timsTOF", parse_error = "DBI/RSQLite not available"))
  }

  meta <- list(instrument_type = "timsTOF")

  tryCatch({
    # Copy TDF to temp for SQLite compatibility (network/SMB mounts can't do WAL/sync)
    tmp_tdf <- tempfile(fileext = ".tdf")
    file.copy(tdf_path, tmp_tdf)
    db <- DBI::dbConnect(RSQLite::SQLite(), tmp_tdf, flags = RSQLite::SQLITE_RO)
    on.exit({ DBI::dbDisconnect(db); unlink(tmp_tdf) }, add = TRUE)

    # Read GlobalMetadata key-value table
    if ("GlobalMetadata" %in% DBI::dbListTables(db)) {
      gm <- DBI::dbGetQuery(db, "SELECT Key, Value FROM GlobalMetadata")
      gm_map <- setNames(gm$Value, gm$Key)

      meta$instrument_model    <- gm_map[["InstrumentName"]] %||% NA_character_
      meta$instrument_serial   <- gm_map[["InstrumentSerialNumber"]] %||% NA_character_
      meta$acquisition_software <- paste0(
        gm_map[["AcquisitionSoftware"]] %||% "",
        if (!is.null(gm_map[["AcquisitionSoftwareVersion"]]))
          paste0(" ", gm_map[["AcquisitionSoftwareVersion"]]) else "")
      meta$sample_name         <- gm_map[["SampleName"]] %||% NA_character_
      meta$acquisition_date    <- gm_map[["AcquisitionDateTime"]] %||% NA_character_
      meta$operator            <- gm_map[["OperatorName"]] %||% NA_character_

      # m/z range
      meta$mz_range_low  <- as.numeric(gm_map[["MzAcqRangeLower"]] %||% NA)
      meta$mz_range_high <- as.numeric(gm_map[["MzAcqRangeUpper"]] %||% NA)

      # Ion mobility range
      meta$im_range_low  <- as.numeric(gm_map[["OneOverK0AcqRangeLower"]] %||% NA)
      meta$im_range_high <- as.numeric(gm_map[["OneOverK0AcqRangeUpper"]] %||% NA)
    }

    # Frames table — scan count, RT range, dia-PASEF detection
    if ("Frames" %in% DBI::dbListTables(db)) {
      frame_stats <- DBI::dbGetQuery(db,
        "SELECT COUNT(*) AS n_frames,
                MIN(Time) AS rt_start,
                MAX(Time) AS rt_end
         FROM Frames")
      meta$n_frames    <- frame_stats$n_frames
      meta$rt_start_min <- round(frame_stats$rt_start / 60, 1)
      meta$rt_end_min   <- round(frame_stats$rt_end / 60, 1)

      # Count MS1 and MS2 frames by MsMsType
      # MsMsType: 0 = MS1, 8/9 = dia-PASEF MS2
      ms1_count <- DBI::dbGetQuery(db,
        "SELECT COUNT(*) AS n FROM Frames WHERE MsMsType = 0")
      meta$n_ms1_spectra <- ms1_count$n

      dia_count <- DBI::dbGetQuery(db,
        "SELECT COUNT(*) AS n FROM Frames WHERE MsMsType IN (8, 9)")
      meta$n_ms2_spectra <- dia_count$n
      meta$n_dia_frames <- dia_count$n
      meta$acquisition_mode <- if (dia_count$n > 0) "dia-PASEF" else "unknown"

      # Estimate cycle time from MS1 frames
      if (ms1_count$n > 1 && !is.na(frame_stats$rt_start) && !is.na(frame_stats$rt_end)) {
        gradient_sec <- frame_stats$rt_end - frame_stats$rt_start
        meta$cycle_time_sec <- round(gradient_sec / ms1_count$n, 3)
      }
    }

    # Mass analyzer type (always TOF for timsTOF)
    meta$ms1_mass_analyzer <- "TOF"
    meta$ms2_mass_analyzer <- "TOF"

    # MS method name from GlobalMetadata
    if (exists("gm_map") && !is.null(gm_map[["MethodName"]])) {
      meta$ms_method <- gm_map[["MethodName"]]
    }

    # --- LC info from HyStarMetadata.xml (standard UTF-8, reliable) ---
    tdf_dir <- dirname(tdf_path)
    hystar_xml_path <- file.path(tdf_dir, "HyStarMetadata.xml")
    if (file.exists(hystar_xml_path) && requireNamespace("xml2", quietly = TRUE)) {
      tryCatch({
        hxml <- xml2::read_xml(hystar_xml_path)
        hns <- xml2::xml_ns(hxml)
        ns_prefix <- if (length(hns) > 0) names(hns)[1] else NULL

        # Helper to find parameter by ID
        find_param <- function(id) {
          xpath <- if (!is.null(ns_prefix)) {
            sprintf(".//%s:Parameter[@ID='%s']", ns_prefix, id)
          } else {
            sprintf(".//Parameter[@ID='%s']", id)
          }
          node <- xml2::xml_find_first(hxml, xpath, hns)
          if (is.na(node)) return(NA_character_)
          # Value might be in <Value>, <Int>, <Double>, or <DateTime> child
          for (tag in c("Value", "Int", "Double", "DateTime")) {
            child_xpath <- if (!is.null(ns_prefix)) {
              sprintf("./%s:%s", ns_prefix, tag)
            } else {
              sprintf("./%s", tag)
            }
            child <- xml2::xml_find_first(node, child_xpath, hns)
            if (!is.na(child)) return(xml2::xml_text(child))
          }
          NA_character_
        }

        # LC system name and method
        lc_device <- find_param("ConfigurationObject_DeviceName")
        lc_method <- find_param("MethodDataObject_MethodName")
        lc_runtime <- find_param("MethodDataObject_RunTimeMinutes")
        lc_serial <- find_param("ConfigurationObject_SerialNumber")

        if (!is.na(lc_device) && nzchar(lc_device)) meta$lc_system <- lc_device
        if (!is.na(lc_method) && nzchar(lc_method)) {
          # Clean up raw method path (e.g., "C:\ProgramData\Evosep\...\100-samples-per-day_2025-...")
          raw_method <- lc_method
          # Extract just the filename without timestamp suffix
          method_base <- basename(gsub("\\\\", "/", raw_method))
          # Strip date-time suffix (e.g., "_2025-11-24_18-57-08")
          method_base <- sub("_\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}-\\d{2}$", "", method_base)
          meta$lc_method <- method_base

          # EvoSep SPD auto-detection from method name
          # Handles both "100-samples-per-day" and "100 samples per day"
          spd_match <- regmatches(method_base,
            regexpr("(\\d+)[- ]samples[- ]per[- ]day", method_base, ignore.case = TRUE))
          if (length(spd_match) > 0 && nzchar(spd_match)) {
            spd_val <- as.integer(sub("[- ]samples.*", "", spd_match))
            meta$evosep_spd <- spd_val
            # Map SPD to gradient length (from EvoSep website)
            spd_gradients <- c("30" = 44, "60" = 21, "100" = 11.5,
                               "200" = 5.5, "300" = 2.3, "500" = 2.2)
            grad_min <- spd_gradients[as.character(spd_val)]
            if (!is.na(grad_min)) meta$evosep_gradient_min <- as.numeric(grad_min)
            meta$lc_method <- sprintf("EvoSep %d SPD (%s)", spd_val, method_base)
          }
        }
        if (!is.na(lc_runtime) && nzchar(lc_runtime)) meta$lc_runtime_min <- as.numeric(lc_runtime)
        if (!is.na(lc_serial) && nzchar(lc_serial)) meta$lc_serial <- lc_serial

        # Also get timsControl version from the MS plugin
        ms_software <- find_param("MS Control")
        if (!is.na(ms_software) && nzchar(ms_software)) meta$ms_control_software <- ms_software
      }, error = function(e) {
        message("[instrument_meta] HyStarMetadata.xml parse failed: ", e$message)
      })
    }

    # --- DIA window info from .m/diaSettings.diasqlite ---
    method_dirs <- list.dirs(tdf_dir, recursive = FALSE)
    m_dirs <- method_dirs[grepl("[.]m$", method_dirs)]
    if (length(m_dirs) > 0) {
      dia_db_path <- file.path(m_dirs[1], "diaSettings.diasqlite")
      if (file.exists(dia_db_path)) {
        tryCatch({
          # Copy to temp dir for SQLite compatibility (SMB mounts can't do WAL/sync)
          tmp_dia <- tempfile(fileext = ".sqlite")
          file.copy(dia_db_path, tmp_dia)
          dia_db <- DBI::dbConnect(RSQLite::SQLite(), tmp_dia, flags = RSQLite::SQLITE_RO)
          on.exit({ DBI::dbDisconnect(dia_db); unlink(tmp_dia) }, add = TRUE)
          if ("DiaWindowsSpecification" %in% DBI::dbListTables(dia_db)) {
            dia_summary <- DBI::dbGetQuery(dia_db, "
              SELECT COUNT(*) AS n_windows,
                     MIN(IsolationMz) AS min_mz,
                     MAX(IsolationMz) AS max_mz,
                     MIN(IsolationWidth) AS isolation_width
              FROM DiaWindowsSpecification")
            if (dia_summary$n_windows > 0) {
              meta$dia_windows <- dia_summary$n_windows
              meta$dia_mz_range <- sprintf("%.0f-%.0f", dia_summary$min_mz, dia_summary$max_mz)
              meta$dia_isolation_width <- dia_summary$isolation_width
            }
          }
        }, error = function(e) NULL)
      }
    }

    # --- Fallback: LC method from submethods dir ---
    if (is.null(meta$lc_method)) {
      submethods_dir <- file.path(tdf_dir, "submethods")
      if (dir.exists(submethods_dir) && requireNamespace("xml2", quietly = TRUE)) {
        method_files <- list.files(submethods_dir, pattern = "\\.method$", full.names = TRUE)
        if (length(method_files) > 0) {
          tryCatch({
            sxml <- xml2::read_xml(method_files[1])
            method_name <- xml2::xml_text(xml2::xml_find_first(sxml, "//method_name"))
            if (!is.na(method_name) && nzchar(method_name)) {
              meta$lc_method <- method_name
            } else {
              meta$lc_method <- tools::file_path_sans_ext(basename(method_files[1]))
            }
          }, error = function(e) NULL)
        }
      }
    }

    meta
  }, error = function(e) {
    list(instrument_type = "timsTOF", parse_error = e$message)
  })
}

#' Parse timsTOF metadata from a .d directory
#' @param d_path Path to a .d directory containing analysis.tdf
#' @return Named list of instrument metadata, or NULL on failure
parse_timstof_metadata <- function(d_path) {
  tdf_path <- file.path(d_path, "analysis.tdf")
  if (!file.exists(tdf_path)) return(NULL)
  parse_timstof_from_tdf(tdf_path)
}

#' Parse Thermo .raw file metadata via rawrr
#' @param raw_path Path to a .raw file
#' @return Named list of instrument metadata, or NULL on failure
parse_thermo_metadata <- function(raw_path) {
  if (!file.exists(raw_path)) return(NULL)

  meta <- list(instrument_type = "Thermo")

  if (!requireNamespace("rawrr", quietly = TRUE)) {
    meta$parse_error <- "rawrr not available"
    return(meta)
  }

  tryCatch({
    header <- rawrr::readFileHeader(raw_path)

    # Field names from rawrr use spaces: "Instrument model", "Serial number", etc.
    meta$instrument_model    <- header$`Instrument model` %||% NA_character_
    meta$instrument_name     <- header$`Instrument name` %||% NA_character_
    meta$instrument_serial   <- header$`Serial number` %||% NA_character_
    meta$acquisition_date    <- header$`Creation date` %||% NA_character_
    meta$operator            <- header$Operator %||% NA_character_
    meta$sample_name         <- header$`Sample name` %||% NA_character_
    meta$acquisition_software <- header$`Software version` %||% NA_character_
    meta$firmware_version    <- header$`Firmware version` %||% NA_character_

    # Scan counts
    meta$n_scans      <- as.integer(header$`Number of scans` %||% NA)
    meta$n_ms1_spectra <- meta$n_scans - as.integer(header$`Number of ms2 scans` %||% 0)
    meta$n_ms2_spectra <- as.integer(header$`Number of ms2 scans` %||% NA)

    # Mass range — "Mass range" may be a string like "100.0000-2000.0000"
    mass_range <- header$`Mass range`
    if (!is.null(mass_range) && is.character(mass_range)) {
      mr_parts <- strsplit(mass_range, "-")[[1]]
      if (length(mr_parts) == 2) {
        meta$mz_range_low  <- as.numeric(trimws(mr_parts[1]))
        meta$mz_range_high <- as.numeric(trimws(mr_parts[2]))
      }
    }

    # Time range — "Time range" may be a string like "0.00-60.00"
    time_range <- header$`Time range`
    if (!is.null(time_range) && is.character(time_range)) {
      tr_parts <- strsplit(time_range, "-")[[1]]
      if (length(tr_parts) == 2) {
        meta$rt_start_min <- round(as.numeric(trimws(tr_parts[1])), 1)
        meta$rt_end_min   <- round(as.numeric(trimws(tr_parts[2])), 1)
      }
    }

    # Fallback: use readIndex if Time range not available
    if (is.na(meta$rt_end_min %||% NA)) {
      tryCatch({
        idx <- rawrr::readIndex(raw_path)
        if (!is.null(idx) && nrow(idx) > 0) {
          meta$rt_start_min <- round(min(idx$StartTime, na.rm = TRUE), 1)
          meta$rt_end_min   <- round(max(idx$StartTime, na.rm = TRUE), 1)
        }
      }, error = function(e) NULL)
    }

    # Mass resolution
    meta$mass_resolution <- header$`Mass resolution` %||% NA_character_

    # Scan filter analysis — extract acquisition mode and DIA info
    # Thermo scan filters look like: "FTMS + p NSI Full ms [400.0000-1200.0000]"
    # DIA MS2: "FTMS + p NSI Full ms2 500.0000-525.0000@hcd30.00 [100.0000-1500.0000]"
    first_filter <- header$`Scan filter (first scan)` %||% NA_character_
    last_filter <- header$`Scan filter (last scan)` %||% NA_character_
    if (!is.na(first_filter)) {
      meta$scan_filter_first <- first_filter
      # Detect mass analyzer from filter
      if (grepl("FTMS", first_filter)) {
        meta$ms1_mass_analyzer <- "Orbitrap"
        meta$ms2_mass_analyzer <- "Orbitrap"  # May be overridden below
      } else if (grepl("ITMS", first_filter)) {
        meta$ms1_mass_analyzer <- "Ion Trap"
      }
    }
    if (!is.na(last_filter)) {
      meta$scan_filter_last <- last_filter
      if (grepl("ITMS", last_filter) && !grepl("ITMS", first_filter %||% "")) {
        meta$ms2_mass_analyzer <- "Ion Trap"  # Hybrid: Orbitrap MS1, IT MS2
      }
    }

    # Try to determine DIA vs DDA from scan index (if available)
    tryCatch({
      idx <- rawrr::readIndex(raw_path)
      if (!is.null(idx) && nrow(idx) > 0) {
        ms2_scans <- idx[idx$MSOrder == "Ms2", ]
        if (nrow(ms2_scans) > 0) {
          # DIA: MS2 scans typically have 0 charge and no precursor selection
          # DDA: MS2 scans have specific precursor masses
          n_zero_charge <- sum(ms2_scans$charge == 0, na.rm = TRUE)
          if (n_zero_charge > nrow(ms2_scans) * 0.9) {
            meta$acquisition_mode <- "DIA"
          } else {
            meta$acquisition_mode <- "DDA"
          }
        }
      }
    }, error = function(e) NULL)

    # Note: LC method is not accessible via rawrr (RawFileReader API limitation)
    # Users should manually specify LC conditions for Thermo data

    meta
  }, error = function(e) {
    meta$parse_error <- e$message
    meta
  })
}

#' Dispatch metadata extraction based on file extension
#' @param file_path Path to a raw MS file (.d directory or .raw file)
#' @return Named list of instrument metadata, or NULL on failure
parse_raw_file_metadata <- function(file_path) {
  if (!file.exists(file_path)) return(NULL)

  ext <- tolower(tools::file_ext(file_path))
  # .d directories have no extension — check if it's a directory ending in .d
  if (dir.exists(file_path) && grepl("\\.d$", file_path, ignore.case = TRUE)) {
    return(parse_timstof_metadata(file_path))
  }

  switch(ext,
    "raw" = parse_thermo_metadata(file_path),
    "d"   = parse_timstof_metadata(file_path),
    NULL  # unsupported format
  )
}

#' Format instrument metadata as multi-line methodology text
#' @param meta Named list from parse_*_metadata()
#' @return Character string for Methods tab
format_instrument_methods_text <- function(meta) {
  if (is.null(meta)) return("")

  lines <- c()

  # Instrument identification
  model <- meta$instrument_model %||% meta$instrument_type
  if (!is.na(model) && nzchar(model)) {
    line <- paste0("Data were acquired on a ", model)
    if (!is.null(meta$instrument_serial) && !is.na(meta$instrument_serial)) {
      line <- paste0(line, " (S/N: ", meta$instrument_serial, ")")
    }
    line <- paste0(line, ".")
    lines <- c(lines, line)
  }

  # Acquisition software
  if (!is.null(meta$acquisition_software) && nzchar(trimws(meta$acquisition_software))) {
    lines <- c(lines, paste0("Acquisition software: ", trimws(meta$acquisition_software), "."))
  }

  # Scan range
  if (!is.na(meta$mz_range_low %||% NA) && !is.na(meta$mz_range_high %||% NA)) {
    lines <- c(lines, sprintf("Precursor m/z scan range: %.0f\u2013%.0f.",
                               meta$mz_range_low, meta$mz_range_high))
  }

  # Ion mobility (timsTOF)
  if (!is.na(meta$im_range_low %||% NA) && !is.na(meta$im_range_high %||% NA)) {
    lines <- c(lines, sprintf("Ion mobility (1/K0) range: %.2f\u2013%.2f Vs/cm\u00B2.",
                               meta$im_range_low, meta$im_range_high))
  }

  # Mass analyzer (if different for MS1/MS2, e.g., hybrid instruments)
  ms1_ana <- meta$ms1_mass_analyzer
  ms2_ana <- meta$ms2_mass_analyzer
  if (!is.null(ms1_ana) && !is.null(ms2_ana) && !is.na(ms1_ana) && !is.na(ms2_ana)) {
    if (ms1_ana == ms2_ana) {
      lines <- c(lines, paste0("Mass analyzer: ", ms1_ana, "."))
    } else {
      lines <- c(lines, sprintf("Mass analyzers: MS1 %s, MS2 %s.", ms1_ana, ms2_ana))
    }
  }

  # Mass resolution (Orbitrap)
  if (!is.null(meta$mass_resolution) && !is.na(meta$mass_resolution) && nzchar(meta$mass_resolution)) {
    lines <- c(lines, paste0("Mass resolution: ", meta$mass_resolution, "."))
  }

  # Acquisition mode
  if (!is.null(meta$acquisition_mode) && meta$acquisition_mode != "unknown") {
    lines <- c(lines, paste0("Acquisition mode: ", meta$acquisition_mode, "."))
  }

  # DIA window scheme
  if (!is.null(meta$dia_windows) && !is.na(meta$dia_windows)) {
    dia_line <- sprintf("%d DIA windows", meta$dia_windows)
    if (!is.null(meta$dia_isolation_width) && !is.na(meta$dia_isolation_width))
      dia_line <- paste0(dia_line, sprintf(" (%.0f Da isolation width)", meta$dia_isolation_width))
    if (!is.null(meta$dia_mz_range) && nzchar(meta$dia_mz_range))
      dia_line <- paste0(dia_line, sprintf(", covering m/z %s", meta$dia_mz_range))
    lines <- c(lines, paste0(dia_line, "."))
  }

  # LC system and method
  lc_parts <- c()
  if (!is.null(meta$lc_system) && nzchar(meta$lc_system))
    lc_parts <- c(lc_parts, meta$lc_system)
  if (!is.null(meta$lc_method) && nzchar(meta$lc_method))
    lc_parts <- c(lc_parts, sprintf("method: %s", meta$lc_method))
  if (!is.null(meta$lc_runtime_min) && !is.na(meta$lc_runtime_min))
    lc_parts <- c(lc_parts, sprintf("%.0f min run time", meta$lc_runtime_min))
  # EvoSep SPD info (only add gradient if not already in method name)
  if (!is.null(meta$evosep_spd) && !is.na(meta$evosep_spd) &&
      !is.null(meta$evosep_gradient_min) && !is.na(meta$evosep_gradient_min)) {
    # Only add standalone SPD line if method name doesn't already contain SPD info
    method_has_spd <- any(grepl("SPD|samples per day", lc_parts, ignore.case = TRUE))
    if (!method_has_spd) {
      lc_parts <- c(lc_parts, sprintf("%d samples per day (%.1f min gradient)",
                                       meta$evosep_spd, meta$evosep_gradient_min))
    } else {
      # Just add gradient length
      lc_parts <- c(lc_parts, sprintf("%.1f min gradient", meta$evosep_gradient_min))
    }
  }
  if (length(lc_parts) > 0) {
    lines <- c(lines, paste0("LC separation: ", paste(lc_parts, collapse = ", "), "."))
  }

  # Gradient length (from Frames RT range — total acquisition time)
  if (!is.na(meta$rt_end_min %||% NA)) {
    lines <- c(lines, sprintf("Total acquisition time: %.0f minutes.", meta$rt_end_min))
  }

  # Spectra counts
  n_ms1 <- meta$n_ms1_spectra %||% meta$n_scans
  n_ms2 <- meta$n_ms2_spectra
  if (!is.null(n_ms1) && !is.na(n_ms1) && !is.null(n_ms2) && !is.na(n_ms2)) {
    lines <- c(lines, sprintf("MS1 spectra: %s, MS2 spectra: %s.",
                               format(n_ms1, big.mark = ","), format(n_ms2, big.mark = ",")))
  }

  # Cycle time
  if (!is.null(meta$cycle_time_sec) && !is.na(meta$cycle_time_sec)) {
    lines <- c(lines, sprintf("Cycle time: %.3f s.", meta$cycle_time_sec))
  }

  if (length(lines) == 0) return("")
  paste(lines, collapse = " ")
}

#' Format instrument metadata as a single-paragraph prompt for Claude/AI export
#' @param meta Named list from parse_*_metadata()
#' @return Character string
format_instrument_for_prompt <- function(meta) {
  if (is.null(meta)) return("")

  parts <- c()
  model <- meta$instrument_model %||% meta$instrument_type
  if (!is.na(model) && nzchar(model)) parts <- c(parts, paste0("Instrument: ", model))
  if (!is.na(meta$mz_range_low %||% NA) && !is.na(meta$mz_range_high %||% NA))
    parts <- c(parts, sprintf("m/z range: %.0f-%.0f", meta$mz_range_low, meta$mz_range_high))
  if (!is.na(meta$im_range_low %||% NA) && !is.na(meta$im_range_high %||% NA))
    parts <- c(parts, sprintf("1/K0: %.2f-%.2f", meta$im_range_low, meta$im_range_high))
  if (!is.null(meta$acquisition_mode) && meta$acquisition_mode != "unknown")
    parts <- c(parts, paste0("Mode: ", meta$acquisition_mode))
  if (!is.null(meta$dia_windows) && !is.na(meta$dia_windows)) {
    dia_str <- sprintf("%d DIA windows", meta$dia_windows)
    if (!is.null(meta$dia_isolation_width) && !is.na(meta$dia_isolation_width))
      dia_str <- paste0(dia_str, sprintf(" (%.0f Da)", meta$dia_isolation_width))
    parts <- c(parts, dia_str)
  }
  if (!is.null(meta$lc_system) && nzchar(meta$lc_system)) {
    lc_str <- paste0("LC: ", meta$lc_system)
    if (!is.null(meta$lc_method) && nzchar(meta$lc_method))
      lc_str <- paste0(lc_str, " (", meta$lc_method, ")")
    parts <- c(parts, lc_str)
  }
  if (!is.null(meta$lc_runtime_min) && !is.na(meta$lc_runtime_min))
    parts <- c(parts, sprintf("LC run: %.0f min", meta$lc_runtime_min))
  else if (!is.na(meta$rt_end_min %||% NA))
    parts <- c(parts, sprintf("Acq time: %.0f min", meta$rt_end_min))
  n_ms1 <- meta$n_ms1_spectra %||% meta$n_scans
  n_ms2 <- meta$n_ms2_spectra
  if (!is.null(n_ms1) && !is.na(n_ms1) && !is.null(n_ms2) && !is.na(n_ms2))
    parts <- c(parts, sprintf("MS1: %s, MS2: %s spectra",
                               format(n_ms1, big.mark = ","), format(n_ms2, big.mark = ",")))
  if (!is.null(meta$cycle_time_sec) && !is.na(meta$cycle_time_sec))
    parts <- c(parts, sprintf("Cycle: %.3f s", meta$cycle_time_sec))
  if (!is.null(meta$acquisition_software) && nzchar(trimws(meta$acquisition_software)))
    parts <- c(parts, paste0("Software: ", trimws(meta$acquisition_software)))

  if (length(parts) == 0) return("")
  paste(parts, collapse = " | ")
}

# ==============================================================================
#  ThermoRawFileParser support — for remote Thermo .raw metadata on HPC
# ==============================================================================

#' Parse ThermoRawFileParser JSON metadata output
#'
#' ThermoRawFileParser outputs CV-term annotated JSON with sections:
#' FileProperties, InstrumentProperties, MSData, ScanSettings, SampleProperties
#' Each entry is: {accession, cvLabel, name, value}
#'
#' @param json_path Path to the JSON metadata file
#' @return Named list of instrument metadata (same schema as parse_thermo_metadata)
parse_thermorawfileparser_json <- function(json_path) {
  if (!file.exists(json_path)) return(NULL)

  meta <- list(instrument_type = "Thermo")

  tryCatch({
    json_text <- readLines(json_path, warn = FALSE)
    parsed <- jsonlite::fromJSON(paste(json_text, collapse = "\n"), simplifyVector = TRUE)

    # Helper: extract value by accession or name from a CV-term array
    cv_val <- function(section, accession = NULL, name = NULL) {
      if (is.null(section) || nrow(section) == 0) return(NA_character_)
      if (!is.null(accession)) {
        row <- section[section$accession == accession, , drop = FALSE]
        if (nrow(row) > 0) return(as.character(row$value[1]))
      }
      if (!is.null(name)) {
        row <- section[grepl(name, section$name, ignore.case = TRUE), , drop = FALSE]
        if (nrow(row) > 0) return(as.character(row$value[1]))
      }
      NA_character_
    }

    inst <- parsed$InstrumentProperties
    ms <- parsed$MSData
    scan <- parsed$ScanSettings
    sample <- parsed$SampleProperties
    fp <- parsed$FileProperties

    # Instrument identification
    meta$instrument_model <- cv_val(inst, accession = "MS:1000494")  # Thermo Scientific instrument model
    if (is.na(meta$instrument_model))
      meta$instrument_model <- cv_val(inst, accession = "MS:1000496")  # instrument attribute/name
    meta$instrument_serial <- cv_val(inst, accession = "MS:1000529")  # instrument serial number
    meta$acquisition_software <- cv_val(inst, accession = "NCIT:C111093")  # software version
    meta$firmware_version <- cv_val(inst, accession = "AFR:0001259")  # firmware version

    # File info
    meta$acquisition_date <- cv_val(fp, accession = "NCIT:C69199")  # content creation date

    # Scan counts
    n_ms1 <- cv_val(ms, accession = "PRIDE:0000481")  # number of MS1 spectra
    n_ms2 <- cv_val(ms, accession = "PRIDE:0000482")  # number of MS2 spectra
    meta$n_ms1_spectra <- as.integer(n_ms1 %||% NA)
    meta$n_ms2_spectra <- as.integer(n_ms2 %||% NA)
    meta$n_scans <- as.integer(cv_val(scan, accession = "PRIDE:0000478") %||% NA)

    # m/z range
    mz_min <- cv_val(ms, accession = "PRIDE:0000476")  # MS min MZ
    mz_max <- cv_val(ms, accession = "PRIDE:0000477")  # MS max MZ
    if (!is.na(mz_min)) meta$mz_range_low <- as.numeric(mz_min)
    if (!is.na(mz_max)) meta$mz_range_high <- as.numeric(mz_max)

    # RT range
    rt_min <- cv_val(ms, accession = "PRIDE:0000474")  # MS min RT
    rt_max <- cv_val(ms, accession = "PRIDE:0000475")  # MS max RT
    if (!is.na(rt_min)) meta$rt_start_min <- round(as.numeric(rt_min), 1)
    if (!is.na(rt_max)) meta$rt_end_min <- round(as.numeric(rt_max), 1)

    # Mass resolution
    mass_res <- cv_val(scan, accession = "MS:1000011")
    if (!is.na(mass_res) && nzchar(mass_res)) meta$mass_resolution <- mass_res

    # Charge range
    charge_min <- cv_val(ms, accession = "PRIDE:0000472")
    charge_max <- cv_val(ms, accession = "PRIDE:0000473")
    if (!is.na(charge_min) && !is.na(charge_max))
      meta$charge_range <- paste0(charge_min, "-", charge_max)

    # Fragmentation types
    frag_types <- cv_val(scan, name = "fragmentation")
    if (!is.na(frag_types) && nzchar(frag_types)) meta$fragmentation_type <- frag_types

    # MS scan range (from ScanSettings)
    scan_range <- cv_val(scan, accession = "PRIDE:0000479")
    if (!is.na(scan_range) && nzchar(scan_range)) meta$scan_range <- scan_range

    # Sample info
    sample_name <- cv_val(sample, name = "sample name")
    if (!is.na(sample_name) && nzchar(sample_name)) meta$sample_name <- sample_name

    # Detect DIA vs DDA from charge distribution
    # DIA typically has mostly charge 0 MS2 scans; also check if n_ms2 >> n_ms1
    if (!is.na(meta$n_ms1_spectra) && !is.na(meta$n_ms2_spectra) && meta$n_ms1_spectra > 0) {
      ms2_ratio <- meta$n_ms2_spectra / meta$n_ms1_spectra
      # Heuristic: DIA typically has consistent MS2/MS1 ratio
      if (ms2_ratio > 5) {
        meta$acquisition_mode <- "DIA"
      } else {
        meta$acquisition_mode <- "DDA"
      }
    }

    # Detect mass analyzer from instrument model name
    model <- meta$instrument_model %||% ""
    if (grepl("Orbitrap|Exploris|Eclipse|Astral|Lumos|QE|Q Exactive|Tribrid", model, ignore.case = TRUE)) {
      meta$ms1_mass_analyzer <- "Orbitrap"
      meta$ms2_mass_analyzer <- "Orbitrap"
      # Tribrid/Eclipse/Lumos can use IT for MS2
      if (grepl("Tribrid|Eclipse|Lumos|Fusion", model, ignore.case = TRUE)) {
        meta$ms2_mass_analyzer <- "Orbitrap"  # default; may be IT depending on method
      }
    }

    # Note: ThermoRawFileParser does NOT extract LC method/gradient info
    # (confirmed: https://github.com/compomics/ThermoRawFileParser/issues/86)

    meta
  }, error = function(e) {
    meta$parse_error <- e$message
    meta
  })
}

#' Run ThermoRawFileParser remotely via SSH to extract metadata from a .raw file
#'
#' Requires ThermoRawFileParser to be installed on the remote system
#' (available as BioContainers: biocontainers/thermorawfileparser)
#'
#' @param cfg SSH config list (host, user, key_path, port, etc.)
#' @param remote_raw_path Full path to the .raw file on the remote system
#' @param ssh_exec_fn The ssh_exec function to use (default: ssh_exec from helpers_search.R)
#' @param scp_download_fn The scp_download function to use
#' @return Named list of instrument metadata, or NULL on failure
run_thermorawfileparser_ssh <- function(cfg, remote_raw_path, ssh_exec_fn = ssh_exec,
                                        scp_download_fn = scp_download) {
  tryCatch({
    # Create temp output dir on remote
    remote_tmp <- paste0("/tmp/.delimp_trfp_", Sys.getpid())
    basename_raw <- basename(remote_raw_path)
    basename_noext <- tools::file_path_sans_ext(basename_raw)
    remote_json <- file.path(remote_tmp, paste0(basename_noext, "-metadata.json"))

    # Check if ThermoRawFileParser is available on the remote system
    # Try multiple common installation locations (PATH, module, conda)
    check_cmd <- paste(
      sprintf("mkdir -p %s", shQuote(remote_tmp)),
      "&&",
      "(",
      "which ThermoRawFileParser.sh 2>/dev/null",
      "|| which thermorawfileparser 2>/dev/null",
      "|| which ThermoRawFileParser 2>/dev/null",
      "|| echo 'NOT_FOUND'",
      ")"
    )
    check_result <- ssh_exec_fn(cfg, check_cmd)
    stdout_lines <- check_result$stdout[nzchar(trimws(check_result$stdout))]
    trfp_path <- if (length(stdout_lines) > 0) trimws(tail(stdout_lines, 1)) else ""

    # If not in PATH, try with login shell (loads modules from .bash_profile)
    if (grepl("NOT_FOUND", trfp_path) || !nzchar(trfp_path)) {
      check_result2 <- ssh_exec_fn(cfg, paste(
        "which ThermoRawFileParser.sh 2>/dev/null",
        "|| which thermorawfileparser 2>/dev/null",
        "|| echo 'NOT_FOUND'"
      ), login_shell = TRUE)
      stdout_lines2 <- check_result2$stdout[nzchar(trimws(check_result2$stdout))]
      trfp_path <- if (length(stdout_lines2) > 0) trimws(tail(stdout_lines2, 1)) else ""
    }

    if (grepl("NOT_FOUND", trfp_path) || !nzchar(trfp_path)) {
      message("[instrument_meta] ThermoRawFileParser not found on remote system")
      return(NULL)
    }

    # Run ThermoRawFileParser: -i=input.raw -m=0 (JSON metadata) -o=output_dir
    # -f=4 means no spectral output (metadata only), -m=0 means JSON metadata
    trfp_cmd <- sprintf(
      "%s -i=%s -m=0 -o=%s -f=4 2>&1",
      shQuote(trfp_path),
      shQuote(remote_raw_path),
      shQuote(remote_tmp)
    )
    run_result <- ssh_exec_fn(cfg, trfp_cmd, timeout = 120)

    if (run_result$status != 0) {
      message("[instrument_meta] ThermoRawFileParser failed: ",
              paste(run_result$stdout, collapse = "\n"))
      # Cleanup remote temp
      ssh_exec_fn(cfg, sprintf("rm -rf %s", shQuote(remote_tmp)), timeout = 10)
      return(NULL)
    }

    # Download the JSON metadata file
    local_json <- file.path(tempdir(), paste0("trfp_meta_", basename_noext, ".json"))
    dl <- scp_download_fn(cfg, remote_json, local_json)

    # Cleanup remote temp
    ssh_exec_fn(cfg, sprintf("rm -rf %s", shQuote(remote_tmp)), timeout = 10)

    if (dl$status != 0 || !file.exists(local_json)) {
      message("[instrument_meta] Failed to download ThermoRawFileParser JSON")
      return(NULL)
    }

    # Parse the JSON
    meta <- parse_thermorawfileparser_json(local_json)
    unlink(local_json)

    meta
  }, error = function(e) {
    message("[instrument_meta] ThermoRawFileParser SSH failed: ", e$message)
    NULL
  })
}

# ==============================================================================
#  TIC Extraction & Chromatography QC — timsTOF
# ==============================================================================

#' Extract TIC trace from a timsTOF analysis.tdf file
#'
#' Reads MS1 frames (MsMsType = 0) from the Frames table and returns
#' retention time + summed ion current for each frame.
#'
#' @param tdf_path Path to analysis.tdf file (direct path, not .d directory)
#' @return data.frame with columns: rt_sec, rt_min, tic. NULL on failure.
extract_tic_timstof <- function(tdf_path) {
  if (!file.exists(tdf_path)) return(NULL)
  if (!requireNamespace("DBI", quietly = TRUE) ||
      !requireNamespace("RSQLite", quietly = TRUE)) {
    return(NULL)
  }

  tryCatch({
    # Copy TDF to temp for SQLite compatibility (network/SMB mounts can't do WAL/sync)
    tmp_tdf <- tempfile(fileext = ".tdf")
    file.copy(tdf_path, tmp_tdf)
    db <- DBI::dbConnect(RSQLite::SQLite(), tmp_tdf, flags = RSQLite::SQLITE_RO)
    on.exit({ DBI::dbDisconnect(db); unlink(tmp_tdf) }, add = TRUE)

    if (!("Frames" %in% DBI::dbListTables(db))) return(NULL)

    # Detect TIC column — older TDF schemas use different names
    frame_cols <- DBI::dbListFields(db, "Frames")
    tic_col <- if ("SummedIntensities" %in% frame_cols) {
      "SummedIntensities"
    } else if ("AccumulatedIntensity" %in% frame_cols) {
      "AccumulatedIntensity"
    } else if ("MaxIntensity" %in% frame_cols) {
      "MaxIntensity"
    } else {
      # Last resort: find any numeric intensity-like column
      intensity_cols <- frame_cols[grepl("ntensit", frame_cols, ignore.case = TRUE)]
      if (length(intensity_cols) > 0) intensity_cols[1] else NULL
    }

    if (is.null(tic_col)) {
      message("[tic] No intensity column found in Frames table. Columns: ",
              paste(frame_cols, collapse = ", "))
      return(NULL)
    }

    # Also handle older MsMsType values — some TDF versions don't have this column
    ms1_filter <- if ("MsMsType" %in% frame_cols) {
      "WHERE MsMsType = 0"
    } else if ("ScanMode" %in% frame_cols) {
      "WHERE ScanMode = 0"
    } else {
      ""  # no filter — take all frames
    }

    query <- sprintf("SELECT Time AS rt_sec, %s AS tic FROM Frames %s ORDER BY Time",
                      tic_col, ms1_filter)
    tic_df <- DBI::dbGetQuery(db, query)

    if (nrow(tic_df) == 0) return(NULL)

    tic_df$rt_min <- tic_df$rt_sec / 60
    tic_df
  }, error = function(e) {
    message("[tic] Failed to extract TIC from ", tdf_path, ": ", e$message)
    NULL
  })
}

#' Compute per-run TIC shape metrics
#'
#' Derives summary statistics from a single TIC trace for QC assessment:
#' AUC, peak RT, gradient boundaries, baseline ratio, late signal, asymmetry.
#'
#' @param tic_df data.frame from extract_tic_timstof() with rt_min, tic columns
#' @param run_name Character string identifying this run
#' @return Named list of metrics, or NULL if tic_df is invalid
compute_tic_metrics <- function(tic_df, run_name) {
  if (is.null(tic_df) || nrow(tic_df) < 10) {
    return(list(run = run_name, valid = FALSE))
  }

  # Sort by RT
  tic_df <- tic_df[order(tic_df$rt_min), ]

  # Smooth with moving average (window = 2% of total points)
  n <- nrow(tic_df)
  win <- max(3, round(n * 0.02))
  if (win %% 2 == 0) win <- win + 1  # ensure odd window
  smoothed <- stats::filter(tic_df$tic, rep(1/win, win), sides = 2)
  smoothed[is.na(smoothed)] <- tic_df$tic[is.na(smoothed)]  # fill edges

  # Trapezoid AUC
  dt <- diff(tic_df$rt_min)
  avg_tic <- (smoothed[-n] + smoothed[-1]) / 2
  total_auc <- sum(dt * avg_tic, na.rm = TRUE)

  # Peak

  peak_idx <- which.max(smoothed)
  peak_rt <- tic_df$rt_min[peak_idx]
  peak_tic <- smoothed[peak_idx]

  # Gradient boundaries (10% of peak threshold)
  threshold <- peak_tic * 0.10
  above <- which(smoothed >= threshold)
  ramp_rt <- if (length(above) > 0) tic_df$rt_min[min(above)] else tic_df$rt_min[1]
  tail_rt <- if (length(above) > 0) tic_df$rt_min[max(above)] else tic_df$rt_min[n]
  gradient_width <- tail_rt - ramp_rt

  # Baseline ratio: median of lowest 10% intensities vs peak
  sorted_tic <- sort(smoothed)
  baseline_n <- max(1, round(n * 0.10))
  baseline_median <- median(sorted_tic[seq_len(baseline_n)])
  baseline_ratio <- if (peak_tic > 0) baseline_median / peak_tic else 0

  # Late signal ratio: signal in last 20% of gradient vs total
  rt_range <- range(tic_df$rt_min)
  late_start <- rt_range[2] - 0.20 * (rt_range[2] - rt_range[1])
  late_idx <- which(tic_df$rt_min >= late_start)
  late_signal_ratio <- if (total_auc > 0 && length(late_idx) > 1) {
    late_dt <- diff(tic_df$rt_min[late_idx])
    late_avg <- (smoothed[late_idx[-length(late_idx)]] + smoothed[late_idx[-1]]) / 2
    sum(late_dt * late_avg, na.rm = TRUE) / total_auc
  } else {
    0
  }

  # Asymmetry: left-width vs right-width at 50% peak height
  half_max <- peak_tic * 0.50
  left_half <- which(smoothed[seq_len(peak_idx)] >= half_max)
  right_half <- which(smoothed[peak_idx:n] >= half_max)
  left_width <- if (length(left_half) > 0) peak_rt - tic_df$rt_min[min(left_half)] else 0
  right_width <- if (length(right_half) > 0) tic_df$rt_min[peak_idx - 1 + max(right_half)] - peak_rt else 0
  asymmetry <- if (left_width > 0 && right_width > 0) left_width / right_width else 1

  list(
    run = run_name,
    valid = TRUE,
    total_auc = total_auc,
    peak_rt_min = round(peak_rt, 2),
    peak_tic = peak_tic,
    ramp_rt_min = round(ramp_rt, 2),
    tail_rt_min = round(tail_rt, 2),
    gradient_width_min = round(gradient_width, 2),
    baseline_ratio = round(baseline_ratio, 4),
    late_signal_ratio = round(late_signal_ratio, 4),
    asymmetry = round(asymmetry, 3)
  )
}

#' Compute shape similarity of TIC traces vs median trace
#'
#' Interpolates all traces onto a common RT grid and computes Pearson r
#' of each trace against the median trace.
#'
#' @param all_traces_list Named list of data.frames (each with rt_min, tic columns)
#' @return data.frame with columns: run, shape_r. NULL if fewer than 2 traces.
compute_shape_similarity <- function(all_traces_list) {
  valid <- Filter(function(x) !is.null(x) && nrow(x) >= 10, all_traces_list)
  if (length(valid) < 2) {
    # Single file — no comparison possible, return perfect score
    if (length(valid) == 1) {
      return(data.frame(run = names(valid), shape_r = 1.0, stringsAsFactors = FALSE))
    }
    return(NULL)
  }

  # Common RT grid: intersection of all RT ranges, 500 points
  rt_mins <- sapply(valid, function(x) min(x$rt_min))
  rt_maxs <- sapply(valid, function(x) max(x$rt_min))
  common_start <- max(rt_mins)
  common_end <- min(rt_maxs)
  if (common_end <= common_start) return(NULL)

  grid_rt <- seq(common_start, common_end, length.out = 500)

  # Interpolate each trace onto the grid
  interp_matrix <- sapply(valid, function(x) {
    stats::approx(x$rt_min, x$tic, xout = grid_rt, rule = 2)$y
  })

  # Median trace
  median_trace <- apply(interp_matrix, 1, median)

  # Pearson correlation of each trace vs median
  cors <- apply(interp_matrix, 2, function(col) {
    if (sd(col) == 0 || sd(median_trace) == 0) return(0)
    cor(col, median_trace, use = "complete.obs")
  })

  data.frame(
    run = names(cors),
    shape_r = round(unname(cors), 4),
    stringsAsFactors = FALSE
  )
}

#' Diagnose a single run based on its metrics and the cohort
#'
#' Uses MAD-based thresholds for robust outlier detection.
#'
#' @param metrics Named list from compute_tic_metrics() for this run
#' @param all_metrics_df data.frame of all runs' metrics (for cohort comparison)
#' @param shape_r Numeric shape similarity (Pearson r) for this run
#' @return list(status = "pass"/"warn"/"fail", flags = character())
diagnose_run <- function(metrics, all_metrics_df, shape_r = 1.0) {
  flags <- character()
  status <- "pass"

  if (!metrics$valid) {
    return(list(status = "fail", flags = "Invalid TIC extraction"))
  }

  # Shape deviation
  if (!is.na(shape_r)) {
    if (shape_r < 0.90) {
      flags <- c(flags, sprintf("Shape deviation: r=%.3f (fail <0.90)", shape_r))
      status <- "fail"
    } else if (shape_r < 0.95) {
      flags <- c(flags, sprintf("Shape deviation: r=%.3f (warn <0.95)", shape_r))
      if (status != "fail") status <- "warn"
    }
  }

  # Cohort-based checks (need at least 3 runs for MAD)
  valid_metrics <- all_metrics_df[all_metrics_df$valid, ]
  if (nrow(valid_metrics) >= 3) {
    # RT shift: peak RT > 3 MAD from median
    med_peak_rt <- median(valid_metrics$peak_rt_min, na.rm = TRUE)
    mad_peak_rt <- mad(valid_metrics$peak_rt_min, na.rm = TRUE)
    if (mad_peak_rt > 0) {
      rt_dev <- abs(metrics$peak_rt_min - med_peak_rt) / mad_peak_rt
      if (rt_dev > 3) {
        flags <- c(flags, sprintf("RT shift: peak at %.1f min (median %.1f, %.1f MAD)",
                                   metrics$peak_rt_min, med_peak_rt, rt_dev))
        if (status != "fail") status <- "warn"
      }
    }

    # Loading anomaly: AUC vs median
    med_auc <- median(valid_metrics$total_auc, na.rm = TRUE)
    if (med_auc > 0) {
      auc_ratio <- metrics$total_auc / med_auc
      if (auc_ratio > 3 || auc_ratio < 0.3) {
        flags <- c(flags, sprintf("Loading anomaly: AUC %.1fx median (fail)", auc_ratio))
        status <- "fail"
      } else if (auc_ratio > 2 || auc_ratio < 0.5) {
        flags <- c(flags, sprintf("Loading anomaly: AUC %.1fx median (warn)", auc_ratio))
        if (status != "fail") status <- "warn"
      }
    }

    # Narrow gradient
    med_gradient <- median(valid_metrics$gradient_width_min, na.rm = TRUE)
    if (med_gradient > 0 && metrics$gradient_width_min < 0.70 * med_gradient) {
      flags <- c(flags, sprintf("Narrow gradient: %.1f min (median %.1f)",
                                 metrics$gradient_width_min, med_gradient))
      if (status != "fail") status <- "warn"
    }

    # File size outlier
    sizes <- valid_metrics$size_mb
    if (!is.null(sizes) && !all(is.na(sizes))) {
      med_size <- median(sizes, na.rm = TRUE)
      this_size <- metrics$size_mb
      if (!is.null(this_size) && !is.na(this_size) && med_size > 0) {
        size_ratio <- this_size / med_size
        if (size_ratio < 0.3 || size_ratio > 3) {
          flags <- c(flags, sprintf("File size: %d MB (%.1fx median %d MB, fail)",
                                     round(this_size), size_ratio, round(med_size)))
          status <- "fail"
        } else if (size_ratio < 0.5 || size_ratio > 2) {
          flags <- c(flags, sprintf("File size: %d MB (%.1fx median %d MB, warn)",
                                     round(this_size), size_ratio, round(med_size)))
          if (status != "fail") status <- "warn"
        }
      }
    }
  }

  # Absolute thresholds (independent of cohort)
  # Late elution
  if (metrics$late_signal_ratio > 0.15) {
    flags <- c(flags, sprintf("Late elution: %.0f%% signal in last 20%% of gradient",
                               metrics$late_signal_ratio * 100))
    if (status != "fail") status <- "warn"
  }

  # Elevated baseline
  if (metrics$baseline_ratio > 0.10) {
    flags <- c(flags, sprintf("Elevated baseline: %.1f%% of peak",
                               metrics$baseline_ratio * 100))
    if (status != "fail") status <- "warn"
  }

  list(status = status, flags = flags)
}

#' Normalize TIC intensities to 0-1 range for overlay plotting
#'
#' @param tic_df data.frame with tic column
#' @return Same data.frame with added tic_norm column
normalize_tic <- function(tic_df) {
  if (is.null(tic_df) || nrow(tic_df) == 0) return(tic_df)
  tic_range <- range(tic_df$tic, na.rm = TRUE)
  denom <- tic_range[2] - tic_range[1]
  tic_df$tic_norm <- if (denom > 0) (tic_df$tic - tic_range[1]) / denom else 0
  tic_df
}
