# figures.R -- the three figures of the article (notebook 2).
#
#   fig_dataset_prisma()  Figure 1: the 122 series of the GEO query, the inclusion criteria,
#                         the nine series analysed (six from the query, three human series
#                         identified separately) and the series held out for validation.
#   fig_workflow()        Figure 3: the analysis workflow.
#   fig_interventions()   Figure 2: three injury pathways, the signature genes on them and
#                         candidate interventions from the literature (sources in the caption).
#                         It stops if a pathway or gene drawn as part of the signature is no
#                         longer supported by the enrichment table or the signature.
# The numbers drawn in Figures 1 and 3 are passed in from the result tables.

suppressWarnings(suppressMessages({ library(ggplot2) }))

fig_dataset_prisma <- function(n_identified, n_included, file, organisms = NULL,
                               searched_on = NULL, coverage = NULL, query = NULL) {
  stopifnot("the coverage table (results/geo_corpus_coverage.csv) is required" = !is.null(coverage),
            "organisms must be named by accession" = !is.null(names(organisms)))
  common <- c("Homo sapiens" = "human", "Mus musculus" = "mouse", "Rattus norvegicus" = "rat")
  corpus   <- coverage[coverage$role == "corpus", ]
  in_query <- corpus$accession[corpus$returned_by_query]
  by_hand  <- corpus$accession[!corpus$returned_by_query]
  held_out <- coverage$accession[coverage$role == "held-out"]
  stopifnot(length(in_query) + length(by_hand) == n_included, all(corpus$accession %in% names(organisms)))
  DOT <- "  \u00b7  "
  species_of <- function(acc) unname(common[organisms[acc]])
  gse_order  <- function(acc) acc[order(as.numeric(sub("^GSE", "", acc)))]

  INK <- "#22303F"; MUTE <- "#5A6B7B"; RULE <- "#C2CEDA"
  FILL <- "#F5F8FB"; FILL_IN <- "#E9F3EC"; EDGE_IN <- "#4E8C63"
  FILL_SIDE <- "#FCF8F0"; EDGE_SIDE <- "#AD8A52"
  HEAD <- 5.4; LEAD <- 7.0; BODY <- 5.6      # text sizes (mm)

  # One text line: words, size, weight, colour.
  ln <- function(text, size = BODY, face = "plain", col = MUTE) list(text = text, size = size, face = face, col = col)
  query_lines <- c("(kidney OR renal)[Title]",
                   "AND   (ischemi* OR reperfusion)[Title]",
                   "AND   expression profiling by array /",
                   "by high-throughput sequencing",
                   "AND   gse[Filter]")
  boxes <- list(
    ident = list(head = "IDENTIFICATION", col = MUTE, fill = FILL, edge = RULE, lines = c(
      list(ln(sprintf("%s series", format(n_identified, big.mark = ",")), LEAD, "bold", INK),
           ln("returned by the declared Entrez query over GEO")),
      lapply(query_lines, ln),
      if (!is.null(searched_on)) list(ln(sprintf("queried %s", searched_on))) else list())),
    criteria = list(head = "INCLUSION CRITERIA", col = MUTE, fill = FILL, edge = RULE, lines = list(
      ln("renal ischemia-reperfusion injury"),
      ln("and reference groups,"),
      ln("verified sample by sample in the GEO metadata"),
      ln(paste("kidney tissue", "whole-transcriptome platform", sep = DOT)))),
    included = list(head = "INCLUDED", col = EDGE_IN, fill = FILL_IN, edge = EDGE_IN, lines = c(
      list(ln(sprintf("%d series", n_included), LEAD, "bold", INK),
           ln(sprintf("%d from the query", length(in_query)))),
      lapply(intersect(c("mouse", "rat", "human"), species_of(in_query)), function(sp)
        ln(sprintf("%s: %s", sp, paste(gse_order(in_query[species_of(in_query) == sp]), collapse = DOT)))),
      list(ln(sprintf("+ %d %s series identified separately", length(by_hand),
                      paste(unique(species_of(by_hand)), collapse = " and ")))))),
    held = list(head = "HELD OUT", col = MUTE, fill = FILL, edge = RULE, lines = list(
      ln(paste(held_out, collapse = DOT), LEAD, "bold", INK),
      ln("returned by the declared query,"),
      ln("reserved for external evaluation"))),
    separate = list(head = "IDENTIFIED SEPARATELY", col = EDGE_SIDE, fill = FILL_SIDE, edge = EDGE_SIDE, lines = list(
      ln(sprintf("%d human transplant series", length(by_hand)), LEAD, "bold", INK),
      ln("not returned by the declared query"),
      ln(paste(gse_order(by_hand), collapse = DOT)),
      ln("their titles carry neither"),
      ln("\"ischemia\" nor \"reperfusion\""))))

  # Sizes in inches, measured on the text itself.
  pt <- function(size) size * .pt
  line_h <- function(size) pt(size) * 1.42 / 72.27
  text_w <- function(l) {
    grDevices::pdf(NULL); on.exit(grDevices::dev.off())
    g <- grid::textGrob(l$text, gp = grid::gpar(fontsize = pt(l$size), fontface = l$face, fontfamily = "sans"))
    grid::convertWidth(grid::grobWidth(g), "in", valueOnly = TRUE)
  }
  head_line <- function(b) ln(b$head, HEAD, "bold", b$col)
  PADX <- 0.42; PADY <- 0.24
  box_w <- function(b) max(vapply(c(list(head_line(b)), b$lines), text_w, numeric(1))) + 2 * PADX
  box_h <- function(b) sum(vapply(c(list(head_line(b)), b$lines), function(l) line_h(l$size), numeric(1))) +
                       0.10 + 2 * PADY
  wl <- max(sapply(boxes[c("ident", "criteria", "included")], box_w))
  wr <- max(sapply(boxes[c("held", "separate")], box_w))
  h  <- sapply(boxes, box_h)
  GAPX <- 0.95; GAPY <- 0.62
  xl <- 0; xr <- wl / 2 + GAPX + wr / 2

  # Rows: identification | criteria + held out | included + identified separately.
  row_h <- c(h["ident"], max(h["criteria"], h["held"]), max(h["included"], h["separate"]))
  yc <- numeric(3); top <- 0
  for (i in 1:3) { yc[i] <- top - row_h[i] / 2; top <- top - row_h[i] - GAPY }
  pos <- list(ident = c(xl, yc[1], wl), criteria = c(xl, yc[2], wl), included = c(xl, yc[3], wl),
              held = c(xr, yc[2], wr), separate = c(xr, yc[3], wr))

  h[c("criteria", "held")] <- row_h[2]; h[c("included", "separate")] <- row_h[3]
  rects <- do.call(rbind, lapply(names(boxes), function(k) data.frame(
    xmin = pos[[k]][1] - pos[[k]][3] / 2, xmax = pos[[k]][1] + pos[[k]][3] / 2,
    ymin = pos[[k]][2] - h[[k]] / 2, ymax = pos[[k]][2] + h[[k]] / 2,
    fill = boxes[[k]]$fill, edge = boxes[[k]]$edge, stringsAsFactors = FALSE)))
  texts <- do.call(rbind, lapply(names(boxes), function(k) {
    b <- boxes[[k]]; ls <- c(list(head_line(b)), b$lines)
    hh <- vapply(ls, function(l) line_h(l$size), numeric(1)); hh[1] <- hh[1] + 0.10
    y_top <- pos[[k]][2] + sum(hh) / 2
    data.frame(x = pos[[k]][1], y = y_top - cumsum(hh) + hh / 2,
               label = vapply(ls, `[[`, "", "text"), size = vapply(ls, `[[`, 0, "size"),
               face = vapply(ls, `[[`, "", "face"), col = vapply(ls, `[[`, "", "col"),
               stringsAsFactors = FALSE)
  }))
  bottom <- function(k) pos[[k]][2] - h[[k]] / 2
  top_of <- function(k) pos[[k]][2] + h[[k]] / 2
  down <- data.frame(x = xl, y = c(bottom("ident"), bottom("criteria")),
                     yend = c(top_of("criteria"), top_of("included")))
  ARR <- arrow(length = unit(0.26, "cm"), type = "closed")

  p <- ggplot() +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = rects$fill, color = rects$edge, linewidth = 0.8) +
    geom_segment(data = down, aes(x = x, xend = x, y = y, yend = yend),
                 arrow = ARR, color = INK, linewidth = 0.7) +
    # criteria -> held out: set aside (dashed)
    annotate("segment", x = xl + wl / 2, xend = xr - wr / 2, y = yc[2], yend = yc[2],
             arrow = ARR, color = MUTE, linewidth = 0.6, linetype = "22") +
    # identified separately -> included (solid)
    annotate("segment", x = xr - wr / 2, xend = xl + wl / 2, y = yc[3], yend = yc[3],
             arrow = ARR, color = INK, linewidth = 0.7) +
    geom_text(data = texts, aes(x = x, y = y, label = label), size = texts$size,
              fontface = texts$face, color = texts$col, family = "sans")

  M <- 0.30
  xlim <- c(xl - wl / 2 - M, xr + wr / 2 + M)
  ylim <- c(min(rects$ymin) - M, max(rects$ymax) + M)
  p <- p + coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0),
          plot.background = element_rect(fill = "white", color = NA))
  ggsave(file, p, width = diff(xlim), height = diff(ylim), dpi = 300)
  p
}

fig_workflow <- function(file, n_sig = NULL, n_terms = NULL, n_path = NULL, n_drug = NULL,
                         n_datasets = NULL, composition = NULL) {
  n <- function(x, alt) if (is.null(x)) alt else as.character(x)

  INK <- "#22303F"; MUTE <- "#5A6B7B"; RULE <- "#C2CEDA"
  FILL <- "#F5F8FB"; FILL_SIG <- "#E9F3EC"; EDGE_SIG <- "#4E8C63"
  FILL_ORA <- "#FCF6E6"; EDGE_ORA <- "#B8922F"

  HEAD_SIZE <- 6.0
  LEAD_SIZE <- 7.4
  BODY_SIZE <- 6.0

  boxes <- list(
    list(head = "DATASETS", head_col = MUTE, fill = FILL, edge = RULE,
         lead = sprintf("%s public renal-IRI series", n(n_datasets, "9")), lead_size = LEAD_SIZE,
         body = n(composition, "3 human   |   3 mouse   |   3 rat"), body_size = BODY_SIZE),
    list(head = "DIFFERENTIAL EXPRESSION", head_col = MUTE, fill = FILL, edge = RULE,

         lead = "within each dataset: injured tissue versus reference group",
         lead_size = LEAD_SIZE,
         body = "limma   |   DESeq2", body_size = BODY_SIZE),
    list(head = "RANK AGGREGATION", head_col = MUTE, fill = FILL, edge = RULE,
         lead = "RobustRankAggreg", lead_size = LEAD_SIZE,
         body = "rodent genes mapped to their human orthologs", body_size = BODY_SIZE),
    list(head = "CONSENSUS SIGNATURE", head_col = EDGE_SIG, fill = FILL_SIG, edge = EDGE_SIG,
         lead = sprintf("%s genes", n(n_sig, "n")), lead_size = LEAD_SIZE,
         body = "measured in \u2265 5 of 9 series   |   adjusted p < 0.05   |   same direction in \u2265 half",
         body_size = BODY_SIZE),
    list(head = "ENRICHMENT", head_col = EDGE_ORA, fill = FILL_ORA, edge = EDGE_ORA,
         lead = sprintf("%s enriched terms", n(n_terms, "n")), lead_size = LEAD_SIZE,
         body = c(sprintf("%s pathways   |   %s drug gene sets", n(n_path, "n"), n(n_drug, "n")),
                  "WebGestalt over KEGG, Reactome, DrugBank and GLAD4U"),
         body_size = BODY_SIZE))

  LHF <- 1.30
  line_in <- function(size, lh = LHF) size * .pt * lh / 72.27
  PAD <- 0.17
  hgt <- function(b) max(0.52, (line_in(HEAD_SIZE, 1.55) + line_in(b$lead_size, 1.70) +
                                length(b$body) * line_in(b$body_size)) / 2 + PAD)
  bh <- vapply(boxes, hgt, numeric(1))
  GAP <- 0.46
  yc <- numeric(length(boxes)); top <- 0
  for (i in seq_along(boxes)) { yc[i] <- top - bh[i]; top <- yc[i] - bh[i] - GAP }
  bw <- 4.85

  p <- ggplot() +
    geom_segment(data = data.frame(y = (yc - bh)[-length(yc)], yend = (yc + bh)[-1]),
                 aes(x = 0, xend = 0, y = y, yend = yend),
                 arrow = arrow(length = unit(0.24, "cm"), type = "closed"),
                 color = INK, linewidth = 0.6) +
    geom_rect(data = data.frame(y = yc, h = bh),
              aes(xmin = -bw, xmax = bw, ymin = y - h, ymax = y + h),
              fill = vapply(boxes, function(b) b$fill, character(1)),
              color = vapply(boxes, function(b) b$edge, character(1)), linewidth = 0.7)

  for (i in seq_along(boxes)) {
    b <- boxes[[i]]
    cur <- yc[i] + bh[i] - PAD - line_in(HEAD_SIZE) / 2
    p <- p + annotate("text", x = 0, y = cur, label = b$head, size = HEAD_SIZE,
                      fontface = "bold", color = b$head_col, family = "sans")
    cur <- cur - line_in(HEAD_SIZE, 1.55)
    p <- p + annotate("text", x = 0, y = cur, label = b$lead, size = b$lead_size,
                      fontface = "bold", color = INK, family = "sans")
    cur <- cur - line_in(b$lead_size, 1.70)
    p <- p + annotate("text", x = 0, y = cur + line_in(b$body_size) / 2,
                      label = paste(b$body, collapse = "\n"), size = b$body_size,
                      color = MUTE, family = "sans", vjust = 1, lineheight = LHF)
  }

  ylim <- c(min(yc - bh) - 0.20, max(yc + bh) + 0.20)

  TOP_MATTER <- 1.10
  p <- p +
    coord_cartesian(xlim = c(-bw - 0.25, bw + 0.25), ylim = ylim) +
    theme_void(base_size = 13) +
    theme(plot.margin = margin(22, 24, 20, 24),
          plot.background = element_rect(fill = "white", color = NA))
  ggsave(file, p, width = 10.9, height = diff(ylim) + TOP_MATTER, dpi = 300)
  p
}

fig_interventions <- function(file, ora, sig) {
  if ("signature" %in% names(sig)) sig <- sig[sig$signature %in% c(TRUE, "TRUE"), ]
  terms <- c("TNF signaling pathway", "NF-kappa B signaling pathway", "IL-17 signaling pathway")
  stopifnot("an enriched term drawn in green is missing from the ORA table" =
              all(terms %in% ora$description))
  stopifnot("CXCL2 is no longer an up-regulated signature gene" =
              identical(sig$direction[sig$human == "CXCL2"], "up"))
  stopifnot("CXCL2 is not in all three terms drawn upstream of it" =
              all(vapply(terms, function(t) "CXCL2" %in%
                           strsplit(ora$overlap_genes[ora$description == t][1], ";")[[1]], logical(1))))

  INK <- "#22303F"; MUTE <- "#5A6B7B"; LINE <- "#6E7F90"
  GREEN <- "#2E8B57"; GREEN_FILL <- "#E3F1E9"
  NEUT  <- "#9FB0C0"; NEUT_FILL  <- "#F1F4F7"
  GOLD  <- "#B8922F"; GOLD_FILL  <- "#FBF2DA"; GOLD_TXT <- "#6B5214"
  END_FILL <- "#2F3E4E"
  PT <- function(pt) pt / ggplot2::.pt

  rr <- function(x0, x1, y0, y1, r, n = 12) {
    a <- seq(0, pi / 2, length.out = n)
    m <- rbind(cbind(x1 - r + r * cos(a), y1 - r + r * sin(a)),
               cbind(x0 + r - r * sin(a), y1 - r + r * cos(a)),
               cbind(x0 + r - r * cos(a), y0 + r - r * sin(a)),
               cbind(x1 - r + r * sin(a), y0 + r - r * cos(a)))
    data.frame(x = m[, 1], y = m[, 2])
  }
  STY <- list(sig = list(fill = GREEN_FILL, col = GREEN, lty = "solid", txt = INK),
              lit = list(fill = NEUT_FILL,  col = NEUT,  lty = "solid", txt = INK),
              int = list(fill = GOLD_FILL,  col = GOLD,  lty = "22",    txt = GOLD_TXT),
              end = list(fill = END_FILL,   col = END_FILL, lty = "solid", txt = "white"))
  node <- function(cx, cy, w, h, label, style, size = 22, parse = FALSE, face = "plain", tx = cx) {
    s <- STY[[style]]
    list(geom_polygon(data = rr(cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2, min(h / 2, 0.22)),
                      aes(x, y), fill = s$fill, colour = s$col, linetype = s$lty, linewidth = 0.9),
         annotate("text", x = tx, y = cy, label = label, parse = parse, size = PT(size),
                  colour = s$txt, fontface = face, family = "sans", lineheight = 0.92))
  }
  arr <- function(x0, y0, x1, y1, col = LINE, lty = "solid", lw = 0.9)
    annotate("segment", x = x0, y = y0, xend = x1, yend = y1, colour = col, linetype = lty,
             linewidth = lw, arrow = arrow(length = unit(0.14, "in"), type = "closed"))
  tbar <- function(x, y0, y1, col = GOLD) list(
    annotate("segment", x = x, xend = x, y = y0, yend = y1, colour = col, linewidth = 1.0,
             linetype = "22"),
    annotate("segment", x = x - 0.17, xend = x + 0.17, y = y1, yend = y1, colour = col,
             linewidth = 1.8))
  down <- function(x, y) annotate("segment", x = x, xend = x, y = y + 0.14, yend = y - 0.14,
                                  colour = INK, linewidth = 0.9,
                                  arrow = arrow(length = unit(0.08, "in"), type = "closed"))

  H <- 0.60; L1 <- 5.15; L2 <- 3.75; L3 <- 2.35
  TOPI <- 6.62; CITY <- 3.05; NAMY <- 1.15; LEGY <- 0.40; EX0 <- 10.10; EX1 <- 11.80

  p <- ggplot() + coord_fixed(xlim = c(0, 12.2), ylim = c(0, 7.12), expand = FALSE) +
    theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA),
                         plot.margin = margin(6, 6, 6, 6))

  p <- p +
    node(1.35, L1 + 0.42, 1.50, H - 0.04, "IL-17", "sig") +
    node(1.35, L1 - 0.42, 1.50, H - 0.04, "TNF", "sig") +
    arr(2.10, L1 + 0.42, 3.00, L1 + 0.10) + arr(2.10, L1 - 0.42, 3.00, L1 - 0.10) +
    node(3.85, L1, 1.70, H, 'NF*"-"*kappa*B', "sig", parse = TRUE) +
    arr(4.70, L1, 5.40, L1) +
    node(6.10, L1, 1.40, H, "CXCL2", "sig") +
    arr(6.80, L1, 7.50, L1) +
    node(8.55, L1, 2.10, H, "Neutrophils", "lit") +
    arr(9.60, L1, EX0, L1) +
    node(1.35, TOPI, 2.20, 0.58, "Anti-IL-17A", "int", 21, face = "bold") +
    tbar(1.35, TOPI - 0.29, L1 + 0.42 + (H - 0.04) / 2 + 0.03) +
    node(7.15, TOPI, 2.00, 0.58, "Anti-CXCR2", "int", 21, face = "bold") +
    tbar(7.15, TOPI - 0.29, L1 + 0.10)

  p <- p +
    node(1.45, L2, 1.70, H, "Succinate", "lit") +
    arr(2.30, L2, 3.70, L2) +
    annotate("text", x = 3.00, y = L2 + 0.25, label = "reperfusion", size = PT(15),
             colour = MUTE, fontface = "italic", family = "sans") +
    node(5.20, L2, 3.00, H, "Mitochondrial ROS", "lit") +
    arr(6.70, L2, EX0, L2) +
    node(8.40, CITY, 1.60, 0.55, "Citrate", "int", 21, face = "bold") +
    tbar(8.40, CITY + 0.275, L2 - 0.10)

  p <- p +
    node(1.45, L3, 1.70, H, "QPRT", "lit", tx = 1.27) + down(1.95, L3) +
    arr(2.30, L3, 3.20, L3) +
    node(3.95, L3, 1.50, H, 'NAD^"+"', "lit", parse = TRUE, tx = 3.78) + down(4.38, L3) +
    arr(4.70, L3, EX0, L3) +
    node(3.95, NAMY, 2.30, 0.55, "Nicotinamide", "int", 21, face = "bold") +
    arr(3.95, NAMY + 0.275, 3.95, L3 - H / 2 - 0.05, col = GOLD, lty = "22", lw = 1.0)

  p <- p + node(mean(c(EX0, EX1)), mean(c(L1 + 0.45, L3 - 0.45)), EX1 - EX0,
                (L1 + 0.45) - (L3 - 0.45), "Kidney\ninjury", "end", 24, face = "bold")

  key <- function(x, style, label) list(
    geom_polygon(data = rr(x, x + 0.55, LEGY - 0.15, LEGY + 0.15, 0.15), aes(x, y),
                 fill = STY[[style]]$fill, colour = STY[[style]]$col,
                 linetype = STY[[style]]$lty, linewidth = 0.8),
    annotate("text", x = x + 0.72, y = LEGY, label = label, hjust = 0, size = PT(16),
             colour = INK, family = "sans"))
  p <- p + key(1.55, "sig", "From the signature") + key(4.45, "lit", "From the literature") +
    key(7.55, "int", "Candidate intervention")

  ggsave(file, p, width = 12.2, height = 7.12, dpi = 300, device = ragg::agg_png, bg = "white")
  p
}
