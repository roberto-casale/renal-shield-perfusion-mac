# figures.R -- the three figures of the article (notebook 2).
#
#   fig_dataset_prisma()  Figure 1: identification and inclusion of the nine series.
#   fig_workflow()        Figure 3: the analysis workflow.
#   fig_interventions()   Figure 2: three injury pathways, the signature genes on them and
#                         candidate interventions from the literature (sources in the caption).
#                         It stops if a pathway or gene drawn as part of the signature is no
#                         longer supported by the enrichment table or the signature.
# The numbers drawn in Figures 1 and 3 are passed in from the result tables.

suppressWarnings(suppressMessages({ library(ggplot2) }))

fig_dataset_prisma <- function(n_identified, n_included, file, organisms = NULL,
                               searched_on = NULL, coverage = NULL, query = NULL) {
  common <- c("Homo sapiens" = "human", "Mus musculus" = "mouse",
              "Rattus norvegicus" = "rat")
  composition <- if (is.null(organisms)) "" else {
    tb <- table(ifelse(organisms %in% names(common), common[organisms], organisms))
    paste(sprintf("%d %s", as.integer(tb), names(tb)), collapse = "   |   ")
  }
  n_from_query <- NA_integer_; n_by_hand <- NA_integer_
  hand_acc <- character(0); ho_acc <- character(0)
  if (!is.null(coverage) && nrow(coverage)) {
    cor <- coverage[coverage$role == "corpus", ]
    n_from_query <- sum(cor$returned_by_query)
    hand_acc <- cor$accession[!cor$returned_by_query]
    n_by_hand <- length(hand_acc)
    ho_acc <- coverage$accession[coverage$role == "held-out"]
  }

  INK <- "#22303F"; MUTE <- "#5A6B7B"; RULE <- "#C2CEDA"
  FILL <- "#F5F8FB"; FILL_IN <- "#E9F3EC"; EDGE_IN <- "#4E8C63"
  FILL_SIDE <- "#FCF8F0"; EDGE_SIDE <- "#AD8A52"

  qlines <- c("(kidney OR renal)[Title]   AND   (ischemi* OR reperfusion)[Title]",
              "AND   expression profiling by array / by high-throughput sequencing",
              "AND   gse[Filter]")
  if (!is.null(searched_on)) qlines <- c(qlines, "", sprintf("queried %s", searched_on))

  boxes <- list(
    list(side = FALSE, head = "IDENTIFICATION", head_col = MUTE,
         lead = sprintf("%s series returned by the declared Entrez query over GEO",
                        format(n_identified, big.mark = ",")),
         lead_size = 5.4, body = qlines, body_size = 4.4,
         fill = FILL, edge = RULE),
    list(side = FALSE, head = "INCLUSION CRITERIA", head_col = MUTE,
         lead = NULL,

         body = c("renal ischemia-reperfusion injury",
                  "injury and control groups verified sample by sample in the GEO metadata"),
         body_size = 4.8, fill = FILL, edge = RULE),
    list(side = FALSE, head = "INCLUDED", head_col = EDGE_IN,
         lead = sprintf("%d datasets analysed", n_included), lead_size = 6.2,
         body = composition, body_size = 4.8, fill = FILL_IN, edge = EDGE_IN))

  sides <- list()
  if (!is.na(n_by_hand) && n_by_hand > 0)
    sides[[1]] <- list(anchor = 1L, head = "IDENTIFIED SEPARATELY", head_col = EDGE_SIDE,
                       lead = sprintf("%d human transplant series", n_by_hand),
                       lead_size = 4.7,
                       body = c("not returned by the declared query",
                                paste(hand_acc, collapse = "   |   "), "",
                                "their titles carry neither",
                                '"ischemia" nor "reperfusion"'),
                       body_size = 4.3, fill = FILL_SIDE, edge = EDGE_SIDE)
  if (length(ho_acc))
    sides[[length(sides) + 1L]] <-
      list(anchor = 3L, head = "HELD OUT", head_col = MUTE,
           lead = paste(ho_acc, collapse = "   |   "), lead_size = 4.7,
           body = c("returned by the declared query,",
                    "reserved for external validation"),
           body_size = 4.3, fill = FILL, edge = RULE)

  HEAD_SIZE <- 4.7
  LHF <- 1.30
  line_in <- function(size, lh = LHF) size * .pt * lh / 72.27
  PAD <- 0.16
  hgt <- function(b) {
    h <- line_in(HEAD_SIZE, 1.55) +
         (if (is.null(b$lead)) 0 else line_in(b$lead_size, 1.70)) +
         (if (length(b$body) && any(nzchar(b$body)))
            length(b$body) * line_in(b$body_size) else 0)
    max(0.50, h / 2 + PAD)
  }
  bh <- vapply(boxes, hgt, numeric(1))
  sh <- if (length(sides)) vapply(sides, hgt, numeric(1)) else numeric(0)

  GAP <- 0.52
  ytop <- 0
  yc <- numeric(3)
  for (i in seq_along(boxes)) {
    yc[i] <- ytop - bh[i]
    ytop <- yc[i] - bh[i] - GAP
  }
  bw <- 3.70; sw <- 2.35; xs <- 0; xside <- bw + 0.85 + sw

  rects <- data.frame(
    xmin = xs - bw, xmax = xs + bw, ymin = yc - bh, ymax = yc + bh,
    fill = vapply(boxes, function(b) b$fill, character(1)),
    edge = vapply(boxes, function(b) b$edge, character(1)), stringsAsFactors = FALSE)

  p <- ggplot() +
    geom_segment(data = data.frame(y = (yc - bh)[-3], yend = (yc + bh)[-1]),
                 aes(x = xs, xend = xs, y = y, yend = yend),
                 arrow = arrow(length = unit(0.24, "cm"), type = "closed"),
                 color = INK, linewidth = 0.6) +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = rects$fill, color = rects$edge, linewidth = 0.7)

  if (length(sides)) {
    sy <- vapply(sides, function(b) yc[b$anchor], numeric(1))
    p <- p +
      geom_segment(data = data.frame(y = sy),
                   aes(x = xs + bw, xend = xside - sw, y = y, yend = y),
                   arrow = arrow(length = unit(0.19, "cm"), type = "closed"),
                   color = MUTE, linewidth = 0.45, linetype = "22") +
      geom_rect(data = data.frame(y = sy, h = sh,
                                  fill = vapply(sides, function(b) b$fill, character(1)),
                                  edge = vapply(sides, function(b) b$edge, character(1)),
                                  stringsAsFactors = FALSE),
                aes(xmin = xside - sw, xmax = xside + sw, ymin = y - h, ymax = y + h),
                fill = vapply(sides, function(b) b$fill, character(1)),
                color = vapply(sides, function(b) b$edge, character(1)), linewidth = 0.6)
  }

  draw <- function(p, b, x, ymid, halfh) {
    cur <- ymid + halfh - PAD - line_in(HEAD_SIZE) / 2
    p <- p + annotate("text", x = x, y = cur, label = b$head, size = HEAD_SIZE,
                      fontface = "bold", color = b$head_col, family = "sans")
    cur <- cur - line_in(HEAD_SIZE, 1.55)
    if (!is.null(b$lead)) {
      p <- p + annotate("text", x = x, y = cur, label = b$lead, size = b$lead_size,
                        fontface = "bold", color = INK, family = "sans")
      cur <- cur - line_in(b$lead_size, 1.70)
    }
    if (length(b$body) && any(nzchar(b$body))) {
      p <- p + annotate("text", x = x, y = cur + line_in(b$body_size) / 2,
                        label = paste(b$body, collapse = "\n"),
                        size = b$body_size, color = MUTE, family = "sans",
                        vjust = 1, lineheight = LHF)
    }
    p
  }
  for (i in seq_along(boxes)) p <- draw(p, boxes[[i]], xs, yc[i], bh[i])
  for (i in seq_along(sides))
    p <- draw(p, sides[[i]], xside, yc[sides[[i]]$anchor], sh[i])

  ylim <- c(min(c(yc - bh, if (length(sh)) sy - sh)) - 0.20,
            max(c(yc + bh, if (length(sh)) sy + sh)) + 0.20)

  TOP_MATTER <- 1.32

  p <- p +
    coord_cartesian(xlim = c(xs - bw - 0.25, xside + sw + 0.25), ylim = ylim) +
    theme_void(base_size = 13) +
    theme(plot.margin = margin(22, 24, 20, 24),
          plot.background = element_rect(fill = "white", color = NA))
  ggsave(file, p, width = 13.8, height = diff(ylim) + TOP_MATTER, dpi = 300)
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
    list(head = "META-ANALYSIS", head_col = MUTE, fill = FILL, edge = RULE,
         lead = "RobustRankAggreg", lead_size = LEAD_SIZE,
         body = "rodent genes mapped to their human orthologs", body_size = BODY_SIZE),
    list(head = "CONSENSUS SIGNATURE", head_col = EDGE_SIG, fill = FILL_SIG, edge = EDGE_SIG,
         lead = sprintf("%s genes", n(n_sig, "n")), lead_size = LEAD_SIZE,
         body = "measured in at least 5 of 9   |   FDR < 0.05   |   consistent direction",
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
  p <- p + key(1.55, "sig", "In the signature") + key(4.45, "lit", "From the literature") +
    key(7.55, "int", "Candidate intervention")

  ggsave(file, p, width = 12.2, height = 7.12, dpi = 300, device = ragg::agg_png, bg = "white")
  p
}
