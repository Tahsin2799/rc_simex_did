# Build paper 1's review draft and its journal files: the all-in-one document (main text, appendices and
# supplement in one PDF), and the main text and the supplement as two PDFs that read each other's labels (xr).
#   Rscript paper1_rc_simex_did/submission/R/build_submission.R    (from the repository root, after build_exhibits.R;
#   manuscript/R/... in the code release)
# main.Rmd and supplement.Rmd are generated here from the header of paper_all_in_one.Rmd, so the three documents
# share one preamble; edit the header there, the text in sections/main.Rmd and the appendices in appendix/.
suppressPackageStartupMessages({
  library(rmarkdown)
})
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
SUB <- if (length(script)) dirname(dirname(script)) else "paper1_rc_simex_did/submission"

# ---- Generate the two journal files ---------------------------------------------
aio <- readLines(file.path(SUB, "paper_all_in_one.Rmd"), encoding = "UTF-8")
yaml_end <- which(aio == "---")[2]
if (is.na(yaml_end)) cli::cli_abort("No closing YAML fence in {.file paper_all_in_one.Rmd}.")
yaml <- aio[seq_len(yaml_end - 1L)]
setup_at <- grep("^```\\{r setup", aio)
setup <- aio[setup_at:(setup_at + which(aio[-seq_len(setup_at)] == "```")[1])]
children <- grep("^```\\{r child = \"appendix/", aio, value = TRUE)
with_xr <- \(y, other, extra = character()) c(y, "  - \\usepackage{xr}", paste0("  - \\externaldocument{", other, "}"), extra, "---")
drop_abstract <- function(y) {
  a <- grep("^abstract:", y)
  if (!length(a)) return(y)
  nxt <- a + which(!grepl("^\\s", y[-seq_len(a)]))[1]
  y[-(a:(nxt - 1L))]
}
title_at <- grep("^title:", yaml)
supp_yaml <- drop_abstract(yaml)
supp_yaml <- append(supp_yaml, sub("^title:", "subtitle:", yaml[title_at]), after = grep("^title:", supp_yaml))
supp_yaml[grep("^title:", supp_yaml)] <- 'title: "Supplementary material"'
writeLines(c(with_xr(yaml, "supplement"), "", setup, "", "```{r child = \"sections/main.Rmd\"}", "```", "",
             "\\bibliography{references}"), file.path(SUB, "main.Rmd"), useBytes = TRUE)
writeLines(c(with_xr(supp_yaml, "main", c("  - \\renewcommand{\\thetable}{S\\arabic{table}}",
                                           "  - \\renewcommand{\\thefigure}{S\\arabic{figure}}")),
             "", setup, "", "\\appendix", "", unlist(lapply(children, \(ch) c(ch, "```", ""))),
             "\\bibliography{references}"), file.path(SUB, "supplement.Rmd"), useBytes = TRUE)

# ---- Render ---------------------------------------------------------------------
# Knit each journal file once, then compile the pair twice with the auxiliary files kept, so the cross-document
# references resolve; the all-in-one document is a single self-contained render.
for (f in c("supplement.Rmd", "main.Rmd")) render(file.path(SUB, f), quiet = TRUE)
compile <- \(tex) withr::with_dir(SUB, tinytex::latexmk(tex, engine = "pdflatex", bib_engine = "bibtex", clean = FALSE))
for (pass in seq_len(2)) for (tex in c("supplement.tex", "main.tex")) compile(tex)
render(file.path(SUB, "paper_all_in_one.Rmd"), quiet = TRUE, clean = FALSE)

# ---- Scorecard ------------------------------------------------------------------
pages <- \(pdf) as.integer(sub(".*:\\s*", "", grep("^Pages:", system2("pdfinfo", pdf, stdout = TRUE), value = TRUE)))
undef <- \(log) sum(grepl("undefined", readLines(log, warn = FALSE), fixed = TRUE, useBytes = TRUE))
docs <- c("main", "supplement", "paper_all_in_one")
writeLines(sprintf("%-17s %3d pages, %d lines mentioning undefined references", docs,
                   vapply(file.path(SUB, paste0(docs, ".pdf")), pages, integer(1)),
                   vapply(file.path(SUB, paste0(docs, ".log")), undef, integer(1))))
