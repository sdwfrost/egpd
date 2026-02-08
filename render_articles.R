#!/usr/bin/env Rscript
# Render all vignettes to GitHub-flavored markdown in articles/

rmd_files <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)

if (length(rmd_files) == 0L) {
  stop("No .Rmd files found in vignettes/")
}

dir.create("articles", showWarnings = FALSE)

for (rmd in rmd_files) {
  message("Rendering ", rmd, " ...")
  rmarkdown::render(
    rmd,
    output_format = rmarkdown::github_document(html_preview = FALSE),
    output_dir = "articles",
    quiet = TRUE
  )
}

md_files <- list.files("articles", pattern = "\\.md$")
message("Done. Rendered ", length(md_files), " article(s): ",
        paste(md_files, collapse = ", "))
