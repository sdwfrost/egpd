#!/usr/bin/env Rscript
# Render all vignettes to GitHub-flavored markdown in articles/

rmd_files <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)

if (length(rmd_files) == 0L) {
  stop("No .Rmd files found in vignettes/")
}

dir.create("articles", showWarnings = FALSE)
articles_dir <- normalizePath("articles")

for (rmd in rmd_files) {
  message("Rendering ", rmd, " ...")

  # Pre-process: override collapse/comment for GFM-friendly output
  txt <- readLines(rmd)
  txt <- sub("collapse\\s*=\\s*TRUE", "collapse = FALSE", txt)
  txt <- sub('comment\\s*=\\s*"#>"', 'comment = ""', txt)
  tmp <- tempfile(fileext = ".Rmd")
  writeLines(txt, tmp)

  fmt <- rmarkdown::github_document(html_preview = FALSE)
  rmarkdown::render(
    tmp,
    output_format = fmt,
    output_dir = articles_dir,
    quiet = TRUE
  )

  # Rename output from temp filename to original vignette name
  tmp_md <- file.path(articles_dir,
                      sub("\\.Rmd$", ".md", basename(tmp)))
  out_md <- file.path(articles_dir,
                      sub("\\.Rmd$", ".md", basename(rmd)))
  # Also rename figure directory
  tmp_figs <- file.path(articles_dir,
                        sub("\\.Rmd$", "_files", basename(tmp)))
  out_figs <- file.path(articles_dir,
                        sub("\\.Rmd$", "_files", basename(rmd)))
  if (file.exists(out_md)) file.remove(out_md)
  file.rename(tmp_md, out_md)
  if (dir.exists(tmp_figs)) {
    if (dir.exists(out_figs)) unlink(out_figs, recursive = TRUE)
    file.rename(tmp_figs, out_figs)
    # Fix figure directory references inside the md
    md_txt <- readLines(out_md)
    md_txt <- gsub(basename(tmp_figs), basename(out_figs), md_txt,
                   fixed = TRUE)
    writeLines(md_txt, out_md)
  }
  unlink(tmp)
}

# Post-process: fix absolute figure paths to relative and remove <!-- --> tags
md_files <- list.files(articles_dir, pattern = "\\.md$", full.names = TRUE)
abs_prefix <- paste0(articles_dir, "/")
for (md in md_files) {
  txt <- readLines(md)
  txt <- gsub(abs_prefix, "", txt, fixed = TRUE)
  txt <- gsub("<!-- -->", "", txt, fixed = TRUE)
  writeLines(txt, md)
}

message("Done. Rendered ", length(md_files), " article(s): ",
        paste(basename(md_files), collapse = ", "))
