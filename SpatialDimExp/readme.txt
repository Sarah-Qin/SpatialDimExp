
# Step1: File -> Open Project in New Session.
# Step2: open the readme.txt, run the codes inside.
library(devtools)
library(roxygen2)

devtools::document()  # to write the .Rd documents.

devtools::check()  # it takes most of the time.

# devtools::document(roclets = c('rd', 'collate', 'namespace'))
devtools::build(binary = TRUE, args = c('--preclean'))

devtools::build_manual(pkg = ".", path = NULL)

