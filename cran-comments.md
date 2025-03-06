## Resubmission 2

This is a resubmission. In response to CRAN requests I have:

* Removed single quotes from title and description as these did not refer to package names.
* Added authors(year) before the doi in the description.
* Replaced all T and F with TRUE and FALSE throughout the package.
* Removed the \dontrun{} examples from the documentation.
* Placed small data objects in inst/extdata/ and removed the data/ folder altogether as it was only used for examples.
* Used system.file() to access the files in extdata for the examples, tests, vignette and readme.  

With these changes, the R CMD checks now suggest the following possibly misspelled words in DESCRIPTION: Omics, omics, et, al.
I believe these are all commonly used terms and appropriate for use in this context.
 
## Resubmission
This is a resubmission. In this version I have:

* Removed use of function shorthand for compatibility with earlier R versions.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
