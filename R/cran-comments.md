This is a new release. R CMD was run locally (windows), and using GitHub Actions
on Ubuntu, Mac-OS and Windows for R-latest and R-devel.

Version 1.0.3 prevents the issue leading to the error on 
r-devel-linux-x86_64-fedora-gcc for 1.0.2, which was caused by an error in the
data.table package that has been fixed in data.table 1.18.0.

Additionally, R CMD was run with `_R_CHECK_DEPENDS_ONLY_ = True`


# R CMD check results

0 errors | 0 warnings | 2 notes

❯ checking CRAN incoming feasibility ... [12s] NOTE
  Maintainer: 'Alex Zwanenburg <alexander.zwanenburg@nct-dresden.de>'
  
  New submission
  
  Package was archived on CRAN
  
  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2026-01-15 as issues were not addressed
      in time.

❯ checking for future file timestamps ... NOTE
  unable to verify current time
