# tseriesEntropy changelog


## Changes in version 0.7-2

- The Fortran sources have been cleaned and improved to comply with the new CRAN checks. The most notable change 
is the adoption of portable kind parameters, as recommended by the Fortran best practices.


## Changes in version 0.7-1

- Added a Markdown README file. 

- Removed old compiler directives in Fortran sources: they were relevant only for Windows compilers (CVF and Intel). 

- In `Srho.test.ts` and `Srho.test.ts.p`, when  `ci.type='mbb'` it is now possible to choose the moving bootstrap block length through the switch `blag`.


## Changes in version 0.7-0

- The workhorse function `Srho.func` and its Fortran engine have been rewritten as to allow generic non-diagonal bandwidth selectors. Also, the Fortran code has been vectorized and now the `hcubature` routine can take advantage from this. The routines for computing the bivariate bandwidths have been amended as to output a matrix. The old routines are still in the source code but not exported and they will be removed in a future release. A new option `bdiag` has been added to `Srho.ts` and related files. If `bdiag=FALSE`, then the bandwitdh matrix is unstructured.

- Fixed two bugs in the Fortran routine `Srhosum`.

- The help files of the parallel routines have been merged with that of the non-parallel ones and improved.

- Minor fixes to the show method of `Srho-class`.

- Bug fix: in `Srho.test` the slot `@call` now contains the correct call. 

- In `surrogate.AR` and `surrogate.ARs` `order.max` is set to `NULL`, which is the default in `stats::ar`.

- Added the horizontal line at 0 in the plot of `Srho-class` consistently with that of `Srho.test-class`.  

- In `surrogate.SA` changed the default value `eps.SA` to 0.05.

- Changed the links to the papers as to use doi instead of url.  

## Changes in version 0.6-1
    
- The integration routine name `adaptIntegrate` has been changed to its alias `hcubature` to comply with the change in the `cubature` package. 
The documentation has been amended accordingly.

- The ellipsis has been added to all the `Srho`-related routined as to allow passing additional arguments, typically to the `hcubature` routine. Incidentally, the option maxpts that was passed to maxEval has been removed.   

- Fixed a buglet in `Srho.ts.uni`.

- The Fortran code has been made cleaner and more portable by setting the machine-independent precision REAL64 through the module ISO_FORTRAN_ENV.
    
- The permutation/bootstrap versions of the test for discrete data and the `surrogate.SA` routine now call the R API C routines for random 
    number generation. The results are now fully reproducible by setting `set.seed`.
    
- Fixed `lag.max` in the documentation.

- Fixed a buglet in plot methods for S4 classes `Srho` and `Srho.test`: the switch 'main' was ignored.

- Adjusted margins in plot method for S4 class `Srho.test`.

- `nslaves` has been replaced with `nwork` throughout the package.

## Changes in version 0.6-0
    
- Added bandwidth selection methods to `Srho.ts` from package `ks`: least square cross validation (lscv), smoothed cross validation (scv), 
plugin (pi). The documentation of `Srho.ts` has been updated accordingly.

- In `Srho.ts.uni` and `Srho.ts.biv` the bandwidths of the marginal densities are computed only once (this assumes stationarity).

- Added `Srho.test.ts.p`, the parallel version of `Srho.test.ts`.

- Added the url of the Biometrika paper to the docs that had not been updated.

- Added the reference to the Biometrika paper in CITATION.

- Minor fixes to the documentation.

- All the native routines have been registered.
    

## Changes in version 0.5-12
    
- First CRAN version. Made the package compliant with the latest CRAN policy.

- Added the url of the Biometrika paper to the docs.

## Changes in version 0.5-11
    
- Moved the package parallel from Suggests to Imports.

- Changelog converted to markdown

- Minor cosmetic changes to documentation files.

## Changes in version 0.5-10
    
- Added the switch "nor" to `Srho`, `Srho.test` and related internal routines. This allows the choice over the normalization for Srho in the 
case of integer/categorical time series. Added a relevant example to the help page of `Srho`. The information has been inserted in the show 
method of `Srho` as a note.

- The Fortran subroutines SS* for computing Srho for integer/categorical time series have been modified as to use a unique workhorse 
function SRHOBIVA.

- Fortran subroutines SS,SS2,SSB have been renamed SSUNI,SSUNI2,SSUNIB for consistency.

- Simplified the routines `Srho.test.uni` and `Srho.test.biv` as to use directly the switch "stationary" to call the associated Fortran 
routines.

- Minor improvements to the docs of `Srho` and `Srho.test`.

## Changes in version 0.5-9
    
- Reverted to the native Fortran dnorm for performance reasons

- The Fortran subroutines have been reorganized in a module, many interfaces have been removed.

- Minor modifications to Srho.test.ts.Rd

## Changes in version 0.5-8
    
- The Fortran subroutines have been cleaned for portability and to remove some obsolescent constructs.

- The Fortran routines now use the C function dnorm from the R API through a C wrapper.

- Updated the references in the documentation files.

## Changes in version 0.5-7
    
- Added several parameters that allow finer control of the plots of the classes Srho and Srho.test.

## Changes in version 0.5-6
    
 - Added a robustification constant to the acf of the original series in the f90 subroutine SURROGATEACF (line 928 of the fortran source 
  code). This should solve the problem of high false positive rates.

## Changes in version 0.5-5
    
- Modified the f90 function MEVA (variance) for floating point precision accuracy.

- Parameter eps becomes tol for consistency with `adaptIntegrate` of package cubature. Docs amended accordingly.

- Paramter maxpts is passed to maxEval in `adaptIntegrate`, defaults to 0 (unlimited). Docs amended accordingly.

- Paramter minpts is suppressed. Docs amended accordingly.

## Changes in version 0.5-4
    
- Rewritten the parallel routines (`Srho.test.AR.p`, `Trho.test.AR.p`, `Trho.test.SA.p`) as to use the package parallel. The package no longer 
depends upon Rmpi and should work out of the box across platforms (including Windoze).

- Set B = 100 in most tests.

- Added a CITATION file.

- Truncated the lines of several Rd files as to avoid line truncation and notes from --as-cran.

- Fixed a bug where the parallel tests did not return the significant lags.

## Changes in version 0.5-3
    
- Moved the initial checks of `Srho.uni.*`, `Srho.biv.*` at the level of `Srho`, those of `Srho.test.uni`, `Srho.test.biv.*` at the level of 
`Srho.test`, those of `Srho.ts.uni`, `Srho.ts.biv.*` at the level of `Srho.ts`.

- `Srho.test.iid` has been renamed `Srho.test.ts` for consistency.

- Fixed p-values names in `Srho.test.ts`.

- Added further examples in Srho.test.ts.Rd.

- Built under R 3.0.2

## Changes in version 0.5-2
   
- Built under R 3.0.1

- Modified files NAMESPACE and DESCRIPTION and the parallel routines to make the package compliant to the latest import/depends directives 
as specified in the manual.

- Cleaned, updated and improved almost all the help files.

## Changes in version 0.5-1
    
- Added the switch na.action = na.exclude in `surrogate.AR` and `surrogate.ARs`.

## Changes in version 0.5-0
    
- Added the computation of the bootstrap p-value for all the tests.

- Corrected a buglet on names that caused a warning on all the tests.

## Changes in version 0.4-0
    
- Added the computation of the bootstrap p-value for `Srho.test`. Added a p.value slot to the class `Srho.test`.

- Commented out line 132 of the f90 code (reset of the random seed) as it caused problems.

- Modified the code of the examples in `Srho` and `Srho.test`.

- Added documentation for `Srho.test.AR` and `Srho.test.AR.p`.

## Changes in version 0.3-0
    
- Restored the routines `Srho.test.AR` and `Srho.test.AR.p` that have been modified as to use also the smoothed sieve which is the default.

## Changes in version 0.2-9
    
- Modified the routine `Srho` for integer/categorical data. Now `Srho` is normalized on the interval [0, 1], the maximum possible value is 
computed analytically and taken as a benchmark. This makes the results of `Srho` comparable among different series.

- Minor improvements of the Fortran routines SRHOBIVA and SS2.

- Fixed some inconsistencies on setMethod("plot") and improved setMethod("show").

## Changes in version 0.2-8
    
- Built under R 2.14.0 with rtools 2.14.

- Fixed some ylab inconsistencies in `Trho.test*`.

- Set na.action=na.pass in acf based estimation of Srho (Srho.cor).

## Changes in version 0.2-7
    
- Added `Surrogate.ARs`, the smoothed sieve bootstrap.

- Added references in Surrogate.ARs.Rd and Surrogate.AR.Rd.

- Added the option smoothed to `Trho.test.AR*`. Modified the documentation accordingly.

## Changes in version 0.2-6
    
- removed `Srho.test.SA`, `Srho.test.AR`, `Srho.test.SA.p` `Srho.test.AR.p`.

- added `Trho.test.SA.p`, `Trho.test.AR.p` with documentation.

- added documentation for `Trho.test.SA`, `Trho.test.AR`.

- added safe.Trho in the computation of `Trho.test.SA`, `Trho.test.AR`.

## Changes in version 0.2-5
    
- fixed bug on setMethod("plot"): it was not possible to set lwd.

- substituted the integration routine `adapt` with `adaptIntegrate` from package "cubature".

- changed the names `Srho.test.**2` in `Trho.test.**` to comply with the notation of the paper.

## Changes in version 0.2-4
    
- modified Srho.test.SA2 and Srho.test.AR2 to use the quadratic divergence.

## Changes in version 0.2-3
    
- added Srho.test.SA2 (Experimental).

- fixed NAMESPACE for "previous import show" warning.

- inserted the missing test type in Srho.test.biv.

- fixed several non utf8 characters in Rd files.

## Changes in version 0.2-2
    
- added Srho.test.AR2 (Experimental).

- set writeout = 200 in safe.Srho.

## Changes in version 0.2-0
    
- fixed and reworked `Srho.test.iid` and its documentation; included the `mbb` option into the routine.

- Rmpi now in the Suggests list (was in the Depends one).

## Changes in version 0.1-5
    
- added `Srho.test.iid`.

- improved and cleaned the documentation.

- modified safe.Srho to deal with the bivariate case.
