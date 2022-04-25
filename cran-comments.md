## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Florian Erhard <florian.erhard@uni-wuerzburg.de>’
                                                                                                                              
New submission

* \dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user. Does not seem necessary. Please unwrap the examples if they are executable in < 5 sec, or replace \dontrun{} with \donttest{}. 

The example for PsiLFC.se takes > 5 sec to execute and is therefore wrapped into \dontrun

