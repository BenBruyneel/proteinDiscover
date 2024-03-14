# proteinDiscover

R Package that provides the ability to read results from .pdResult (.msf) files produced by the Proteome Discoverer software of [Thermo Scientific](https://www.thermoscientific.com/)

[proteinDiscover manual](https://benbruyneel.github.io/proteinDiscover/)

[proteinDiscover news 22/05/30](https://benbruyneel.github.io/proteinDiscover/updates220530/)

[proteinDiscover news 23/10/31](https://benbruyneel.github.io/proteinDiscover/updates231031/)

[proteinDiscoverExtra](https://github.com/BenBruyneel/proteinDiscoverExtra)

Install with the command:

devtools::install_github("BenBruyneel/proteinDiscover")

Update (v0.10.0) : revisited all functions and attempted to make argument names etc more consistent. Please note that this can have consequences for older code using these functions

Update (v0.11.0) : newer versions of Proteome Discoverer (latest version I have seen is v3.1) have a different approach to storing some data. Most of the tables are still in the .pdResult files, but (MS2 & MS3) spectrum information is stored in the .pdResultDetails files. Please note that the current implementation of all spectrum functions in this package is not final.  

Work in progress!

March 13, 2024

