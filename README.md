# proteinDiscover

R Package that provides the ability to read results from .pdResult (.msf) files produced by the Proteome Discoverer software of [Thermo Scientific](https://www.thermoscientific.com/)

[proteinDiscover manual](https://benbruyneel.github.io/proteinDiscover/)

[proteinDiscover news 22/05/30](https://benbruyneel.github.io/proteinDiscover/updates220530/)

[proteinDiscover news 23/10/31](https://benbruyneel.github.io/proteinDiscover/updates231031/)

[proteinDiscoverExtra](https://github.com/BenBruyneel/proteinDiscoverExtra)

Install with the command:

devtools::install_github("BenBruyneel/proteinDiscover")

Latest update (v0.10.0) : revisited all functions and attempted to make argument names etc more consistent. Please note that this can have consequences for older code using these functions

March 11, 2024: Please note that pdResult files coming from Proteome Discoverer 3.1 do not seem to contain the table 'MassSpectrumItems' anymore. Until I figure out where the MS2 spectra information is located, functions working with this table (and the previously related tables) will NOT work. 

Work in progress!

March 11, 2024

