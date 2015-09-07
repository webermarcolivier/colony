# Colony Software

## Remarks

- A lot of documentation and comments are already included in the mathematica notebook.

## List of future improvements

- Fix and test the GUI. Small modifications may be necessary to adapt to the new version of Qt (not much work).
- Implement XML input file support (much more flexible than hard-coded by line text file reading).
- More coherent way of setting options.
    - Set the compilations options in debug.h to be written by Mathematica notebook, to be more coherent and avoid possible sources of user errors (as for example setting a simulation for a growing population and compiling with the fixed population option with undefined TIME_DEPENDENT_PROPENSITIES).
    - rename the file debug.h to compilation_options.h and change \#ifdef statements accordingly.