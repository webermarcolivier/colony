# Colony Software

## Remarks

- A lot of documentation and comments are already included in the mathematica notebook.

### List of future improvements

- Fix and test the GUI. Small modifications may be necessary to adapt to the new version of Qt (not much work).
- Implement XML input file support (much more flexible than hard-coded by line text file reading).
- Optimize memory usage. Currently, one trajectory is kept into memory in the form of a list (one element per time mesh). For very long trajectories (simulation time), this results into increasing memory usage. Possible optimizations:
    + check if some intermediate time steps are kept in memory. This is probably the case because increasing the time mesh interval does not seem to decrease memory usage for the same simulation time.
    + One possible optimization for the code would be to write trajectory output to the file every several time mesh (chunk) and keep in memory only the current trajectory chunk. However, this is a little cumbersome to implement.