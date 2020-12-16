      block data blstio
c
c  define unit numbers for standard input and output
c
c  istdin: standard input (typically terminal)
c  istdou: standard output for diagnostics, etc (typically terminal)
c  istdpr: printed output (terminal or printer)
c
c  For testing purposes, the initial default values are stored 
c  separately in common/cstdio_def/
c
c
c  hp9000 version
c  ************
c
      common/cstdio/ istdin, istdou, istdpr, istder
      common/cstdio_def/ istdin_def, istdou_def, istdpr_def, istder_def
c
      data istdin, istdou, istdpr, istder
     *  /    5,       6,     6,      0    /
      data istdin_def, istdou_def, istdpr_def, istder_def
     *  /       5,          6,         6,           0    /
      end
