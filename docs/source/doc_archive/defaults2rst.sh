#!/bin/bash

defaults2rst(){
    sed 's/^ *\! \?//' "$1" > "$2"
}

# astero
defaults2rst ../astero/defaults/astero_search.defaults source/reference/astero_search.rst
defaults2rst ../astero/defaults/astero_pgstar.defaults source/reference/astero_pgstar.rst

# eos
defaults2rst ../eos/defaults/eos.defaults source/eos/defaults.rst

# kap
defaults2rst ../kap/defaults/kap.defaults source/kap/defaults.rst

# binary
defaults2rst ../binary/defaults/binary_job.defaults source/reference/binary_job.rst
defaults2rst ../binary/defaults/binary_controls.defaults source/reference/binary_controls.rst
defaults2rst ../binary/defaults/pgbinary.defaults source/reference/pgbinary.rst


# star
defaults2rst ../star/defaults/star_job.defaults source/reference/star_job.rst
defaults2rst ../star/defaults/controls.defaults source/reference/controls.rst
defaults2rst ../star/defaults/pgstar.defaults source/reference/pgstar.rst
defaults2rst ../star/defaults/FORMAT source/reference/format.rst
