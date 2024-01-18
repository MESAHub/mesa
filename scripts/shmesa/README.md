**shmesa**: a set of bash utilities providing convenience functions for MESA 

To install, simply add `PATH=$PATH:$MESA_DIR/scripts/shmesa` to your `~/.bashrc` and start a new terminal.

To use, simply call: `shmesa`

```bash
> shmesa
      _     __  __ _____ ____    _
  ___| |__ |  \/  | ____/ ___|  / \
 / __| '_ \| |\/| |  _| \___ \ / _ \
 \__ \ | | | |  | | |___ ___) / ___ \
 |___/_| |_|_|  |_|_____|____/_/   \_\

Usage: shmesa [work|change|defaults|cp|grep|zip|help] [arguments]

Subcommands:
  work      copy the work directory to the current location
  change    change a parameter in the given inlist
  defaults  copy the history/profile defaults to the current location
  cp        copy a MESA directory without copying LOGS, photos, etc.
  grep      search the MESA source code for a given string
  zip       prepare a MESA directory for sharing
  help      display this helpful message
  -h        get additional details about any of the above
```

Want more details about any of the subprograms? Just run any of them with `-h`, for example 
```bash
> shmesa cp -h
Usage: shmesa cp source_dir target_dir
Copies a MESA working directory but without copying LOGS, photos, or .mesa_temp_cache
```

### Some example usage 
You can use this for example to make a grid of models, and e.g. include the output `delta_nu` and `nu_max` in `LOGS` with 
```bash
### Make a grid in M and Z in a few lines of bash

# copy over a new work directory 
shmesa work "grid_dir"
cd grid_dir
./mk

# output asteroseismic quantities in `LOGS/history.data` 
shmesa defaults delta_nu nu_max

# loop through values of M and Z and run MESA 
for M in 1.0 1.5 2.0; do 
    for Z in 0.01 0.02 0.03; do 
        shmesa change inlist_project \
            initial_mass $M \
            Zbase $Z \
            initial_z $Z
        ./star inlist_project 
        mv LOGS "M='$M'_'Z='$Z'" 
    done 
done 

```


### Tests
There is functionality for testing all the subprograms:

```bash
> shmesa test
>> TESTING SHMESA <<

>>> test: shmesa work
SHMESA> DEBUG: Calling work with arguments: shmesa_test_work
gfortran -Wno-uninitialized -fno-range-check -fmax-errors=7  -fprotect-parens -fno-sign-zero ...
<<< success


>>> test: shmesa change
SHMESA> DEBUG: Calling change with arguments: inlist_project pgstar_flag .false.
SHMESA> DEBUG: SHMESA> BACKING UP: inlist_project inlist_project.bak
SHMESA> DEBUG: Calling change with arguments: inlist_project initial_mass 1.2 Z 0.01 Zbase 0.01
SHMESA> DEBUG: SHMESA> BACKING UP: inlist_project inlist_project.bak
<<< success

...

>>> test: shmesa defaults
SHMESA> DEBUG: Calling defaults with arguments: nu_max delta_nu
SHMESA> DEBUG: Calling defaults with arguments: logM
SHMESA> DEBUG: SHMESA> BACKING UP: profile_columns.list
SHMESA> DEBUG: SHMESA> BACKING UP: history_columns.list
<<< success

all done!
```

### Development
There's instructions inside the code for adding a new bash subprogram to this driver: 
```bash
    ### DEVELOPMENT USE 
    ### TEMPLATE FOR NEW SHMESA FUNCTIONS 
    shmesa_template () {
        if shmesa_check_h_flag "$@"; then
            echo "Usage: mesa funcname arg [optional arg]"
            echo "put discription here"
            return 0
        fi

        UNTESTED # comment this out when done (crashes the program to avoid problems)

        # parse required variables 
        local example_var=5 
        if [[ ! -z $1 ]]; then # check if it's set 
            example_var=$1
        #else # uncomment if you want this to be a required variable
        #    echo "Error: "
        #    return 1
        fi
        #shift # and then copy this block to parse the next one, for example 

        # afterwards, update:
        # 1. shmesa_test 
        # 2. shmesa_help 
        # 3. parser at the end of this file 
    }
```

