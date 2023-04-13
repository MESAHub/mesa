#!/bin/bash

#### MESA SHMESA 
#### Command line utilies for MESA 
# provides commands such as `mesa change` and `mesa grep` 
# for usage, source this file (`source shmesa.sh`) and call: mesa help 
# hot tip: add `source $MESA_DIR/scripts/shmesa/shmesa.sh` to your ~/.bashrc 

export MESA_SHMESA_DEBUG=0 # set to 1 for commentary 
export MESA_SHMESA_BACKUP=1 # back up modified files before modification (e.g. to inlist.bak) 

mesa () {
    ( # start a subshell so that the `set` commands are not persistent)
    set -Eeo pipefail # exit if any commands fail 
    #if [[ $MESA_SHMESA_DEBUG -eq 1 ]]; then # print each command before it is executed
    #    set -x
    #else
    #    set +x
    #fi

    if [[ -z $MESA_DIR ]] || [[ ! -d $MESA_DIR ]]; then
        echo "Error: \$MESA_DIR is not set or does not point to a valid directory."
        echo "       \$MESA_DIR=$MESA_DIR"
        echo "       Please download and install MESA:"
        echo "https://docs.mesastar.org"
        exit 1
    fi
    
    ### Define main utilities; parse command line tokens at the end 
    mesa_help () {
         cat << "EOF"
      _     __  __ _____ ____    _    
  ___| |__ |  \/  | ____/ ___|  / \   
 / __| '_ \| |\/| |  _| \___ \ / _ \  
 \__ \ | | | |  | | |___ ___) / ___ \ 
 |___/_| |_|_|  |_|_____|____/_/   \_\
                                      
EOF
        echo "Usage: mesa [work|change|defaults|cp|grep|zip|version|update|help] [arguments]"
        echo
        echo "Subcommands:"
        echo "  work      copy the work directory to the current location"
        echo "  change    change a parameter in the given inlist"
        echo "  defaults  copy the history/profile defaults to the current location"
        echo "  cp        copy a MESA directory without copying LOGS, photos, etc."
        echo "  grep      search the MESA source code for a given string"
        echo "  zip       prepare a MESA directory for sharing"
        #echo "  version   check the version of MESA (and check for updates)//TODO"
        #echo "  update    update MESA to the latest version //TODO"
        echo "  help      display this helpful message"
        echo "  -h        flag for getting additional details about any of the above"
        echo
    }
    
    ### Now define all the subcommands in the above help message 

    mesa_work () {
        if mesa_check_h_flag "$@"; then
            echo "Usage: mesa work [optional: target_name]"
            echo "Copies the MESA star/work directory to the specified location."
            echo "If no target_name is provided, it will copy the directory to the current location."
            return 0
        fi

        local target_dir="."
        if [[ -n $1 ]]; then # check for nonempty arg
            target_dir=$1
        fi

        if [[ -d $target_dir ]]; then
            echo "Error: Target directory '$target_dir' already exists."
            return 1
        fi

        cp -R "$MESA_DIR/star/work" "$target_dir"
    }


    mesa_change () {
        if mesa_check_h_flag "$@"; then
            echo "Usage: mesa change inlist parameter value [parameter value [parameter value]]"
            echo "Modifies one or more parameters in the supplied inlist."
            echo "Uncomments the parameter if it's commented out."
            echo "Creates a backup of the inlist in the corresponding .bak file."
            echo ""
            echo "Examples:"
            echo "  mesa change inlist_project initial_mass 1.3"
            echo "  mesa change inlist_project log_directory 'LOGS_MS'"
            echo "  mesa change inlist_project do_element_diffusion .true."
            echo "  mesa change inlist_project initial_mass 1.3 do_element_diffusion .true."
            return 0
        fi

        if [[ -z $1 || -z $2 || -z $3 ]]; then
            echo "Error: Missing arguments."
            return 1
        fi

        local filename=$1
        shift

        if [[ ! -f "$filename" ]]; then
            echo "Error: '$filename' does not exist or is not a file."
            exit 1
        fi

        # Create a backup of the inlist before making any changes
        backup_copy "$filename" "${filename}.bak"

        local ESCAPE="s#[^^]#[&]#g; s#\^#\\^#g" # sed escape string
        # iterate through parameter,value pairs in the supplied input 
        while [[ -n $1 && -n $2 ]]; do
            local param=$1
            local newval=$2
            shift 2

            local escapedParam=$(sed "$ESCAPE" <<< "$param")
            local search="^\s*\!*\s*$escapedParam\s*=.+$" 
            local replace="    $param = $newval"

            # Check if the parameter is present in the inlist
            if ! grep -q "$param" "$filename"; then
                echo "Error: Parameter '$param' not found in the inlist '$filename'."
                backup_restore "$filename"
                return 1
            fi

            # Update its value
            sed -r -i -e "s#$search#$replace#g" "$filename"
        done
    }


    mesa_defaults () {
        if mesa_check_h_flag "$@"; then
            echo "Usage: mesa defaults [parameter [parameter]]"
            echo "Copies profile_columns.list and history_columns.list to the current location."
            echo "Also uncomments any specified parameters."
            echo "If the files are already in the present directory, just uncomment the specified parameters."
            echo ""
            echo "Example: mesa defaults nu_max delta_nu"
            return 0
        fi

        # Copy the files; create backups if they already exist in the current directory
        if [[ -f profile_columns.list ]]; then
            backup_copy profile_columns.list
        else
            cp "$MESA_DIR/star/defaults/profile_columns.list" .
        fi

        if [[ -f history_columns.list ]]; then
            backup_copy history_columns.list
        else
            cp "$MESA_DIR/star/defaults/history_columns.list" .
        fi

        local ESCAPE="s#[^^]#[&]#g; s#\^#\\^#g" # sed escape string
        # Uncomment the specified parameters
        while [[ -n $1 ]]; do
            local param=$1
            shift

            local escapedParam=$(sed "$ESCAPE" <<< "$param")
            local search="^\s*\!*\s*$escapedParam\s*.+$"
            local replace="    $param"

            if grep -q "$param" profile_columns.list || \
               grep -q "$param" history_columns.list; then
                sed -r -i -e "s#$search#$replace#g" profile_columns.list
                sed -r -i -e "s#$search#$replace#g" history_columns.list
            else
                echo "Warning: Parameter '$param' not found in either profile_columns.list or history_columns.list."
            fi
        done
    }


    mesa_cp () {
        if mesa_check_h_flag "$@"; then
            echo "Usage: mesa cp source_dir target_dir"
            echo "Copies a MESA working directory but without copying"
            echo "LOGS, photos, or .mesa_temp_cache"
            return 0
        fi

        if [[ -z $1 || -z $2 ]]; then
            echo "Error: Missing arguments."
            echo "Usage: mesa cp source_dir target_dir"
            return 1
        fi

        # args: ($1) source directory to be copied from
        #       ($2) target directory to be copied to
        local SOURCE=$1
        local TARGET=$2
        rsync -av --progress $SOURCE/ $TARGET \
            --exclude LOGS \
            --exclude photos \
            --exclude .mesa_temp_cache
    }
    

    mesa_grep () {
        if mesa_check_h_flag "$@"; then
            echo "Usage: mesa grep term [optional: directory or filename]"
            return 0
        fi

        if [[ -z $1 ]]; then
            echo "Error: Missing search term."
            echo "Usage: mesa grep term [optional: directory or filename]"
            return 1
        fi

        local search_term=$1
        local search_dir=$MESA_DIR

        if [[ -n $2 ]]; then
            search_dir=$MESA_DIR/$2
        fi

        grep -r --color=always "$search_term" "$search_dir"
    }


    mesa_zip () {
        if mesa_check_h_flag "$@"; then
            echo "Usage: mesa zip [directory]"
            echo "zips the inlists and models of the specified directory for sharing"
            return 0
        fi

        UNTESTED

        local dir_to_zip='.'
        local zip_name="mesa.zip"
        if [[ ! -z $1 ]]; then
            dir_to_zip=$1
            zip_name="${dir_to_zip}_mesa.zip"
        fi

        if [[ ! $(find "$dir_to_zip" -type f -iname "*inlist*" -o -iname "*.mod") ]]; then
            echo "Warning: No files with 'inlist' in their names or files ending in '.mod' found in the directory '$dir_to_zip'."
        fi

        # Create a zip file containing only inlists and models
        zip -r "$zip_name" "$dir_to_zip" -i '*inlist*' '*.mod'
    }


    ### DEVELOPMENT USE 
    mesa_template () {
        ### TEMPLATE FOR NEW SHMESA FUNCTIONS 
        if mesa_check_h_flag "$@"; then
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
        # 1. mesa_test 
        # 2. mesa_help 
        # 3. parser at the end of this file 
    }


    ## Test shmesa with different subcommands and arguments
    mesa_test () {
        echo
        echo 
        echo 
        echo
        echo ">> TESTING SHMESA <<"
        echo
        echo

        # store current value of MESA_SHMESA_DEBUG and turn on debugging 
        local temp_value=$MESA_SHMESA_DEBUG
        if [[ $# -gt 0 ]]; then # decide whether we want debugging 
            MESA_SHMESA_DEBUG=$1
        else
            MESA_SHMESA_DEBUG=1
        fi 
        
        # mesa [work|change|defaults|cp|grep|zip|version|update|help]
        local SHMESA_TEST_DIR="shmesa_test_work"

        # mesa work
        run_mesa_test "work" \
            "mesa work $SHMESA_TEST_DIR && \
            cd $SHMESA_TEST_DIR && \
            ./mk"

        # mesa change 
        run_mesa_test "change" \
            "mesa change inlist_project \
                         pgstar_flag .false. && \
             mesa change inlist_project \
                         initial_mass 1.2 \
                         Z 0.01 \
                         Zbase 0.01"
        
        # mesa defaults
        run_mesa_test "defaults" \
            "mesa defaults nu_max delta_nu && \
             mesa defaults logM"

        # mesa cp
        run_mesa_test "cp" \
            "cd .. &&\
             mesa cp '$SHMESA_TEST_DIR' '$SHMESA_TEST_DIR'2"

        # mesa grep
        # TODO

        # mesa zip
        # TODO

        # mesa version
        # TODO

        # mesa update
        # TODO

        MESA_SHMESA_DEBUG=$temp_value
        echo "all done!"
    }


    ## Some helper functions used above for debugging and backing up files 
    
    debug_print () {
        if [[ -n $MESA_SHMESA_DEBUG && $MESA_SHMESA_DEBUG -ne 0 ]]; then
            echo "SHMESA> DEBUG: $@"
        fi
    }
    
    backup_copy () {
        if [[ $MESA_SHMESA_BACKUP -ne 0 && ! -z $1 ]]; then
            debug_print "SHMESA> BACKING UP: $@"
            cp "$1" "$1".bak
        else
            debug_print "SHMESA> NOT BACKING UP: $@ (\$MESA_SHMESA_BACKUP=$MESA_SHMESA_BACKUP)"
        fi
    }
    
    backup_restore () {
        if [[ $MESA_SHMESA_BACKUP -ne 0 && ! -z $1 ]]; then
            debug_print "SHMESA> RESTORING: $@"
            cp "$1".bak "$1"
        else
            debug_print "SHMESA> NOT RESTORING: $@ (\$MESA_SHMESA_BACKUP=$MESA_SHMESA_BACKUP)"
        fi
    }

    mesa_check_h_flag () {
        for arg in "$@"; do
            if [[ $arg == "-h" ]]; then
                return 0
            fi
        done
        return 1
    }

    run_mesa_test() {
        local submodule_name=$1
        local test_command=$2

        echo
        echo
        echo ">>> test: mesa $submodule_name"
        echo
        eval "$test_command"
        echo
        echo "<<< success"
        echo
        echo
    }

    UNTESTED () {
        echo "WARNING: UNTESTED"
        exit 1
    }


    #############################
    ### Parse command line tokens
    ###
    if [[ -z $1 ]]; then
        mesa_help
        return 1
    fi

    local subcommand=$1
    shift

    debug_print "Calling $subcommand with arguments: $@"

    case "$subcommand" in
        work)
            mesa_work "$@"
            ;;
        change)
            mesa_change "$@"
            ;;
        defaults)
            mesa_defaults "$@"
            ;;
        cp)
            mesa_cp "$@"
            ;;
        grep)
            mesa_grep "$@"
            ;;
        zip)
            mesa_zip "$@"
            ;;
        test)
            mesa_test "$@"
            ;;
        help)
            mesa_help
            ;;
        "" | *[[:space:]]*)
            mesa_help
            ;;
        *)
            
            echo "Invalid subcommand: $subcommand"
            # TODO: "The most similar command is: ..."
            mesa_help
            return 1
        ;;
    esac
    ) # close subshell
}
