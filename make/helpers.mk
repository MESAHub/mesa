define clean-comp
$(shell perl -e 'use Cwd; use File::Spec; print File::Spec->abs2rel(Cwd::realpath("$1"))' )
endef
