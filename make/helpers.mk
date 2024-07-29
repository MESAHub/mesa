define clean-comp
$(shell realpath -s --relative-to=$${PWD} $1)
endef
