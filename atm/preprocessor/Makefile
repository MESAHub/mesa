bt-settl-create ck-read plot table-create table-merge tau100-create wd-tau25-create:
	$(MAKE) --no-print-directory -f $@.mk

all: ck-read plot table-create table-merge tau100-create wd-tau25-create

clean:
	$(MAKE) --no-print-directory -f ck-read.mk clean
	$(MAKE) --no-print-directory -f table-create.mk clean
	$(MAKE) --no-print-directory -f table-merge.mk clean
	$(MAKE) --no-print-directory -f tau100-create.mk clean
	$(MAKE) --no-print-directory -f wd-tau25-create.mk clean
	$(MAKE) --no-print-directory -f plot.mk clean
	$(MAKE) --no-print-directory -f bt-settl-create.mk clean

.DEFAULT_GOAL = all
.PHONY: all clean gen-data bt-settl-create ck-read plot table-create table-merge tau100-create wd-tau25-create
