include ../common/common.mk

CONTRIBS = beo-genpbs.pl prunedata.pl quick_binary_plot.sh \
	quick_binary_plot_noe.sh quick_cluster_plot.sh \
	quick_energy_plot.sh snap2data.pl

# default super-target
all: $(CONTRIBS)

install: $(CONTRIBS)
	mkdir -p $(PREFIX)/bin/
	install -m 0755 $^ $(PREFIX)/bin/

clean:
