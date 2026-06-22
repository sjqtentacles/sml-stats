# sml-stats build
#
#   make            build the test binary with MLton (default)
#   make test       build + run tests under MLton
#   make test-poly  run tests under Poly/ML (use-and-run; no link step)
#   make all-tests  run the suite under both compilers
#   make example    build + run examples/demo.sml
#   make clean      remove build artifacts
#
# Layout B (dependent): own sources live in src/; sml-prng and sml-specfun are
# vendored under lib/ and loaded first, in dependency order.

MLTON      ?= mlton
POLY       ?= poly
BIN        := bin
PRNGDIR    := lib/github.com/sjqtentacles/sml-prng
SPECFUNDIR := lib/github.com/sjqtentacles/sml-specfun
TEST_MLB   := test/sources.mlb
SRCS       := $(wildcard $(PRNGDIR)/* $(SPECFUNDIR)/* src/* test/*.sml) $(TEST_MLB)

.PHONY: all test poly test-poly all-tests example clean

all: $(BIN)/test-mlton

example: $(BIN)/demo
	./$(BIN)/demo

$(BIN)/demo: $(SRCS) examples/demo.sml examples/sources.mlb | $(BIN)
	$(MLTON) -output $@ examples/sources.mlb

$(BIN)/test-mlton: $(SRCS) | $(BIN)
	$(MLTON) -output $@ $(TEST_MLB)

test: $(BIN)/test-mlton
	$(BIN)/test-mlton

# Poly/ML has no native .mlb support; the suite runs at top level and exits on
# its own. Load the vendored sml-prng first, then the stats sources, then the
# test driver, in dependency order.
poly test-poly:
	printf 'use "$(PRNGDIR)/prng.sig";\nuse "$(PRNGDIR)/prng.sml";\nuse "$(SPECFUNDIR)/specfun.sig";\nuse "$(SPECFUNDIR)/specfun.sml";\nuse "src/stats.sig";\nuse "src/stats.sml";\nuse "test/harness.sml";\nuse "test/support.sml";\nuse "test/test_descriptive.sml";\nuse "test/test_distributions.sml";\nuse "test/test_regression.sml";\nuse "test/test_ttest.sml";\nuse "test/test_correlation.sml";\nuse "test/test_chisquare.sml";\nuse "test/test_ftest.sml";\nuse "test/entry.sml";\nuse "test/main.sml";\n' | $(POLY) -q --error-exit

all-tests: test test-poly

$(BIN):
	mkdir -p $(BIN)

clean:
	rm -f $(BIN)/test-mlton $(BIN)/demo
