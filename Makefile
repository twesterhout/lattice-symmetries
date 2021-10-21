PREFIX = $(PWD)/prefix
PROJECT_NAME = lattice-symmetries-haskell
LIBRARY_NAME = $(subst -,_,$(PROJECT_NAME))

GHC_VERSION ?= $(shell ghc --numeric-version)
CABAL_BUILD_DIR = $(shell dirname $$(find dist-newstyle/build -name 'libHS$(PROJECT_NAME)*.a'))
CABAL_AUTOGEN_DIR = $(CABAL_BUILD_DIR)/global-autogen

HS_LDFLAGS = $(shell cat "$(CABAL_AUTOGEN_DIR)/HS_LIBRARY_PATHS_LIST" | sed 's/^/-L"/;s/$$/"/' | tr '\n' ' ')
HS_LDFLAGS += $(shell cat "$(CABAL_AUTOGEN_DIR)/HS_LIBRARIES_LIST" | sed 's/^/-l/' | tr '\n' ' ')
C_LDFLAGS = $(shell cat "$(CABAL_AUTOGEN_DIR)/EXTRA_LIBRARIES_LIST" | sed 's/^/-l/' | tr '\n' ' ')

$(shell mkdir -p build)

all: shared static main

shared: build/lib$(LIBRARY_NAME).so
static: build/lib$(LIBRARY_NAME).a

$(CABAL_BUILD_DIR)/cbits/init.o:
	cabal v2-build

build/api.txt: cbits/init.c
	cat $< \
		| tr '\n' ' ' \
		| sed -E 's/.*ls_hs_symbol_table\s*\[\]\s*=\s*\{([^}]*)\};.*/\1/' \
		| tr -d '& ' \
		| tr ',' '\n' \
		> $@

build/lib$(LIBRARY_NAME).so: $(CABAL_BUILD_DIR)/cbits/init.o build/api.txt $(CABAL_AUTOGEN_DIR)/HS_LIBRARY_PATHS_LIST $(CABAL_AUTOGEN_DIR)/HS_LIBRARIES_LIST $(CABAL_AUTOGEN_DIR)/EXTRA_LIBRARIES_LIST
	ghc --make -no-hs-main -shared \
		-pgmc $(CC) -pgml $(CC) \
		-optl -Wl,--retain-symbols-file=build/api.txt \
		$< -o $@ \
		-L"$(CABAL_BUILD_DIR)" \
		$(HS_LDFLAGS) $(C_LDFLAGS)

build/lib$(LIBRARY_NAME).a: $(CABAL_BUILD_DIR)/cbits/init.o $(CABAL_AUTOGEN_DIR)/HS_LIBRARY_PATHS_LIST $(CABAL_AUTOGEN_DIR)/HS_LIBRARIES_LIST $(CABAL_AUTOGEN_DIR)/EXTRA_LIBRARIES_LIST
	$(LD) -o $@ --relocatable --whole-archive -L"$(CABAL_BUILD_DIR)" $(HS_LDFLAGS)

main: main.c
	gcc -o $@ -I cbits $< -L build -l$(LIBRARY_NAME)

clean:
	rm -rf main build/
