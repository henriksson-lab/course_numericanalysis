.PHONY: all render preview clean gitadd help figs demos typecheck edit-slides

# Default: render the site, then stage every source file in git.
all: render gitadd

render: demos
	quarto render
	@touch docs/.nojekyll

preview:
	quarto preview

# Re-generate script-made figures (not run by quarto render).
figs:
	@for f in */make-*-figs.R */make-*-figs.py; do \
	  [ -e "$$f" ] || continue; \
	  echo "$$f"; d=$$(dirname $$f); b=$$(basename $$f); \
	  case $$b in *.R) (cd $$d && Rscript $$b);; *.py) (cd $$d && python3 $$b);; esac; \
	done

# Build the interactive demos (TypeScript, demos/src/) into demos/dist/.
# One-time setup: npm i   (installs esbuild + typescript into node_modules/).
demos:
	demos/build.sh

# Type-check the demo sources (esbuild only strips types, it does not check them).
typecheck:
	node_modules/.bin/tsc --noEmit -p demos

# Edit slide text and SVG labels in the browser (tools/slide-editor/).
# Prints a URL with an access token; PORT=... and SLIDE_EDIT_TOKEN=... override.
PORT ?= 8737
edit-slides:
	node tools/slide-editor/server.mjs --port $(PORT)

# Remove the rendered site and Quarto's cache.
clean:
	rm -rf docs .quarto

# Stage every source file AND the built HTML (docs/) so the committed site is
# always current. The .quarto/ cache is never added (it is .gitignored).
gitadd:
	git add -- .

help:
	@echo "Targets:"
	@echo "  make           render the site, then git add all source files"
	@echo "  make render    render the whole site to docs/"
	@echo "  make preview   start a live-reload preview server"
	@echo "  make demos     build the interactive demos into demos/dist/ (needs npm i once)"
	@echo "  make typecheck type-check the demo TypeScript sources"
	@echo "  make figs      re-generate figures made by <unit>/make-*-figs.{R,py}"
	@echo "  make edit-slides  edit slide text and figure labels in the browser (port 8737)"
	@echo "  make gitadd    stage all source files in git"
	@echo "  make clean     remove docs/ and .quarto/"
