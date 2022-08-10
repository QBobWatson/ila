
# Developer's Getting Started Guide

Overview:
* [Fetch the repositories](#fetch-the-repositories)
* [Install nixpkgs](#install-nixpkgs)
* [Build the site](#build-the-site)
* [Editing XML](#editing-xml)
* [Resources](#resources)

The command-line arguments given in this guide assume you're running a Linux distribution, MacOS, or [WSL2](https://docs.microsoft.com/en-us/windows/wsl/install) on a Windows machine.  (See [Install nixpkgs](#install-nixpkgs) below.)


## Fetch the repositories

First decide where you want to put all the files.  I have my files in `~/projects/`.  Change to that directory and run:
```
~/projects$ git clone https://github.com/QBobWatson/ila.git
~/projects$ cd ila
~/projects/ila$ git submodule update --init --recursive
```
This will put this repository in `ila/`, and will clone the submodules `mathbook`, `mathbook-assets`, and `mathbox`, which contain support files needed to build the book.


## Install nixpkgs

The build system has a large number of dependencies, and is sensitive to versioning.  It also requires a patched version of Inkscape.  In order to distribute a reproducible build environment, ILA is packaged as a [Nix flake](https://nixos.wiki/wiki/Flakes).  In order to use nix flakes, you first need to install [nixpkgs](https://nixos.org/download.html).  This works natively on Linux and MacOS, and also works under Windows using [WSL2](https://docs.microsoft.com/en-us/windows/wsl/install).  Follow the installation instructions linked above; it is not important whether you use a single-user or multi-user installation.

Now enable nix flakes using the instructions [here](https://nixos.wiki/wiki/Flakes#Non-NixOS).  Test the build environment using:
```
~/projects/ila$ nix develop
```
The first time you run this command, you will have to wait a long time as it pulls in several gigabytes of dependencies and rebuilds Inkscape.  Please be patient.  If all goes well, you'll be dropped into a shell from which you can build the book.


## Build the site

All build commands are run in the development shell obtained by running `nix develop`.  Alternatively, a build command, for example `scons --production`, can be run from outside the development shell by running `nix develop -c scons --production`.

This project uses the build system [scons](https://scons.org).  Before building the first time, build the dependencies:
```
~/projects/ila$ scons subpackages
```
You only need to rebuild the dependencies when they are updated.  Now you can build the book:
```
~/projects/ila$ scons
```
Beware that the book can take a very long time to build the first time, especially on a laptop: first some 4,000 html files are generated, then each one has to be preprocessed to insert math symbols, figures, etc.  If all goes well, the result of the build can be found in `ila/build`.  In order to view the output, start a web server in the build directory with:
```
~/projects/ila$ (cd build && python3 -m http.server)
```
Now point your browser at `http://localhost:8000/` to see the results of the build.

The build system accepts several options, which you pass to `scons`:
* `--build-pdf` Build the pdf version of the book in addition to the html.
* `--minify` Minify generated css and js files (for production builds).
* `--scratch` Empty the build directory before building (in case some files disappeared from the build; for production builds).
* `--production` Synonym for `--build-pdf --minify --scratch`.
* `--delete-cache` The build system maintains a cache from previous builds in `ila/cache`.  Among other things, this contains hundreds of megabytes of cached LaTeX output converted to `svg` format.  Use this option to regenerate the cache (this takes a long time).  Implies `--scratch`.
* `--theme THEME` Build with a specified visual theme.  Valid options for `THEME` are `gt` and `duke`.  Defaults to `duke`.


## Editing XML

The book is written in Robert Beezer's [PreTeXt](https://pretextbook.org), formerly known as Mathbook XML.  The source for the book is contained in the `xml` files in `ila/src`.  See the resources below for documentation on the source format.  Note that I've made many modifications to stock PreTeXt, so you'll want to poke around the source to see what's possible.  Importantly, all LaTeX code actually gets compiled by `pdflatex`, so you're not limited by what MathJax supports.

It will probably save you time in the long run to obtain and learn to use a good XML editor.  Beezer recommends [XML Copy Editor](http://xml-copy-editor.sourceforge.net/).  One advantage of a smart XML editor is that it knows what tags are allowed where.  I've adapted Beezer's XML schema files for this purpose; they are contained in `ila/mathbook/build` (after building dependencies with `scons subpackages`).  Use `pretext.rnc`, `pretext.rng`, or `pretext.xsd`, in that order of preference, depending on what kind of schema your editor supports.  You'll have to figure out how to tell your editor to use those schemas.

I use Emacs's `nxml-mode`, which came pre-installed with my Emacs distribution.  I don't think it's worth learning to use Emacs just to edit xml.  If you already use Emacs, then be sure to open the xml files in `nxml-mode`; it should automatically read the appropriate schema file from the file `ila/src/schemas.xml`.

Note that the build system will refuse to compile malformed `xml` files.  Instead it will exit with a completely useless error message from `xmllint`.  This is another good reason to use a smart xml editor.  See Beezer's documentation below for a description of the PreTeXt format, and look at `ila/mathbook/schema/pretext.xml` for the definitive reference with my additions.

## Resources

* Beezer's documentation on Mathbook XML:
    https://pretextbook.org
* Author's guide:
    https://pretextbook.org/doc/author-guide/html/
* FCLA source:
    https://github.com/rbeezer/fcla


