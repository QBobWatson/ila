
# Developer's Getting Started Guide

The command-line arguments given in this guide assume you're using a UNIX-like system (e.g., a Mac).  Building under Windows should be similar.

Overview: 
* Fetch the repositories
* Setup the build system
* Build the site
* Editing XML
* Resources


## Fetch the repositories

First decide where you want to put all the files.  I have my files in `~/projects/`.  Change to that directory and run:
```
~/projects$ git clone https://github.com/QBobWatson/ila.git
~/projects$ cd ila
~/projects/ila$ git submodule update --init --recusrive
```
This will put this repository in `ila/`, and will clone the submodules `mathbook`, `mathbook-assets`, and `mathbox`, which contain support files needed to build the book.

## Setup the build system

The build system has a large number of dependencies, and is sensitive to versioning.  It also requires a patched version of Inkscape.  For this reason, I've packaged everything needed to build the book into a [Vagrant](https://www.vagrantup.com/) box.  This is a prepackaged virtual machine that can be launched from any Unix, Mac, or Windows computer.  It has two prerequisites:
* VirtualBox: the underlying virtual machine software.  [Download](https://www.virtualbox.org/wiki/Downloads).
* Vagrant: the program that manages the virtual machine.  [Download](https://www.vagrantup.com/downloads.html).

(If you know what you are doing, you can also use the `libvirt` provider, which is faster on a Linux system.)

To install the build environment, change into `ila/build-environment`, and type:
```
ila/build-environment$ vagrant up --provider virtualbox
```
The box will now configure itself.  It installs an enormous number of packages, including all of `TeXLive`, so be patient.  You may want to tweak `Vagrantfile` to configure how much memory and CPU resources to give to the virtual machine.

If you want to poke around the virtual machine, use `vagrant ssh`.  To stop it, type `vagrant halt`.  All `vagrant` commands must be run from the `ila/build-environment` directory.

## Build the site

All build commands are run inside the Vagrant virtual machine.  Before building, run `vagrant ssh` from `ila/build-environment`; this will give you a prompt inside the virtual machine.  Then `cd` into the `/base` directory, which mirrors the `ila` directory on your home system.

This project uses the build system [scons](https://scons.org).  Before building the first time, build the dependencies:
```
/base$ scons subpackages
```
You only need to rebuild the dependencies when they are updated.  Now you can build the book:
```
/base$ scons
```
Beware that the book can take a very long time to build the first time, especially on a laptop: first some 4,000 html files are generated, then each one has to be preprocessed to insert math symbols, figures, etc.  If all goes well, the result of the build can be found in `/base/build`, which is `ila/build` on your host machine.  The virtual machine runs a webserver; point your browser at `http://localhost:8081/` to see the results of the build.  (Opening the file `build/index.html` directly in a browser doesn't work with knowls for technical reasons.)

The build system maintains a cache from previous builds in `ila/staging`.  Among other things, this contains hundreds of megabytes of cached LaTeX output converted to `svg` format.  To build from scratch again, simply `rm -r ila/staging`.

The build system accepts several options, which you pass to `scons`:
* `--production` Makes a version suitable for publishing.  The main difference is minification of CSS and Javascript files.
* `--build-pdf` Also build the pdf version of the book.  Implied by `--production`.
* `--theme THEME` Build with a specified visual theme.  Valid options for `THEME` are `gt` and `duke`.  Defaults to `gt`.
* `--variant VARIANT` Build a variant of the book.  Currently the only variant is `1553`.  The default variant is `default`.

## Editing XML

The book is written in Robert Beezer's [PreTeXt](https://pretextbook.org), formerly known as Mathbook XML.  The source for the book is contained in the `xml` files in `ila/src`.  See the resources below for documentation on the source format.  Note that I've made many modifications to stock PreTeXt, so you'll want to poke around the source to see what's possible.  Importantly, all LaTeX code actually gets compiled by `pdflatex`, so you're not limited by what MathJax supports.

It will probably save you time in the long run to obtain and learn to use a good XML editor.  Beezer recommends [XML Copy Editor](http://xml-copy-editor.sourceforge.net/).  One advantage of a smart XML editor is that it knows what tags are allowed where.  I've adapted Beezer's XML schema files for this purpose; they are contained in `ila/mathbook/build` (after building dependencies with `scons subpackages`).  Use `pretext.rnc`, `pretext.rng`, or `pretext.xsd`, in that order of preference, depending on what kind of schema your editor supports.  You'll have to figure out how to tell your editor to use those schemas.

I use Emacs's `nxml-mode`, which came pre-installed with my Emacs distribution.  I don't think it's worth learning to use Emacs just to edit xml.  If you already use Emacs, then be sure to open the xml files in `nxml-mode`; it should automatically read the appropriate schema file from the file `ila/src/schemas.xml`.

Note that the build system will refuse to compile malformed `xml` files.  See Beezer's documentation below for a description of the PreTeXt format, and look at `ila/mathbook/schema/pretext.xml` for the definitive reference with my additions.

## Resources

* Beezer's documentation on Mathbook XML:
    https://pretextbook.org
* Author's guide:
    https://pretextbook.org/doc/author-guide/html/
* FCLA source:
    https://github.com/rbeezer/fcla


