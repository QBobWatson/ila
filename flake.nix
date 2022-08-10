{
  description = "Interactive Linear Algebra";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-22.05";
  };

  outputs = { self, nixpkgs }@inputs:
    let
      pkgs = import nixpkgs {
        overlays = [
          (self: super: {
            inkscape = super.inkscape.overrideAttrs (old: {
              patches = (old.patches or []) ++ [ ./build-environment/inkscape.patch ];
            });

            # Need to use updated Gemfile.  Presumably this will be fixed
            # in nixpkgs at some point; then ./build-environment/compass can be
            # deleted.
            compass = super.bundlerApp {
              pname = "compass";
              gemdir = ./build-environment/compass;
              exes = [ "compass" ];

              passthru.updateScript = super.bundlerUpdateScript "compass";
            };
          })
        ];
        system = "x86_64-linux";
      };

      python-pkgs = pkgs.python3.withPackages (p: (with p; [
        beautifulsoup4
        cssutils
        fontforge
        lxml
        Mako
        pdfrw
        pycairo
        pygobject3
        pyyaml
      ]));

    in {
      devShell.x86_64-linux = pkgs.mkShell {

        buildInputs = with pkgs; [

          compass
          inkscape
          libxml2
          nodejs
          nodePackages.npm
          nodePackages.coffee-script
          poppler_gi
          scons

          python-pkgs
        ];

        packages = with pkgs; [
          fontforge
          git
          jing-trang
          libxslt
          texlive.combined.scheme-full
        ];
      };

      shellHook = ''
        PYTHONPATH=${python-pkgs}/${python-pkgs.sitePackages}:$PYTHONPATH
      '';
    };
}
