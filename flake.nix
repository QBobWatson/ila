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
        poppler-qt5
        pycairo
        pyyaml
      ]));

    in {
      devShell.x86_64-linux = pkgs.mkShell {

        buildInputs = with pkgs; [

          inkscape
          libxml2
          nodePackages.npm
          nodePackages.coffee-script
          scons

          python-pkgs
        ];

        packages = with pkgs; [
          fontforge
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
