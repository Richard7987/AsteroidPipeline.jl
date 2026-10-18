{
  description = "AsteroidPipeline.jl — and asteroidpipeline, its terminal app";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
        py = pkgs.python3Packages;

        # The terminal app (pyproject.toml), with Julia on its PATH instead
        # of juliapkg (not in nixpkgs; only needed where Julia is missing).
        app = py.buildPythonApplication {
          pname = "asteroidpipeline";
          version = (builtins.fromTOML (builtins.readFile ./Project.toml)).version;
          pyproject = true;
          src = ./.;
          build-system = [ py.hatchling ];
          dependencies = [ py.textual py.textual-image py.pillow py.pyfiglet ];
          pythonRemoveDeps = [ "juliapkg" ];
          makeWrapperArgs = [
            "--prefix" "PATH" ":" (pkgs.lib.makeBinPath [ pkgs.julia-bin pkgs.wl-clipboard pkgs.xclip ])
          ];
          doCheck = false;
          pythonImportsCheck = [ "asteroid_ui" ];
          meta.mainProgram = "asteroidpipeline";
        };
      in
      {
        packages.default = app;
        packages.asteroidpipeline = app;
        apps.default = { type = "app"; program = "${app}/bin/asteroidpipeline"; };

        devShells.default = pkgs.mkShell {
          packages = [
            pkgs.julia-bin
            (pkgs.python3.withPackages (p: [ p.textual p.textual-image p.pillow p.pyfiglet p.build p.hatchling ]))
          ];

          shellHook = ''
            export JULIA_PROJECT="$PWD"
            export PYTHONPATH="$PWD/ui''${PYTHONPATH:+:$PYTHONPATH}"
          '';
        };
      });
}
