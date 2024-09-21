{
  description = "gjk";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable"; # or whatever vers
  };

  outputs = { self, nixpkgs, ... }@inputs:
    let
      system = "aarch64-darwin"; # your version
      pkgs = nixpkgs.legacyPackages.${system};
    in {
      devShells.${system}.default = pkgs.mkShell {
        packages = with pkgs; [
          zig
          zls
          (python3.withPackages (python-pkgs: [
            python-pkgs.hypothesis
            python-pkgs.pytest
            python-pkgs.numpy
          ]))
        ];

        shellHook = ''
          unset NIX_CFLAGS_COMPILE
          unset NIX_LDFLAGS
        '';
      };
    };
}
