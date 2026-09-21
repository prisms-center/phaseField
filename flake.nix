{
  description = "C++ partial differential equation framework";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-parts.url = "github:hercules-ci/flake-parts";

    dealii.url = "git+https://codeberg.org/landinjm/dealii-flake.git";
  };

  outputs = {flake-parts, ...} @ inputs:
    flake-parts.lib.mkFlake {inherit inputs;} {
      systems = [
        "aarch64-darwin"
        "aarch64-linux"
        "x86_64-darwin"
        "x86_64-linux"
      ];

      perSystem = {
        pkgs,
        system,
        ...
      }: let
        config = {
          allowUnfree = true;
          cpuArch = "NATIVE";
          cudaSupport = false;
          cudaArch = "ADA89";
          rocmSupport = false;
          rocmArch = "";
        };
      in {
        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            # Main requirements
            cmake
            gnumake
            ninja
            gcc
            inputs.dealii.packages.${system}.default

            # Pre-commit
            llvmPackages_18.clang-tools
            pre-commit

            # Documentation
            doxygen
            graphviz
          ];
        };
      };
    };
}
