{
  description = "C++ partial differential equation framework";

  inputs = {
    # NOTE: Using 26.05 to support intel macs
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-26.05";
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
      }: {
        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            pkg-config

            # Main packages
            cmake
            gnumake
            ninja
            gcc
            mpi
            inputs.dealii.packages.${system}.default

            # Optional packages
            python3
            vtk

            # Pre-commit
            llvmPackages_18.clang-tools
            pre-commit

            # Documentation
            doxygen
            graphviz
          ];

          shellHook = ''
            # Grab submodules
            git submodule update --init --recursive

            # Use half the available jobs so we don't run out of memory
            jobs=$(( $(nproc) / 2 ))
            [ "$jobs" -lt 1 ] && jobs=1
            export CMAKE_BUILD_PARALLEL_LEVEL=$jobs

            export PRISMS_PF_DIR="$PWD/install"

            alias prisms_build="cmake --workflow --preset debugrelease"
            alias prisms_test="cmake --workflow --preset test"
            alias prisms_docs="cmake --workflow --preset docs"
            alias prisms_prm_format="./contrib/utilities/prm_format.sh"
            alias prisms_copyright="./contrib/utilities/update_copyright.sh"

            prisms_serve_docs() {
              local docs_dir="build/docs/doc/$(git branch --show-current)"

              if [ ! -d "$docs_dir" ]; then
                echo "Error: Documentation directory not found: $docs_dir" >&2
                return 1
              fi

              (
                cd "$docs_dir" || exit 1
                echo "Serving docs at http://localhost:8000"
                python3 -m http.server 8000
              )
            }

            prisms_help() {
              echo "Available development commands:"
              echo "  prisms_build        Build debug and release configurations"
              echo "  prisms_test         Build and run the tests"
              echo "  prisms_docs         Build the docs"
              echo "  prisms_serve_docs   Serve docs at localhost:8000"
              echo "  prisms_prm_format   Format .prm files"
              echo "  prisms_copyright    Update copyright headers"
              echo ""
              echo "Run prisms_help to display this list again."
            }

            echo "Dev shell loaded!"
            prisms_help
          '';
        };
      };
    };
}
