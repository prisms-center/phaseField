{
  lib,
  stdenv,
  fetchFromGitHub,
  cmake,
  ninja,
  python3,
  rev ? "v2.15.0",
  hash ? "sha256-O5CmG9TPPbe5TwIhRwMw5f6USa/wOSn23XEKv9PkYjc=",
  kokkos,
  mpi,
  debug ? false,
}:
stdenv.mkDerivation {
  pname = "caliper";
  version = rev;

  src = fetchFromGitHub {
    owner = "LLNL";
    repo = "Caliper";
    inherit rev hash;
  };

  nativeBuildInputs = [cmake ninja python3];

  cmakeBuildType = lib.optional debug "Debug";

  buildInputs = [kokkos mpi];

  cmakeFlags = [
    (lib.cmakeBool "WITH_MPI" true)
  ];
}
