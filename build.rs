use std::env;
use std::path::PathBuf;

fn main() {
    // Canonicalize the path as `rustc-link-search` requires an absolute path.
    let lib_dir_path = PathBuf::from("lammps/build")
        .canonicalize()
        .expect("cannot canonicalize path");
    // Tell cargo to look for libraries in the specified directory
    println!("cargo:rustc-link-search={}", lib_dir_path.display());

    // Tell cargo to tell rustc to link the `lammps` library. Cargo will
    // automatically know it must look for a `liblammps.a` file.
    println!("cargo:rustc-link-lib=static=lammps");
    println!("cargo:rustc-link-lib=dylib=stdc++");
    println!("cargo:rustc-link-lib=dylib=mpi");
    println!("cargo:rustc-link-lib=dylib=gomp");
    println!("cargo:rustc-link-lib=dylib=png");
    println!("cargo:rustc-link-lib=dylib=jpeg");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("lammps/build/includes/lammps/library.h")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var_os("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
