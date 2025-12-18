use std::{env, path::PathBuf};

fn main() {
    // Build lammps with cmake
    let lammps_build_dir = cmake::Config::new("lammps/cmake")
        .define("BUILD_MPI", "no")
        // TODO: is OpenMP needed?
        .define("BUILD_OMP", "yes")
        .define("WITH_JPEG", "no")
        .define("WITH_PNG", "no")
        .define("WITH_FFMPEG", "no")
        .define("WITH_GZIP", "no")
        .define("WITH_CURL", "no")
        .build();

    // Tell cargo to look for libraries in the specified directory
    println!(
        "cargo:rustc-link-search={}",
        lammps_build_dir.join("build").display()
    );

    // Tell cargo to tell rustc to link the `lammps` library. Cargo will
    // automatically know it must look for a `liblammps.a` file.
    println!("cargo:rustc-link-lib=static=lammps");
    println!("cargo:rustc-link-lib=static=gomp");
    println!("cargo:rustc-link-lib=static=stdc++");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header(
            lammps_build_dir
                .join("build/includes/lammps/library.h")
                .to_str()
                .unwrap(),
        )
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_dir = PathBuf::from(env::var_os("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_dir.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
