use std::ptr;

mod lammps;

fn main() {
    println!("Hello, world!");
    unsafe {
        let simulation = lammps::lammps_open_no_mpi(0, ptr::null_mut(), ptr::null_mut());
        lammps::lammps_close(simulation);
    }
}
