use std::{ffi::CString, ptr, slice};

mod lammps;

fn main() {
    println!("Hello, world!");
    unsafe {
        let simulation = lammps::lammps_open_no_mpi(0, ptr::null_mut(), ptr::null_mut());

        lammps::lammps_command(simulation, CString::new("units metal").unwrap().as_ptr());
        lammps::lammps_command(simulation, CString::new("dimension 3").unwrap().as_ptr());
        lammps::lammps_command(simulation, CString::new("boundary p p p").unwrap().as_ptr());
        lammps::lammps_command(
            simulation,
            CString::new("atom_style atomic").unwrap().as_ptr(),
        );
        lammps::lammps_command(
            simulation,
            CString::new("region box block 0 20 0 20 0 20")
                .unwrap()
                .as_ptr(),
        );
        lammps::lammps_command(
            simulation,
            CString::new("create_box 1 box").unwrap().as_ptr(),
        );
        lammps::lammps_command(
            simulation,
            CString::new("create_atoms 1 single 10 10 10")
                .unwrap()
                .as_ptr(),
        );
        lammps::lammps_command(
            simulation,
            CString::new("create_atoms 1 single 15 10 10")
                .unwrap()
                .as_ptr(),
        );
        lammps::lammps_command(simulation, CString::new("mass 1 40").unwrap().as_ptr());
        lammps::lammps_command(
            simulation,
            CString::new("pair_style lj/cut 7").unwrap().as_ptr(),
        );
        lammps::lammps_command(
            simulation,
            CString::new("pair_coeff 1 1 0.01 3.3").unwrap().as_ptr(),
        );
        lammps::lammps_command(simulation, CString::new("fix 1 all nve").unwrap().as_ptr());
        lammps::lammps_command(simulation, CString::new("timestep 0.001").unwrap().as_ptr());
        lammps::lammps_command(simulation, CString::new("run 1000").unwrap().as_ptr());

        let n_atoms = lammps::lammps_get_natoms(simulation) as usize;
        println!("{}", n_atoms);

        let coords = lammps::lammps_extract_atom(simulation, CString::new("x").unwrap().as_ptr())
            as *const *const [f64; 3];
        for &atom_coords in slice::from_raw_parts(coords, n_atoms).into_iter() {
            println!("{:?}", *atom_coords);
        }

        lammps::lammps_close(simulation);
    }
}
