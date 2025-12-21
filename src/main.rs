use lammps::Lammps;
use std::slice;

mod lammps;

fn main() {
    println!("Hello, world!");
    let mut simulation = Lammps::open();

    simulation.command("units metal");
    simulation.command("dimension 3");
    simulation.command("boundary p p p");
    simulation.command("atom_style atomic");
    simulation.command("region box block 0 20 0 20 0 20");
    simulation.command("create_box 1 box");
    simulation.command("create_atoms 1 single 10 10 10");
    simulation.command("create_atoms 1 single 15 10 10");
    simulation.command("mass 1 40");
    simulation.command("pair_style lj/cut 7");
    simulation.command("pair_coeff 1 1 0.01 3.3");
    simulation.command("fix 1 all nve");
    simulation.command("timestep 0.001");
    simulation.command("run 1000");

    let n_atoms = simulation.get_natoms();
    println!("{}", n_atoms);

    let coords = simulation.extract_atom("x") as *const *const [f64; 3];
    unsafe {
        for &atom_coords in slice::from_raw_parts(coords, n_atoms).into_iter() {
            println!("{:?}", *atom_coords);
        }
    }
}
