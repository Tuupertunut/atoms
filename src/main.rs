use kiss3d::{
    light::Light,
    nalgebra::{self, Point3, Translation3},
    window::Window,
};
use lammps::Lammps;
use std::{iter, slice};

mod lammps;

#[kiss3d::main]
async fn main() {
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

    let box_corners = simulation.extract_box();
    println!("{:?}", box_corners);
    let n_atoms = simulation.get_natoms();
    println!("{}", n_atoms);

    unsafe {
        let positions = slice::from_raw_parts(
            simulation.extract_atom("x") as *const *const [f64; 3],
            n_atoms,
        );
        let velocities = slice::from_raw_parts(
            simulation.extract_atom("v") as *const *const [f64; 3],
            n_atoms,
        );
        let atom_types =
            slice::from_raw_parts(simulation.extract_atom("type") as *const i32, n_atoms);

        for i in 0..n_atoms {
            println!("{}", atom_types[i]);
            println!("{:?}", *positions[i]);
            println!("{:?}", *velocities[i]);
        }
    }

    let mut window = Window::new_with_size("Atoms", 1000, 800);
    window.set_light(Light::StickToCamera);

    let mut box_cuboid = window.add_cube(0., 0., 0.);
    box_cuboid.set_surface_rendering_activation(false);
    box_cuboid.set_lines_width(1.);

    let mut atom_spheres = iter::repeat_with(|| window.add_sphere(1.))
        .take(n_atoms)
        .collect::<Vec<_>>();

    while window.render().await {
        simulation.command("run 50");

        let (box_low, box_high) = simulation.extract_box();
        let (box_low, box_high) = (Point3::from(box_low), Point3::from(box_high));

        let scale = (box_high - box_low).cast::<f32>();
        box_cuboid.set_local_scale(scale.x, scale.y, scale.z);

        box_cuboid.set_local_translation(
            Translation3::from(nalgebra::center(&box_low, &box_high)).cast::<f32>(),
        );

        let positions = unsafe {
            slice::from_raw_parts(
                simulation.extract_atom("x") as *const *const [f64; 3],
                n_atoms,
            )
            .iter()
            .map(|&ptr| *ptr)
        };

        for (atom_sphere, atom_pos) in atom_spheres.iter_mut().zip(positions) {
            atom_sphere.set_local_translation(Translation3::from(atom_pos).cast::<f32>());
        }
    }
}
