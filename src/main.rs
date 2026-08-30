use kiss3d::{
    camera::{Camera3d, OrbitCamera3d},
    color::{self, Color},
    egui::{Grid, Panel, ProgressBar, Slider},
    event::{Action, Key, MouseButton, WindowEvent},
    glamx::{DVec2, DVec3, Vec3},
    light::Light,
    scene::SceneNode3d,
    window::Window,
};
use lammps::Lammps;
use parry3d::{
    query::{Ray, RayCast},
    shape::Ball,
};

mod lammps;

/// Jmol atom colors for H,C,N,O
fn atom_color(atom_type: i32) -> Color {
    let (r, g, b) = match atom_type {
        1 => (0xFF, 0xFF, 0xFF),
        2 => (0x90, 0x90, 0x90),
        3 => (0x30, 0x50, 0xF8),
        4 => (0xFF, 0x0D, 0x0D),
        _ => unreachable!(),
    };
    return Color::new(r as f32 / 256., g as f32 / 256., b as f32 / 256., 1.);
}

/// Covalent atom radii for H,C,N,O from Wikipedia table, times 1.5 because it looks good
/// https://en.wikipedia.org/wiki/Covalent_radius
fn atom_radius(atom_type: i32) -> f32 {
    let radius = match atom_type {
        1 => 0.31,
        2 => 0.76,
        3 => 0.71,
        4 => 0.66,
        _ => unreachable!(),
    };
    return radius * 1.5;
}

/// Deletes a single atom which is closest to given position. Lammps doesn't have any native way to
/// delete just one atom, only all atoms in a region, so we have to implement it ourselves using
/// atom IDs. Unsafe when given atom count is too high.
unsafe fn delete_closest_atom(simulation: &mut Lammps, atom_count: usize, pos: DVec3) {
    let ids = unsafe { simulation.extract_atom_id(atom_count) };
    let positions = unsafe { simulation.extract_atom_x(atom_count) };
    let (closest_id, _) = ids
        .zip(positions.map(|atom_pos| DVec3::from(atom_pos).distance_squared(pos)))
        .min_by(|(_, dist_a), (_, dist_b)| dist_a.partial_cmp(dist_b).unwrap())
        .unwrap();
    simulation.command(&format!("group temp id {}", closest_id));
    simulation.command("delete_atoms group temp");
    simulation.command("group temp delete");
}

#[kiss3d::main]
async fn main() {
    // Initialize window
    let mut window = Window::new_with_size("Atoms", 1000, 800).await;

    let mut camera = OrbitCamera3d::new(Vec3::new(10., 10., -20.), Vec3::new(10., 10., 10.));
    camera.set_dist_step(0.99);

    let mut scene = SceneNode3d::empty();
    let mut light = scene.add_light(Light::default());

    // Initialize simulation
    let mut simulation = Lammps::open(&["-log", "none"]);
    simulation.command("units real");
    simulation.command("dimension 3");
    simulation.command("boundary p p p");
    simulation.command("atom_style charge");
    simulation.command("region box block 0 20 0 20 0 20");
    simulation.command("create_box 4 box");
    // Using H,C,N,O mass and reaxff parameters
    simulation.command("mass 1 1.008");
    simulation.command("mass 2 12");
    simulation.command("mass 3 14");
    simulation.command("mass 4 15.999");
    // Hack: arbitrary values to prevent crashing, use kokkos instead
    simulation.command("pair_style reaxff NULL safezone 5 mincap 300 minhbonds 300");
    simulation.command("pair_coeff * * ffield.reax.chon2019 H C N O");
    simulation.command("fix 2 all qeq/reaxff 1 0 10 1e-6 reaxff");
    simulation.command("timestep 0.2");

    // Initialize simulation control parameters
    let mut simulation_running = true;

    let mut ensemble_changed = true;

    let mut thermostat_enabled = false;
    let mut thermostat_temperature = 300.;
    let mut barostat_enabled = false;
    let mut barostat_pressure = 300.;

    // Initialize simulation structures in 3D scene
    let mut box_cuboid = scene.add_cube(0., 0., 0.);
    box_cuboid.set_surface_rendering_activation(false);
    box_cuboid.set_lines_width(1., false);

    let mut atom_spheres = Vec::<SceneNode3d>::new();

    // Initialize template sphere, the temporary sphere to display when adding/deleting atoms
    let mut template_sphere = scene.add_sphere(0.);
    template_sphere.set_surface_rendering_activation(false);
    template_sphere.set_lines_width(0.5, false);

    let mut template_add_mode = false;
    let mut template_distance = 30.;
    let mut template_atom_type = 1;
    let mut selection = Option::<Vec3>::None;

    // Initialize button drag monitoring
    let mut last_button1_pressed_pos = DVec2::ZERO;
    let mut button1_pressed_without_dragging = false;
    let mut last_button2_pressed_pos = DVec2::ZERO;
    let mut button2_pressed_without_dragging = false;

    // Run the main render loop
    while window.render_3d(&mut scene, &mut camera).await {
        // Check for 3D UI interaction events
        for mut event in window.events().iter() {
            match event.value {
                // Handle 3D UI interaction
                WindowEvent::MouseButton(MouseButton::Button1, Action::Release, _)
                    if !window.is_egui_capturing_mouse() =>
                {
                    if button1_pressed_without_dragging {
                        if template_add_mode {
                            // Add new atom
                            let template_pos = template_sphere.position().as_dvec3().to_array();
                            simulation.create_atom(template_atom_type, template_pos, [0., 0., 0.]);
                            atom_spheres.push(scene.add_sphere(0.));

                            template_add_mode = false;
                        } else {
                            if let Some(selected_pos) = selection {
                                // Delete atom
                                // Safety:
                                // We always push created atoms and pop deleted atoms from
                                // atom_spheres, and lammps never changes the atom count by itself,
                                // so the length of atom_spheres is also the number of atoms in the
                                // simulation.
                                unsafe {
                                    delete_closest_atom(
                                        &mut simulation,
                                        atom_spheres.len(),
                                        selected_pos.as_dvec3(),
                                    );
                                }
                                atom_spheres.pop().unwrap().remove();
                            }
                        }
                    }

                    button1_pressed_without_dragging = false;
                }
                WindowEvent::MouseButton(MouseButton::Button2, Action::Release, _)
                    if !window.is_egui_capturing_mouse() =>
                {
                    if button2_pressed_without_dragging {
                        template_add_mode = !template_add_mode;
                    }

                    button2_pressed_without_dragging = false;
                }
                WindowEvent::Scroll(_, y_offset, _) if !window.is_egui_capturing_mouse() => {
                    if template_add_mode {
                        event.inhibited = true;
                        template_distance *= f32::powf(1.01, y_offset as f32);
                    }
                }
                WindowEvent::Key(Key::Space, Action::Press, _)
                    if !window.is_egui_capturing_keyboard() =>
                {
                    simulation_running = !simulation_running;
                }

                WindowEvent::Key(Key::Key1, Action::Press, _)
                    if !window.is_egui_capturing_keyboard() =>
                {
                    template_atom_type = 1;
                }
                WindowEvent::Key(Key::Key2, Action::Press, _)
                    if !window.is_egui_capturing_keyboard() =>
                {
                    template_atom_type = 2;
                }
                WindowEvent::Key(Key::Key3, Action::Press, _)
                    if !window.is_egui_capturing_keyboard() =>
                {
                    template_atom_type = 3;
                }
                WindowEvent::Key(Key::Key4, Action::Press, _)
                    if !window.is_egui_capturing_keyboard() =>
                {
                    template_atom_type = 4;
                }

                // Update button drag monitoring. Some actions should only happen when buttons are
                // clicked without dragging so we keep track of it.
                WindowEvent::MouseButton(MouseButton::Button1, Action::Press, _)
                    if !window.is_egui_capturing_mouse() =>
                {
                    let cursor_pos = window.cursor_pos().unwrap();
                    last_button1_pressed_pos = DVec2::from(cursor_pos);
                    button1_pressed_without_dragging = true;
                }
                WindowEvent::MouseButton(MouseButton::Button2, Action::Press, _)
                    if !window.is_egui_capturing_mouse() =>
                {
                    let cursor_pos = window.cursor_pos().unwrap();
                    last_button2_pressed_pos = DVec2::from(cursor_pos);
                    button2_pressed_without_dragging = true;
                }
                WindowEvent::CursorPos(cursor_x, cursor_y, _) => {
                    let moved_pos = DVec2::new(cursor_x, cursor_y);
                    if button1_pressed_without_dragging
                        && last_button1_pressed_pos.distance_squared(moved_pos) >= f64::powi(10., 2)
                    {
                        button1_pressed_without_dragging = false;
                    }
                    if button2_pressed_without_dragging
                        && last_button2_pressed_pos.distance_squared(moved_pos) >= f64::powi(10., 2)
                    {
                        button2_pressed_without_dragging = false;
                    }
                }

                _ => {}
            }
        }

        // Draw overlay 2D UI and check for its interaction events
        let temperature = simulation.get_thermo("temp");
        let pressure = simulation.get_thermo("press");

        window.draw_ui(|ctx| {
            Panel::bottom("bottom_panel").show(ctx, |ui| {
                Grid::new("stat_grid").num_columns(2).show(ui, |ui| {
                    let mut bar_width = 0.;

                    let checkbox = ui.checkbox(&mut thermostat_enabled, "Thermostat");
                    if checkbox.changed() {
                        ensemble_changed = true;
                    }

                    ui.vertical(|ui| {
                        // Calculating bar width here because we don't know the available space in
                        // the second column before this point.
                        bar_width = f32::max(0., ui.available_width() - 80.);

                        let max_temperature = 1200.;

                        ui.spacing_mut().slider_width = bar_width;
                        let slider = ui.add(
                            Slider::new(&mut thermostat_temperature, 0.0..=max_temperature)
                                .suffix(" K"),
                        );
                        if slider.changed() {
                            ensemble_changed = true;
                        }

                        ui.add(
                            ProgressBar::new((temperature / max_temperature) as f32)
                                .desired_width(bar_width)
                                .text(format!("{:.2} K", temperature)),
                        );
                    });

                    ui.end_row();

                    let checkbox = ui.checkbox(&mut barostat_enabled, "Barostat");
                    if checkbox.changed() {
                        ensemble_changed = true;
                    }

                    ui.vertical(|ui| {
                        let max_pressure = 1200.;

                        ui.spacing_mut().slider_width = bar_width;
                        let slider = ui.add(
                            Slider::new(&mut barostat_pressure, 0.0..=max_pressure).suffix(" bar"),
                        );
                        if slider.changed() {
                            ensemble_changed = true;
                        }

                        ui.add(
                            ProgressBar::new((pressure / max_pressure) as f32)
                                .desired_width(bar_width)
                                .text(format!("{:.2} bar", pressure)),
                        );
                    });

                    ui.end_row();
                });
            });
        });

        // Check for ensemble changes, nve, nvt, nph, npt
        if ensemble_changed {
            simulation.command("unfix 1");

            if thermostat_enabled && barostat_enabled {
                simulation.command(&format!(
                    "fix 1 all npt temp {} {} 200 iso {} {} 1000",
                    thermostat_temperature,
                    thermostat_temperature,
                    barostat_pressure,
                    barostat_pressure
                ));
            } else if thermostat_enabled {
                simulation.command(&format!(
                    "fix 1 all nvt temp {} {} 200",
                    thermostat_temperature, thermostat_temperature
                ));
            } else if barostat_enabled {
                simulation.command(&format!(
                    "fix 1 all nph iso {} {} 1000",
                    barostat_pressure, barostat_pressure
                ));
            } else {
                simulation.command("fix 1 all nve");
            }

            ensemble_changed = false;
        }

        // Run simulation one frame forward
        if simulation_running {
            simulation.command("run 50");
        }

        // Update box position in 3D
        let (box_low, box_high) = simulation.extract_box();
        let (box_low, box_high) = (
            DVec3::from(box_low).as_vec3(),
            DVec3::from(box_high).as_vec3(),
        );

        let scale = box_high - box_low;
        box_cuboid.set_local_scale(scale.x, scale.y, scale.z);

        box_cuboid.set_position(box_low.midpoint(box_high));

        // Update atom positions in 3D
        // Safety:
        // We always push created atoms and pop deleted atoms from atom_spheres, and lammps never
        // changes the atom count by itself, so the length of atom_spheres is also the number of
        // atoms in the simulation.
        let positions = unsafe { simulation.extract_atom_x(atom_spheres.len()) };
        let types = unsafe { simulation.extract_atom_type(atom_spheres.len()) };

        for (atom_sphere, (atom_pos, atom_type)) in
            atom_spheres.iter_mut().zip(positions.zip(types))
        {
            atom_sphere.set_position(DVec3::from(atom_pos).as_vec3());
            atom_sphere.set_color(atom_color(atom_type));
            let scale = atom_radius(atom_type) * 2.;
            atom_sphere.set_local_scale(scale, scale, scale);
        }

        // Update template sphere in 3D
        if let Some(cursor_pos) = window.cursor_pos() {
            let (origin, direction) =
                camera.unproject(DVec2::from(cursor_pos).as_vec2(), window.size().as_vec2());
            let cursor_ray = Ray::new(origin, direction);

            if template_add_mode {
                selection = None;

                template_sphere.set_visible(true);
                template_sphere.set_lines_color(Some(color::LIME));
                template_sphere.set_position(cursor_ray.point_at(template_distance));
                let scale = atom_radius(template_atom_type) * 2.;
                template_sphere.set_local_scale(scale, scale, scale);
            } else {
                // Ray tracing to find which sphere the cursor is pointing at
                let sphere_selection = atom_spheres
                    .iter()
                    .filter_map(|atom_sphere| {
                        Ball::new(atom_sphere.local_scale().x / 2.)
                            .cast_ray(
                                &atom_sphere.local_transformation(),
                                &cursor_ray,
                                camera.clip_planes().1,
                                true,
                            )
                            .map(|distance| (atom_sphere, distance))
                    })
                    .min_by(|(_, dist_a), (_, dist_b)| dist_a.partial_cmp(dist_b).unwrap());

                selection = sphere_selection.map(|(selected_sphere, _)| selected_sphere.position());

                match sphere_selection {
                    Some((selected_sphere, _)) => {
                        template_sphere.set_visible(true);
                        template_sphere.set_lines_color(Some(color::RED));
                        template_sphere.set_position(selected_sphere.position());
                        let scale = selected_sphere.local_scale().x;
                        template_sphere.set_local_scale(scale, scale, scale);
                    }
                    None => {
                        template_sphere.set_visible(false);
                    }
                }
            }
        }

        // Update light position in 3D
        light.set_position(camera.eye());
    }
}
