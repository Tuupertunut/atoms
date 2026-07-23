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
use std::slice;

mod lammps;

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
    simulation.command("units metal");
    simulation.command("dimension 3");
    simulation.command("boundary p p p");
    simulation.command("atom_style atomic");
    simulation.command("region box block 0 20 0 20 0 20");
    simulation.command("create_box 1 box");
    // Using argon mass and Lennard-Jones parameters
    // https://www.researchgate.net/figure/Lennard-Jones-LJ-potential-parameters-of-different-materials-considered-in-thepresent_tbl2_319412425
    simulation.command("mass 1 40");
    simulation.command("pair_style lj/cut 15");
    simulation.command("pair_coeff 1 1 0.01 3.4");
    simulation.command("timestep 0.001");

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
    let mut template_sphere = scene.add_sphere(1.);
    template_sphere.set_surface_rendering_activation(false);
    template_sphere.set_lines_width(0.5, false);

    let mut template_add_mode = false;
    let mut template_distance = 30.;
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
                            simulation.create_atom(1, template_pos, [0., 0., 0.]);
                            let mut atom_sphere = scene.add_sphere(1.);
                            atom_sphere.set_color(Color::new(
                                0x80 as f32 / 256.,
                                0xD1 as f32 / 256.,
                                0xE3 as f32 / 256.,
                                1.,
                            ));
                            atom_spheres.push(atom_sphere);

                            template_add_mode = false;
                        } else {
                            if let Some(selected_pos) = selection {
                                // Delete atom
                                // Lammps does not directly allow deleting a single atom, so create
                                // a temporary region around it and delete everything inside it. In
                                // practice the region is so small that there is always only that
                                // one atom in it.
                                simulation.command(&format!(
                                    "region temp sphere {} {} {} 0.01",
                                    selected_pos.x, selected_pos.y, selected_pos.z
                                ));
                                simulation.command("delete_atoms region temp");
                                simulation.command("region temp delete");
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
                    "fix 1 all npt temp {} {} 0.1 iso {} {} 1",
                    thermostat_temperature,
                    thermostat_temperature,
                    barostat_pressure,
                    barostat_pressure
                ));
            } else if thermostat_enabled {
                simulation.command(&format!(
                    "fix 1 all nvt temp {} {} 0.1",
                    thermostat_temperature, thermostat_temperature
                ));
            } else if barostat_enabled {
                simulation.command(&format!(
                    "fix 1 all nph iso {} {} 1",
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
        let positions = unsafe {
            slice::from_raw_parts(
                simulation.extract_atom("x") as *const *const [f64; 3],
                atom_spheres.len(),
            )
            .iter()
            .map(|&ptr| *ptr)
        };

        for (atom_sphere, atom_pos) in atom_spheres.iter_mut().zip(positions) {
            atom_sphere.set_position(DVec3::from(atom_pos).as_vec3());
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
            } else {
                selection = atom_spheres
                    .iter()
                    .filter_map(|atom_sphere| {
                        Ball::new(atom_sphere.local_scale().x / 2.)
                            .cast_ray(
                                &atom_sphere.local_transformation(),
                                &cursor_ray,
                                camera.clip_planes().1,
                                true,
                            )
                            .map(|distance| (atom_sphere.position(), distance))
                    })
                    .min_by(|(_, dist_a), (_, dist_b)| dist_a.partial_cmp(dist_b).unwrap())
                    .map(|(selected_pos, _)| selected_pos);

                match selection {
                    Some(selected_pos) => {
                        template_sphere.set_visible(true);
                        template_sphere.set_lines_color(Some(color::RED));
                        template_sphere.set_position(selected_pos);
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
