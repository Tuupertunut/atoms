use itertools::Itertools;
use kiss3d::{
    camera::{Camera3d, OrbitCamera3d},
    color,
    egui::{self, Align2, Button, Color32, FontId, Grid, Panel, ProgressBar, Slider},
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
use rgb::{ComponentMap, Rgb};
use std::io::Write;
use strum::{EnumIter, FromRepr, IntoEnumIterator};
use tempfile::{NamedTempFile, TempPath};

mod lammps;

#[derive(Clone, Copy, PartialEq, Eq, FromRepr, EnumIter)]
enum AtomType {
    H = 1,
    C = 2,
    N = 3,
    O = 4,
}

impl AtomType {
    fn symbol(self) -> &'static str {
        return match self {
            Self::H => "H",
            Self::C => "C",
            Self::N => "N",
            Self::O => "O",
        };
    }

    fn atomic_number(self) -> u32 {
        return match self {
            Self::H => 1,
            Self::C => 6,
            Self::N => 7,
            Self::O => 8,
        };
    }

    /// Mass in atomic mass units, copied from the force field file
    fn mass(self) -> f64 {
        return match self {
            Self::H => 1.008,
            Self::C => 12.,
            Self::N => 14.,
            Self::O => 15.999,
        };
    }

    /// Jmol atom colors
    fn color(self) -> Rgb<u8> {
        let (r, g, b) = match self {
            Self::H => (0xFF, 0xFF, 0xFF),
            Self::C => (0x90, 0x90, 0x90),
            Self::N => (0x30, 0x50, 0xF8),
            Self::O => (0xFF, 0x0D, 0x0D),
        };
        return Rgb::new(r, g, b);
    }

    /// Covalent atom radii in angstroms from Wikipedia table, times 1.5 because it looks good
    /// https://en.wikipedia.org/wiki/Covalent_radius
    fn radius(self) -> f32 {
        let radius = match self {
            Self::H => 0.31,
            Self::C => 0.76,
            Self::N => 0.71,
            Self::O => 0.66,
        };
        return radius * 1.5;
    }

    /// Keyboard key for selecting this atom type
    fn hotkey(self) -> Key {
        return match self {
            Self::H => Key::Key1,
            Self::C => Key::Key2,
            Self::N => Key::Key3,
            Self::O => Key::Key4,
        };
    }
}

/// In femtoseconds
const TIMESTEP: f64 = 0.2;
const STEPS_PER_FRAME: u32 = 50;
/// In angstroms
const INITIAL_BOX_SIZE: f64 = 20.;
/// In femtoseconds
const THERMOSTAT_TIME_SCALE: f64 = 200.;
/// In femtoseconds
const BAROSTAT_TIME_SCALE: f64 = 1000.;

/// In kelvins
const THERMOSTAT_MAX_TEMP: f64 = 10000.;
/// In kelvins
const THERMOSTAT_MIN_NONZERO_TEMP: f64 = 0.1;
/// In atmospheres
const BAROSTAT_MAX_PRESS: f64 = 10000.;
/// In atmospheres
const BAROSTAT_MIN_NONZERO_PRESS: f64 = 0.1;

/// Using CHON-2019 ReaxFF force field from Kowalik et al. https://doi.org/10.1021/acs.jpcb.9b04298
const FORCE_FIELD_FILE: &str = include_str!("../ffield.reax.chon2019");

/// Create a temporary file with the force field file contents in it and return a path to it. The
/// file will be deleted when this path is dropped.
fn force_field_path() -> TempPath {
    let mut temp_file = NamedTempFile::new().expect("temporary file should be created");
    temp_file
        .write_all(FORCE_FIELD_FILE.as_bytes())
        .expect("temporary file should be written");
    return temp_file.into_temp_path();
}

/// Deletes a single atom which is closest to given position. Lammps doesn't have any native way to
/// delete just one atom, only all atoms in a region, so we have to implement it ourselves using
/// atom IDs. This does not search for closest atoms over a periodic boundary. Unsafe when given
/// atom count is too high.
unsafe fn delete_closest_atom(simulation: &mut Lammps, atom_count: usize, pos: DVec3) {
    let ids = unsafe { simulation.extract_atom_id(atom_count) };
    let positions = unsafe { simulation.extract_atom_x(atom_count) };
    let (closest_id, _) = ids
        .zip(positions.map(|atom_pos| DVec3::from(atom_pos).distance_squared(pos)))
        .min_by(|(_, dist_a), (_, dist_b)| dist_a.partial_cmp(dist_b).unwrap())
        .expect("there should be atoms when deleting one");
    simulation.command(&format!("group temp id {}", closest_id));
    simulation.command("delete_atoms group temp");
    simulation.command("group temp delete");
}

/// The simulation box uses periodic boundary, which means it lives in a kind of infinite periodic
/// space where the same box repeats over and over again. This function maps simulation coordinates
/// to a shifted view of the box that is centered around a given center. The coordinates are allowed
/// be outside the simulation box.
fn sim_to_visual_pos(visual_center: Vec3, box_low: DVec3, box_high: DVec3, sim_pos: DVec3) -> Vec3 {
    let box_center = box_low.midpoint(box_high);
    let visual_offset = visual_center.as_dvec3() - box_center;
    let box_size = box_high - box_low;

    let box_relative_sim_pos = sim_pos - box_low;
    // Vector distance from the sim_pos to the box_low (the lowest corner) of the visual box view
    let visual_box_low_distance = visual_offset - box_relative_sim_pos;
    // Ceiling the distance to the next whole number of box lengths. The distance vector will now
    // point inside the visual box.
    let visual_pos_distance = (visual_box_low_distance / box_size).ceil() * box_size;
    let visual_pos = sim_pos + visual_pos_distance;
    return visual_pos.as_vec3();
}

/// Maps coordinates from a shifted box view back to the corresponding real simulation coordinates.
fn visual_to_sim_pos(box_low: DVec3, box_high: DVec3, visual_pos: Vec3) -> DVec3 {
    let box_size = box_high - box_low;

    let box_relative_visual_pos = visual_pos.as_dvec3() - box_low;
    // This has to be Euclidean or flooring modulo, not truncating like %, to work with negative
    // coordinates
    let box_relative_sim_pos = box_relative_visual_pos.rem_euclid(box_size);
    let sim_pos = box_relative_sim_pos + box_low;
    return sim_pos;
}

/// Updates a visual box view center position after the box size has changed.
fn update_visual_center(
    old_box_low: DVec3,
    old_box_high: DVec3,
    new_box_low: DVec3,
    new_box_high: DVec3,
    old_visual_center: Vec3,
) -> Vec3 {
    let old_box_size = old_box_high - old_box_low;
    let new_box_size = new_box_high - new_box_low;
    let box_size_change_factor = new_box_size / old_box_size;
    let box_center = new_box_low.midpoint(new_box_high);

    let old_visual_offset = old_visual_center.as_dvec3() - box_center;
    let new_visual_offset = old_visual_offset * box_size_change_factor;
    let new_visual_center = new_visual_offset + box_center;
    return new_visual_center.as_vec3();
}

#[kiss3d::main]
async fn main() {
    // Initialize window
    let mut window = Window::new_with_size("Atoms", 1000, 800).await;

    let mut camera = OrbitCamera3d::new(Vec3::new(0., 0., -30.), Vec3::ZERO);
    camera.set_dist_step(0.99);

    let mut scene = SceneNode3d::empty();
    let mut light = scene.add_light(Light::default());

    // Initialize simulation
    let mut simulation = Lammps::open(&["-log", "none"]);
    simulation.command("units real");
    simulation.command("dimension 3");
    simulation.command("boundary p p p");
    simulation.command("atom_style charge");
    simulation.command(&format!("timestep {}", TIMESTEP));

    // Initialize simulation box
    simulation.command(&format!(
        "region box block {} {} {} {} {} {}",
        -INITIAL_BOX_SIZE / 2.,
        INITIAL_BOX_SIZE / 2.,
        -INITIAL_BOX_SIZE / 2.,
        INITIAL_BOX_SIZE / 2.,
        -INITIAL_BOX_SIZE / 2.,
        INITIAL_BOX_SIZE / 2.,
    ));
    simulation.command(&format!("create_box {} box", AtomType::iter().len()));

    let mut previous_box_low = DVec3::splat(-INITIAL_BOX_SIZE / 2.);
    let mut previous_box_high = DVec3::splat(INITIAL_BOX_SIZE / 2.);

    // Initialize simulation atom types and force field
    for atom_type in AtomType::iter() {
        simulation.command(&format!("mass {} {}", atom_type as i32, atom_type.mass()));
    }
    // Hack: arbitrary values to prevent crashing, kokkos should be better but it crashes even
    // faster
    simulation.command("pair_style reaxff NULL safezone 5 mincap 300 minhbonds 300");
    {
        // Hack: lammps only takes reaxff force field parameters as a file path, so making a
        // temporary file
        let temp_path = force_field_path();

        simulation.command(&format!(
            "pair_coeff * * {} {}",
            temp_path.display(),
            AtomType::iter()
                .map(|atom_type| atom_type.symbol())
                .join(" "),
        ));
    }
    simulation.command("fix 2 all qeq/reaxff 1 0 10 1e-6 reaxff");

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
    let mut template_atom_type = AtomType::H;
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
                            let sim_template_pos = visual_to_sim_pos(
                                previous_box_low,
                                previous_box_high,
                                template_sphere.position(),
                            );
                            simulation.create_atom(
                                template_atom_type as i32,
                                sim_template_pos.to_array(),
                                [0., 0., 0.],
                            );
                            atom_spheres.push(scene.add_sphere(0.));

                            template_add_mode = false;
                        } else {
                            if let Some(selected_pos) = selection {
                                // Delete atom
                                let sim_selected_pos = visual_to_sim_pos(
                                    previous_box_low,
                                    previous_box_high,
                                    selected_pos,
                                );
                                // Safety:
                                // We always push created atoms and pop deleted atoms from
                                // atom_spheres, and lammps never changes the atom count by itself,
                                // so the length of atom_spheres is also the number of atoms in the
                                // simulation.
                                unsafe {
                                    delete_closest_atom(
                                        &mut simulation,
                                        atom_spheres.len(),
                                        sim_selected_pos,
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

                WindowEvent::Key(key, Action::Press, _) if !window.is_egui_capturing_keyboard() => {
                    for atom_type in AtomType::iter() {
                        if atom_type.hotkey() == key {
                            template_atom_type = atom_type;
                            break;
                        }
                    }
                }

                // Update button drag monitoring. Some actions should only happen when buttons are
                // clicked without dragging so we keep track of it. Button is considered dragged
                // when it moves more than 10 pixels
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
            Panel::top("top_panel").show(ctx, |ui| {
                ui.horizontal(|ui| {
                    for atom_type in AtomType::iter() {
                        // Make the buttons look like periodic table
                        let Rgb { r, g, b } = atom_type.color();
                        let background_color = Color32::from_rgb(r, g, b);

                        let mut button = Button::new(())
                            .selected(atom_type == template_atom_type)
                            .fill(background_color)
                            .min_size(egui::vec2(60., 60.));
                        if atom_type == template_atom_type {
                            button = button.stroke((3., Color32::LIGHT_GRAY));
                        }

                        let button = ui.add(button);

                        let text_color = Color32::from_rgb(0x10, 0x10, 0x10);
                        let symbol_rect = ui.painter().text(
                            button.rect.center(),
                            Align2::CENTER_CENTER,
                            atom_type.symbol(),
                            FontId::proportional(36.),
                            text_color,
                        );
                        ui.painter().text(
                            symbol_rect.center_top(),
                            Align2::CENTER_CENTER,
                            atom_type.atomic_number(),
                            FontId::default(),
                            text_color,
                        );
                        ui.painter().text(
                            symbol_rect.center_bottom(),
                            Align2::CENTER_CENTER,
                            atom_type.mass(),
                            FontId::default(),
                            text_color,
                        );

                        if button.clicked() {
                            template_atom_type = atom_type;
                        }
                    }
                });
            });
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

                        ui.spacing_mut().slider_width = bar_width;
                        let slider = ui.add(
                            Slider::new(&mut thermostat_temperature, 0.0..=THERMOSTAT_MAX_TEMP)
                                .logarithmic(true)
                                .smallest_positive(THERMOSTAT_MIN_NONZERO_TEMP)
                                .suffix(" K"),
                        );
                        if slider.changed() {
                            ensemble_changed = true;
                        }

                        ui.add(
                            ProgressBar::new(
                                ((temperature.log2() - THERMOSTAT_MIN_NONZERO_TEMP.log2())
                                    / (THERMOSTAT_MAX_TEMP.log2()
                                        - THERMOSTAT_MIN_NONZERO_TEMP.log2()))
                                    as f32,
                            )
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
                        ui.spacing_mut().slider_width = bar_width;
                        let slider = ui.add(
                            Slider::new(&mut barostat_pressure, 0.0..=BAROSTAT_MAX_PRESS)
                                .logarithmic(true)
                                .smallest_positive(BAROSTAT_MIN_NONZERO_PRESS)
                                .suffix(" atm"),
                        );
                        if slider.changed() {
                            ensemble_changed = true;
                        }

                        ui.add(
                            ProgressBar::new(
                                ((pressure.log2() - BAROSTAT_MIN_NONZERO_PRESS.log2())
                                    / (BAROSTAT_MAX_PRESS.log2()
                                        - BAROSTAT_MIN_NONZERO_PRESS.log2()))
                                    as f32,
                            )
                            .desired_width(bar_width)
                            .text(format!("{:.2} atm", pressure)),
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
                    "fix 1 all npt temp {} {} {} iso {} {} {}",
                    thermostat_temperature,
                    thermostat_temperature,
                    THERMOSTAT_TIME_SCALE,
                    barostat_pressure,
                    barostat_pressure,
                    BAROSTAT_TIME_SCALE,
                ));
            } else if thermostat_enabled {
                simulation.command(&format!(
                    "fix 1 all nvt temp {} {} {}",
                    thermostat_temperature, thermostat_temperature, THERMOSTAT_TIME_SCALE
                ));
            } else if barostat_enabled {
                simulation.command(&format!(
                    "fix 1 all nph iso {} {} {}",
                    barostat_pressure, barostat_pressure, BAROSTAT_TIME_SCALE
                ));
            } else {
                simulation.command("fix 1 all nve");
            }

            ensemble_changed = false;
        }

        // Run simulation one frame forward
        if simulation_running {
            simulation.command(&format!("run {}", STEPS_PER_FRAME));
        }

        // Update box position in 3D
        let (box_low, box_high) = simulation.extract_box();
        let (box_low, box_high) = (DVec3::from(box_low), DVec3::from(box_high));

        let scale = (box_high - box_low).as_vec3();
        box_cuboid.set_local_scale(scale.x, scale.y, scale.z);

        // Readjust camera focus point after box size has possibly changed
        let new_camera_focus = update_visual_center(
            previous_box_low,
            previous_box_high,
            box_low,
            box_high,
            camera.at(),
        );
        camera.set_at(new_camera_focus);

        // Box visual position is always centered on the camera focus
        box_cuboid.set_position(camera.at());

        (previous_box_low, previous_box_high) = (box_low, box_high);

        // Update atom spheres in 3D
        // Safety:
        // We always push created atoms and pop deleted atoms from atom_spheres, and lammps never
        // changes the atom count by itself, so the length of atom_spheres is also the number of
        // atoms in the simulation.
        let positions = unsafe { simulation.extract_atom_x(atom_spheres.len()) };
        let types = unsafe { simulation.extract_atom_type(atom_spheres.len()) };

        for (atom_sphere, (atom_pos, atom_type)) in
            atom_spheres.iter_mut().zip(positions.zip(types))
        {
            let atom_type = AtomType::from_repr(atom_type as usize)
                .expect("lammps should not return unknown atom types");
            atom_sphere.set_color(atom_type.color().map(|c| c as f32 / 256.).with_alpha(1.));
            let scale = atom_type.radius() * 2.;
            atom_sphere.set_local_scale(scale, scale, scale);

            let visual_atom_pos =
                sim_to_visual_pos(camera.at(), box_low, box_high, DVec3::from(atom_pos));
            atom_sphere.set_position(visual_atom_pos);
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
                let scale = template_atom_type.radius() * 2.;
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
