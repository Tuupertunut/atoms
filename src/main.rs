use kiss3d::{
    camera::{ArcBall, Camera},
    egui::{Grid, ProgressBar, Slider, TopBottomPanel},
    event::{Action, Key, MouseButton, WindowEvent},
    light::Light,
    nalgebra::{self, Point2, Point3, Translation3},
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
    simulation.command("mass 1 40");
    simulation.command("pair_style lj/cut 7");
    simulation.command("pair_coeff 1 1 0.01 3.3");
    simulation.command("fix 1 all nve");
    simulation.command("timestep 0.001");

    let mut n_atoms = simulation.get_natoms();

    let mut window = Window::new_with_size("Atoms", 1000, 800);
    window.set_light(Light::StickToCamera);

    let mut camera = ArcBall::new(Point3::new(10., 10., -20.), Point3::new(10., 10., 10.));
    camera.set_dist_step(0.99);

    let mut box_cuboid = window.add_cube(0., 0., 0.);
    box_cuboid.set_surface_rendering_activation(false);
    box_cuboid.set_lines_width(1.);

    let mut atom_spheres = iter::repeat_with(|| window.add_sphere(1.))
        .take(n_atoms)
        .collect::<Vec<_>>();

    let mut template_sphere = window.add_sphere(1.);
    template_sphere.set_surface_rendering_activation(false);
    template_sphere.set_lines_width(0.5);
    template_sphere.set_visible(false);

    let mut template_distance = 30.;

    let mut last_button1_pressed_pos = Point2::origin();
    let mut button1_pressed_without_dragging = false;
    let mut last_button2_pressed_pos = Point2::origin();
    let mut button2_pressed_without_dragging = false;

    let mut simulation_running = true;

    let mut thermostat_enabled = false;
    let mut thermostat_temperature = 300.;
    let mut barostat_enabled = false;
    let mut barostat_pressure = 300.;

    while window.render_with_camera(&mut camera).await {
        if simulation_running {
            simulation.command("run 50");
        }

        for mut event in window.events().iter() {
            match event.value {
                WindowEvent::MouseButton(MouseButton::Button1, Action::Press, _)
                    if !window.is_egui_capturing_mouse() =>
                {
                    let cursor_pos = window.cursor_pos().unwrap();
                    last_button1_pressed_pos = Point2::new(cursor_pos.0, cursor_pos.1);
                    button1_pressed_without_dragging = true;
                }
                WindowEvent::MouseButton(MouseButton::Button1, Action::Release, _)
                    if !window.is_egui_capturing_mouse() =>
                {
                    if template_sphere.is_visible() && button1_pressed_without_dragging {
                        let template_pos = template_sphere
                            .data()
                            .local_translation()
                            .cast::<f64>()
                            .vector
                            .data
                            .0[0];
                        simulation.create_atom(1, template_pos, [0., 0., 0.]);
                        atom_spheres.push(window.add_sphere(1.));
                        n_atoms += 1;

                        template_sphere.set_visible(false);
                    }
                    button1_pressed_without_dragging = false;
                }
                WindowEvent::MouseButton(MouseButton::Button2, Action::Press, _)
                    if !window.is_egui_capturing_mouse() =>
                {
                    let cursor_pos = window.cursor_pos().unwrap();
                    last_button2_pressed_pos = Point2::new(cursor_pos.0, cursor_pos.1);
                    button2_pressed_without_dragging = true;
                }
                WindowEvent::MouseButton(MouseButton::Button2, Action::Release, _)
                    if !window.is_egui_capturing_mouse() =>
                {
                    if button2_pressed_without_dragging {
                        template_sphere.set_visible(!template_sphere.is_visible());
                    }
                    button2_pressed_without_dragging = false;
                }
                WindowEvent::CursorPos(cursor_x, cursor_y, _) => {
                    let moved_pos = Point2::new(cursor_x, cursor_y);
                    if button1_pressed_without_dragging
                        && nalgebra::distance_squared(&last_button1_pressed_pos, &moved_pos)
                            >= f64::powi(10., 2)
                    {
                        button1_pressed_without_dragging = false;
                    }
                    if button2_pressed_without_dragging
                        && nalgebra::distance_squared(&last_button2_pressed_pos, &moved_pos)
                            >= f64::powi(10., 2)
                    {
                        button2_pressed_without_dragging = false;
                    }
                }
                WindowEvent::Scroll(_, y_offset, _) if !window.is_egui_capturing_mouse() => {
                    if template_sphere.is_visible() {
                        event.inhibited = true;
                        template_distance *= f32::powf(1.01, y_offset as f32);
                    }
                }
                WindowEvent::Key(Key::Space, Action::Press, _)
                    if !window.is_egui_capturing_keyboard() =>
                {
                    simulation_running = !simulation_running;
                }
                _ => {}
            }
        }

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

        if let Some(cursor_pos) = window.cursor_pos() {
            let (origin, direction) = camera.unproject(
                &Point2::new(cursor_pos.0, cursor_pos.1).cast::<f32>(),
                &window.size().cast::<f32>(),
            );
            template_sphere
                .set_local_translation(Translation3::from(origin + direction * template_distance));
        }

        window.draw_ui(|ctx| {
            TopBottomPanel::bottom("bottom_panel").show(ctx, |ui| {
                Grid::new("stat_grid").num_columns(2).show(ui, |ui| {
                    let mut bar_width = 0.;

                    ui.checkbox(&mut thermostat_enabled, "Thermostat");
                    ui.vertical(|ui| {
                        // Calculating bar width here because we don't know the available space in
                        // the second column before this point.
                        bar_width = f32::max(0., ui.available_width() - 80.);

                        ui.spacing_mut().slider_width = bar_width;
                        ui.add(Slider::new(&mut thermostat_temperature, 0.0..=1200.).suffix(" K"));
                        ui.add(
                            ProgressBar::new(240. / 1200.)
                                .desired_width(bar_width)
                                .text(format!("{:.2} K", 240.)),
                        );
                    });
                    ui.end_row();

                    ui.checkbox(&mut barostat_enabled, "Barostat");
                    ui.vertical(|ui| {
                        ui.spacing_mut().slider_width = bar_width;
                        ui.add(Slider::new(&mut barostat_pressure, 0.0..=1200.).suffix(" bar"));
                        ui.add(
                            ProgressBar::new(0. / 1200.)
                                .desired_width(bar_width)
                                .text(format!("{:.2} bar", 0.)),
                        );
                    });
                    ui.end_row();
                });
            });
        });
    }
}
