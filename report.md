This report was written at the time of Atoms 1.0.0.

### Introduction

I wanted to create an application that would show people how atoms really behave, so they could gain an intuition on the laws of molecular dynamics and try things out for themselves. Similar educational simulations have previously been very helpful for me in teaching intuitive understanding of physics, such as Falstad's [gas](https://www.falstad.com/gas/), [wave](https://www.falstad.com/ripple/) and [circuit](https://www.falstad.com/circuit/) simulators, or on a larger scale [Kerbal Space Program](https://www.kerbalspaceprogram.com/) for orbital mechanics and [SpaceEngine](https://spaceengine.org/) for the structure of the universe. I have tried to find a similar application for molecular dynamics for a long time but found none. Only professional simulators like [LAMMPS](https://www.lammps.org/) or [very simple school examples](https://mw.concord.org/nextgen/#interactives/chemistry/chemical-reactions/reaction-hydrogen-oxygen-molecules) seemed to be available. It's time molecular dynamics got its own educational simulator.

### Plan

The plan was to create a game-like application with 3D user interface showing a running [LAMMPS](https://www.lammps.org/) simulation in real time.

1. Technology stack

    - I chose [Rust](https://rust-lang.org/) as the programming language for this project, both because I was already familiar with it as well as because I know it has enough performance to run a heavy simulation.

    - The frontend 3D UI would be implemented using a Rust 3D graphics library [Kiss3d](https://github.com/sebcrozet/kiss3d) which I was already familiar with.

    - The backend simulation engine would be LAMMPS and it would be called through its [C API](https://docs.lammps.org/Library.html#lammps-c-library-api).

2. The application needs a 3D user interface that shows the box, the atoms and allows moving a camera around and in the box and interacting with the simulation.

    - The simulation box would be a cube outline and the atoms would be solid spheres in the 3D space, in a similar way to how [OVITO](https://www.ovito.org/) displays them.

    - The camera could be either a first person camera that the user needs to fly around the box, or a rotating camera that rotates around the center of the box and allows zooming in and out. Kiss3d provides both of these out of the box.

    - In addition to displaying the simulation, the user needs to be able to add new atoms at a selected position and delete old ones through the UI.

    - Temperature and pressure and possibly other measurements should be displayed on the UI, as well as controls for adjusting the thermostat temperature and barostat pressure.

    - There needs to be some way to pause the simulation, so that the user has more time to observe the current configuration of atoms.

    - Other controls may be implemented for selecting box size, time speed, atom type, etc.

3. The simulation is run in LAMMPS in the background.

    - First a noble gas like argon could be used with a Lennard-Jones potential.

    - Later it would support other atom types as well such as hydrogen and oxygen. Chemical bond formation and breaking could be implemented with reactive potentials or ab initio MD.

    - Simulation box boundary would be periodic. A hard boundary where atoms bounce off could also be a good choice, but it's not available in LAMMPS by default.

    - The thermodynamic ensemble would first be NVE, but it would need to change mid-run as thermostat and barostat are turned on.
    
    - The simulation would need to run in short bursts and return new atom positions every time the frontend asks for an new frame. The frame rate of the front end would probably be 60 Hz, meaning that the simulation would need to stop, send positions and start again 60 times per second.

    - The simulation needs to support adding and deleting atoms between each frame.

### Implementation

- First I needed to find out how to call LAMMPS from Rust.

    - As a native compiled language, Rust supports linking C libraries and calling C APIs through Rust's Foreign Function Interface (FFI), but it requires a [bindgen](https://rust-lang.github.io/rust-bindgen/) tool to generate Rust bindings from a C header. I added a bindgen phase to the build script of the Rust application to generate bindings for LAMMPS.

    - To statically link the liblammps library to the application I also needed to compile LAMMPS from source. I concluded that for my purpose of running a simple molecular dynamics simulation I needed zero LAMMPS packages. Thus I was able to compile a really light LAMMPS with all packages and extra features disabled and it had almost zero dependencies. The extra features I disabled were `MPI`, `JPEG`, `PNG`, `FFMPEG`, `GZIP` and `CURL`. I left OpenMP enabled in case LAMMPS would be able to automatically parallelize some workload but I haven't seen any parallelization. I added a LAMMPS compile phase in the build script of the Rust application so that it automatically compiles a custom LAMMPS for linking.

    - The LAMMPS C API was surprisingly narrow. Most functionalities were only available through a generic function sending LAMMPS scripting commands as strings. Only session opening/closing, atom creation and data reading were available as specialized functions. Simulation needed to be initialized with a list of LAMMPS commands just like when running LAMMPS as a script. Luckily there were specialized functions for atom data reading as those would become useful when reading new atom positions every frame. I created a small Rust wrapper (`lammps.rs`) for the LAMMPS API to allow using it in a more Rust-friendly way.

- The application sets up the LAMMPS environment at startup.
    - Logging is disabled to avoid generating LAMMPS log files in every directory where the application is run. Simply running `log none` command is not enough to disable a log file, so it is given as a command line parameter to `Lammps::open`.

    - `atom_style` is `atomic` because this application deals with only single uncharged spherical atoms. `units` are `metal` simply because I was already familiar with them.

    - A 3-dimensional periodic simulation box is created with a default side length 20. The side length may change when the barostat is enabled.

    - [Argon](https://en.wikipedia.org/wiki/Argon) was chosen as the element of the atoms with a mass of 40 u. A Lennard-Jones potential is used with parameters `epsilon = 0.01 eV`, `sigma = 3.4 Å` ([source](https://www.researchgate.net/figure/Lennard-Jones-LJ-potential-parameters-of-different-materials-considered-in-thepresent_tbl2_319412425)) and `cutoff = 15 Å`.

    - A timestep of 1 fs was selected because it seems small enough.
  
- A Kiss3d scene displays box and atoms in 3D.

    - The simulation box is drawn as a cube without faces, so just an outline.

    - The argon atoms are solid spheres with an arbitrary radius of 1 Å which looks good enough. Their color is `#80D1E3` which is the [Jmol](https://jmol.sourceforge.net/jscolors/) standard color for argon.

    - On each frame LAMMPS is asked to run the simulation 50 timesteps forward with command `run 50`. The value 50 was experimentally tested to look good in real time. After running, box and atom coordinates are read from LAMMPS and updated to the box cube and atom spheres in the 3D scene.

    - A Kiss3d ArcBall camera is used, which allows rotating the camera around a center point. The center point is by default at the center of the simulation box. The camera can be rotated by dragging the **left mouse button** and the center point can be moved by dragging the **right mouse button**. Scrolling the **mouse wheel** zooms in and out of the center point.

    - Kiss3d also supports drawing a regular 2D GUI in front of the 3D scene with [egui](https://www.egui.rs/), which is a popular GUI library for Rust. A temperature and pressure meter and thermostat and barostat controls were drawn with egui at the bottom of the screen.

        - Unlike the 3D objects, the egui GUI has no persistent state and is redrawn every frame. This is a feature of egui as it is an [immediate mode](https://en.wikipedia.org/wiki/Immediate_mode_(computer_graphics)) GUI library.

- The UI allows interaction with the simulation, such as deleting and adding atoms at any location with mouse.

    - Kiss3d allows registering listeners for different mouse and keyboard buttons and mouse moves.

    - When an atom is **hovered** with the mouse cursor, a red marker is shown around it. **Left clicking** will delete the atom.

        - The hovered atom is found by casting a ray from the camera in the direction of the cursor and checking every atom whether the ray intersects them. The closest intersected atom is picked.

        - LAMMPS does not directly allow deleting a single atom, so a hack was used. First a very small (0.01 Å) region named `temp` is defined around the atom center, next every atom inside the region is deleted and finally the region itself is deleted.

    - **Right clicking** displays a green atom template sphere that follows the mouse cursor and can be brought closer or farther by scrolling the **mouse wheel**. **Left clicking** with the template visible will add an atom to its position.

    - The simulation can be paused or continued by pressing **space** on the keyboard. Pausing allows adding or deleting atoms precisely. From the LAMMPS perspective pausing just means that no `run` commands are executed on those frames.

    - At the bottom panel, temperature and pressure can be seen in real time at the same time as the thermostat and barostat can be adjusted.

        - The barostat and thermostat controls allow changing the LAMMPS thermodynamic ensemble on the fly. Every time they are adjusted, an `unfix` LAMMPS command is executed followed by `fix nve/nvt/nph/npt` with the relevant parameters.

### Further discussion

The simulator is in a very functional and usable state, but it is by no means complete. There are many further development ideas that could be implemented in the future.

- Supporting multiple atom types and chemical bond formation is an important goal that would make the simulator a lot more useful in forming a molecular dynamics intuition. The current application structure would support this with quite small changes.

    - This could possibly make the simulation unbearably slow, especially if using quantum mechanics. Further testing is needed.

    - A periodic table could be displayed when adding new atoms and the user could select an element there. The Jmol color palette could be used to color atoms based on their type.

- Multithreading the simulation. I don't know which operations in LAMMPS support multithreading but I have heard most simulations could be multithreaded using OpenMP or MPI. This could make them faster if performance becomes an issue.

- A first person camera could be useful in certain cases, such as taking the camera inside a clump of atoms and taking a closer look at what's inside. This would require minimal changes to the application as Kiss3d already has a first person camera type.

- Allowing temperature/pressure scales to adjust dynamically to higher values. The current maximum for both is 1200 units, but it's quite easy to make the current simulation go above that.

    - The scales could maybe be logarithmic to allow new orders of magnitude at the cost of accuracy at higher values. I don't think accuracy at high values is needed anyway so I don't see any downsides.

- Adjusting box size, time speed or thermostat/barostat dampening time manually. There might be some cases where these would be useful.

- Show list of keyboard/mouse controls inside the application. This could make the application easier to use for new users.

### Conclusion

An educational simulator is by nature evaluated very subjectively, so I can only draw subjective conclusions from it. The simulator is in a good shape. It allows me to see argon in solid, liquid and gas phases, although not entirely at the correct temperatures/pressures. Real argon should melt at 84 K, but in the box it seems to melt at around 60 K. It is inaccurate but qualitatively correct. The barostat and thermostat seem to work, but their values seem to fluctuate a lot, so it's usually easier to turn the -stat off when it has reached the desired value. There is a visual bug in the bottom panel when there are no atoms in the box, but that's not really a big issue. I am personally really pleased with the end result. It feels intuitive and easy to use. The code is quite clean and I can't wait to get to develop new features for it.
