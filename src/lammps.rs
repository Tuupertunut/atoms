use lammps_sys::*;
use std::{
    ffi::{CString, c_void},
    ptr,
};

/// Lammps C API bindings generated from library.h
pub mod lammps_sys {
    #![allow(nonstandard_style, dead_code)]
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

/// A semi-safe wrapper for lammps C API
pub struct Lammps {
    pub session: *mut c_void,
}

impl Lammps {
    pub fn open() -> Self {
        return Self {
            session: unsafe { lammps_open_no_mpi(0, ptr::null_mut(), ptr::null_mut()) },
        };
    }

    pub fn command(&mut self, command: &str) {
        unsafe { lammps_command(self.session, CString::new(command).unwrap().as_ptr()) };
    }

    pub fn get_natoms(&self) -> usize {
        return unsafe { lammps_get_natoms(self.session) } as usize;
    }

    pub fn extract_atom(&self, property: &str) -> *const c_void {
        return unsafe {
            lammps_extract_atom(self.session, CString::new(property).unwrap().as_ptr())
        };
    }

    pub fn extract_box(&self) -> ([f64; 3], [f64; 3]) {
        let mut box_bounds = ([0.; 3], [0.; 3]);
        unsafe {
            lammps_extract_box(
                self.session,
                &mut box_bounds.0 as *mut f64,
                &mut box_bounds.1 as *mut f64,
                ptr::null_mut(),
                ptr::null_mut(),
                ptr::null_mut(),
                ptr::null_mut(),
                ptr::null_mut(),
            )
        };
        return box_bounds;
    }

    pub fn create_atom(&mut self, atom_type: i32, position: [f64; 3], velocity: [f64; 3]) {
        unsafe {
            lammps_create_atoms(
                self.session,
                1,
                std::ptr::null(),
                &atom_type,
                &position as *const f64,
                &velocity as *const f64,
                std::ptr::null(),
                0,
            )
        };
    }

    pub fn reset_box(&mut self, box_low: [f64; 3], box_high: [f64; 3]) {
        unsafe {
            lammps_reset_box(
                self.session,
                &box_low as *const f64 as *mut f64,
                &box_high as *const f64 as *mut f64,
                0.,
                0.,
                0.,
            )
        };
    }
}

impl Drop for Lammps {
    fn drop(&mut self) {
        unsafe { lammps_close(self.session) };
    }
}
