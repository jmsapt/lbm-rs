use core::{f64, panic};
use std::{isize, usize};

macro_rules! coords {
    ($x:expr, $y:expr) => {
        (0..$x).flat_map(|x| (0..$y).map(move |y| (x, y)))
    };
    ($z:expr) => {
        (0..$z.0).flat_map(|x| (0..$z.1).map(move |y| (x, y)))
    };
}
const X_COUNT: usize = 800;
const Y_COUNT: usize = 800;
const L_COUNT: usize = 9;
const TAU: f64 = 0.53;
const RHO_0: f64 = 100.0;
const C: f64 = 1.0;

const N: usize = 9;
const WEIGHTS: [f64; 9] = [
    4.0 / 9.0, // 4/9 for 0
    1.0 / 9.0, // 1/9 for [1, 4]
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0, // 1/36 for [5, 8]
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

const DIRS: [(isize, isize); 9] = [
    (0, 0),
    (0, 1),
    (1, 1),
    (1, 0),
    (1, -1),
    (0, -1),
    (-1, -1),
    (-1, 0),
    (-1, 1),
];

// pub trait LBM<T> {
//     fn step(&mut self, buffer: &mut Vec<T>);
//
//     fn collision(&self, buffer: &mut impl IntoIterator<Item = T>);
//
//     fn stream(&self, buffer: &mut impl IntoIterator<Item = T>);
// }

#[derive(Clone)]
struct Cell2D9Q {
    pub velocity: (f64, f64),
    // Cell density vector (func f)
    pub f: [f64; 9],
    pub rho: f64,
}
impl Cell2D9Q {
    fn new() -> Cell2D9Q {
        Self {
            velocity: (0.0, 0.0),
            f: [0.0; 9],
            rho: 0.0,
        }
    }
}

/// Lattice Boltzman Method (LBM) for 2D9Q
pub struct LBM2D {
    mesh: Vec<Cell2D9Q>,
    geometry: Vec<bool>,

    /// Dimensinos given as (Width, Height) => (num columns, num rows)
    dimensions: (usize, usize),

    /// Relaxation Time
    tau: f64,
    /// Total number of timesteps
    n_steps: usize,
    /// Time step
    delta_time: f64,
    /// Inlet Velocity
    inlet_velocity: (f64, f64),
    inlet_rho: f64,
}
impl LBM2D {
    fn index(&self, x: usize, y: usize) -> usize {
        return self.dimensions.0 * y + x;
    }

    /// Construct and intialise the simulation
    /// TODO: currently hardcoded, make adjustments to parameterise this instead
    pub fn new(
        x_count: usize,
        y_count: usize,
        cell_size: f64,
        inlet_c: f64,
        inlet_velocity: (f64, f64),
        inlet_density: f64,
    ) -> Self {
        let mesh = (0..x_count * y_count).map(|_| Cell2D9Q::new()).collect();

        // As defined here
        // <https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods#Lattice_units_conversion>
        let delta_t = cell_size / inlet_c;

        // Hardcoded circle
        let mut geometry = vec![false; X_COUNT * Y_COUNT];
        {
            let (x_0, y_0, r): (isize, isize, isize) = (
                X_COUNT as isize / 2,
                Y_COUNT as isize / 2,
                Y_COUNT as isize / 4,
            );
            coords!(x_count, y_count).for_each(|(x, y)| {
                // circle boundary
                geometry[y * x_count + x] =
                    (x as isize - x_0).pow(2) + (y as isize - y_0).pow(2) < r.pow(2);
            })
        }

        // // hardcode initial conditions
        // coords!(x_count, y_count).for_each(|(x, y)| {
        //     let mut rng = rand::thread_rng();
        //     for l in 0..L_COUNT {
        //         // add random pertubations, and set flow to the right
        //         mesh[y][x].density[l] = 1.0 + 0.01 * rng.gen::<f64>();
        //     }
        //     // TODO
        //     mesh[y][x].density[3] = 2.0;
        // });

        Self {
            mesh,
            geometry,
            dimensions: (x_count, y_count),
            // TODO replace hard coded values
            n_steps: 0,
            delta_time: delta_t,
            tau: 0.5,
            inlet_velocity,
            inlet_rho: inlet_density,
        }
    }

    pub fn step(&mut self) {
        use rayon::prelude::*;

        // 1. Find equillibrium, and collide. This will be collected into a buffer.
        // 2. Stream particles
        //
        // 3. Apply boundary conditions
        let new_timestep = (0..self.dimensions.0)
            .flat_map(|x| (0..self.dimensions.1).map(move |y| (x, y)))
            .collect::<Vec<_>>()
            .into_par_iter()
            .zip(self.mesh.par_iter())
            .map(|(pos, cell)| {
                // perform collision
                let f_eq = self.incompressible_collision(&cell, pos);
                // perform streaming
                let f_stream = self.stream(&self.mesh, pos);

                // sum result & calculate new properties
                let mut f = [0.0; 9];
                let mut rho = 0.0;
                let mut u = (0.0, 0.0);
                for i in 0..9 {
                    // F_i(x_i + v_i \Delta t, t + \Delta t) - F_i(x_i, t)
                    //      = -\frac{F_i(x_i, t) - F_i^eq(x_i, t)}{\tau}
                    let simga = -(f_stream[i] - f_eq[i]) / self.tau;
                    f[i] = f_eq[i] + f_stream[i];

                    // rho = \Sigma F_i
                    rho += f[i];

                    // rho = \Sigma F_i
                    u.0 += f[i] * DIRS[i].0 as f64;
                    u.1 += f[i] * DIRS[i].1 as f64;
                }
                u = (u.0 / rho, u.1 / rho);

                Cell2D9Q {
                    f,
                    rho,
                    velocity: u,
                }
            })
            .collect::<Vec<_>>();

        self.mesh = new_timestep;
        self.n_steps += 1;
    }

    /// Get cell density
    pub fn get_density(&self, x: usize, y: usize) -> f64 {
        let i = self.index(x, y);
        self.mesh[i].rho
    }

    /// Get cell velocity
    pub fn get_velocity(&self, x: usize, y: usize) -> (f64, f64) {
        let i = self.index(x, y);
        self.mesh[i].velocity
    }

    /// Performs incompressible collision, seting new_cell's density distributions to f^eq
    fn incompressible_collision(&self, prev_cell: &Cell2D9Q, pos: (usize, usize)) -> [f64; 9] {
        let mut f_eq = [0.0; 9];

        for (i, u) in DIRS.iter().enumerate() {
            let t_i = WEIGHTS[i];
            let rho = prev_cell.rho;
            let e = DIRS[i];
            let e_dot_u = e.0 * u.0 + e.1 * u.1;
            let u_dot_u = u.0 * u.0 + u.1 * u.1;

            f_eq[i] = nofmt::pls! {
                t_i * (
                    rho + 3.0 * (e_dot_u as f64)
                    + 9.0 / 2.0 * (e_dot_u.pow(2) as f64)
                    - 3.0 / 2.0 * (u_dot_u as f64)
                )
            };
        }

        f_eq
    }

    fn stream(&self, cells: &Vec<Cell2D9Q>, pos: (usize, usize)) -> [f64; 9] {
        let (x, y) = pos;
        let (x_limit, y_limit) = (self.dimensions.0 - 1, self.dimensions.1 - 1);

        let rho = self.mesh[self.index(x, y)].rho;
        let rho_in = self.inlet_rho;
        let f = &self.mesh[self.index(x, y)].f;
        let u_wall = (0.0, 0.0);

        // for each direction:
        // TODO 1. add bounceback for geometry
        // TODO 2. handle corner cases
        if x == 0 {
            // Inlet (left wall)
            let u = self.inlet_velocity;

            fixed_velocity_boundary(1, f, rho_in, u)
        } else if x == x_limit {
            // 0 Pressure case
            // Outlet (right wall)

            // Solve for outlet velocity
            let u_x = 1.0 - (f[0] + f[2] + f[4] + f[1] + f[5] + f[8]) / rho_in;
            let u = (u_x, 0.0);

            fixed_velocity_boundary(4, f, rho_in, u)
        } else if y == 0 {
            // Top Wall
            // Only 2, 5, 6 point back into the fluid
            // noraml is 4
            fixed_velocity_boundary(4, f, rho, u_wall)
        } else if y == y_limit {
            // Bottom Wall
            // Only 2, 5, 6 point back into the fluid
            // normal is 2
            fixed_velocity_boundary(2, f, rho, u_wall)
        } else {
            // Normal case
            nofmt::pls!([
                // Zero vector
                cells[self.index(x, y)].f[0],
                // Carindal Vectors
                cells[self.index(x - 1, y)].f[1],
                cells[self.index(x, y - 1)].f[2],
                cells[self.index(x + 1, y)].f[3],
                cells[self.index(x, y + 1)].f[4],
                // Diagonal Vectors
                cells[self.index(x - 1, y - 1)].f[5],
                cells[self.index(x + 1, y - 1)].f[6],
                cells[self.index(x + 1, y + 1)].f[7],
                cells[self.index(x - 1, y + 1)].f[8],
            ])
        }
    }
}
fn fixed_velocity_boundary(i: usize, f: &[f64; 9], rho: f64, u: (f64, f64)) -> [f64; 9] {
    let (u_x, u_y) = u;
    match i {
        1 => {
            // Left Wall: Only 1, 5, 8 point back into the fluid
            [
                0.0,
                f[3] + 2.0 / 3.0 * rho * u_x,
                0.0,
                0.0,
                0.0,
                f[7] + 1.0 / 2.0 * (f[2] - f[4]) - 1.0 / 2.0 * rho * u_y + 1.0 / 6.0 * rho * u_x,
                0.0,
                0.0,
                f[6] - 1.0 / 2.0 * (f[4] - f[2]) + 1.0 / 2.0 * rho * u_y + 1.0 / 6.0 * rho * u_x,
            ]
        }
        2 => {
            // Top Wall: Only 4, 7, 8 point back into the fluid
            [
                0.0,
                0.0,
                0.0,
                0.0,
                f[2] + 2.0 / 3.0 * rho / u_y,
                0.0,
                0.0,
                f[7] - 1.0 / 2.0 * (f[1] - f[3]) + 1.0 / 2.0 * rho * u_x + 1.0 / 6.0 * rho * u_y,
                f[6] + 1.0 / 2.0 * (f[1] - f[3]) - 1.0 / 2.0 * rho * u_x + 1.0 / 6.0 * rho * u_y,
            ]
        }
        3 => {
            // Right Wall: Only 1, 5, 8 point back into the fluid
            [
                0.0,
                f[3] - 2.0 / 3.0 * rho * u_x,
                0.0,
                0.0,
                0.0,
                f[7] + 1.0 / 2.0 * (f[2] - f[4]) - 1.0 / 2.0 * rho * u_y - 1.0 / 6.0 * rho * u_x,
                0.0,
                0.0,
                f[6] - 1.0 / 2.0 * (f[4] - f[2]) + 1.0 / 2.0 * rho * u_y - 1.0 / 6.0 * rho * u_x,
            ]
        }
        4 => {
            // Bottom Wall: Only 2, 5, 6 point back into the fluid
            [
                0.0,
                0.0,
                f[4] + 2.0 / 3.0 * rho / u_y,
                0.0,
                0.0,
                f[7] - 1.0 / 2.0 * (f[1] - f[3]) + 1.0 / 2.0 * rho * u_x + 1.0 / 6.0 * rho * u_y,
                f[8] + 1.0 / 2.0 * (f[1] - f[3]) - 1.0 / 2.0 * rho * u_x + 1.0 / 6.0 * rho * u_y,
                0.0,
                0.0,
            ]
        }
        _ => panic!("Invalid used of boundary condition function"), // Default: No bounce-back for other directions
    }
}
