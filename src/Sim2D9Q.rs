use core::{f64, panic};
use std::{
    isize, mem,
    ops::{Index, IndexMut},
    usize,
};

use rand::{seq, Rng};
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
struct LBM2D {
    mesh: Vec<Cell2D9Q>,
    buffer: Vec<Cell2D9Q>,
    geometry: Vec<bool>,

    /// Dimensinos given as (Width, Height) => (num columns, num rows)
    dimensions: (usize, usize),

    /// Relaxation Time
    tau: f64,
    /// Lattice Speed
    c_s: f64,

    /// Total time
    time: f64,
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
    fn new(x_count: usize, y_count: usize, delta_t: f64) -> Self {
        let mesh = (0..x_count * y_count).map(|_| Cell2D9Q::new()).collect();
        let buffer = (0..x_count * y_count).map(|_| Cell2D9Q::new()).collect();
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
            buffer,
            geometry,
            dimensions: (x_count, y_count),
            // TODO replace hard coded values
            time: 0.0,
            delta_time: delta_t,
            c_s: 1.0,
            tau: 0.5,
            inlet_velocity: (1.0, 0.0),
            inlet_rho: 1.0,
        }
    }
}
impl LBM2D {
    fn step(&mut self, buffer: &mut Vec<Cell2D9Q>) {
        let (x_count, y_count) = self.dimensions;

        let new_timestep = coords!(self.dimensions)
            .zip(self.buffer.iter())
            // perform collision
            .map(|(pos, cell)| (pos, cell))
            // perform streaming
            .map(|(pos, cell)| (pos, cell))
            .map(|(_, cell)| cell)
            .collect::<Vec<_>>();
        // self.collision(buffer);
        // // performing streaming step
        // self.collision(buffer);

        // swap buffer and mesh
        mem::swap(buffer, &mut self.mesh);
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
        let rho = self.buffer[self.index(x, y)].rho;
        let f = &self.buffer[self.index(x, y)].f;
        let (x_limit, y_limit) = (self.dimensions.0 - 1, self.dimensions.1 - 1);

        // for each direction:
        if x == 0 {
            // Inlet (left wall)
            let (u_x, u_y) = self.inlet_velocity;

            let f_1 = f[3] + 2.0 / 3.0 * rho * u_x;
            let f_5 =
                f[7] + 1.0 / 2.0 * (f[2] - f[4]) - 1.0 / 2.0 * rho * u_y + 1.0 / 6.0 * rho * u_x;
            let f_8 =
                f[6] - 1.0 / 2.0 * (f[4] - f[2]) + 1.0 / 2.0 * rho * u_y + 1.0 / 6.0 * rho * u_x;

            [0.0, f_1, 0.0, 0.0, 0.0, f_5, 0.0, 0.0, f_8]
        } else if x == x_limit {
            // 0 Pressure case
            // Outlet (right wall)
            let rho_in = self.inlet_rho;
            let u_y = 0;

            // Solve for outlet velocity
            let u_x = 1.0 - (f[0] + f[2] + f[4] + f[1] + f[5] + f[8]) / rho_in;

            todo!()
        } else if y == 0 {
            // Top Wall
            // Only 2, 5, 6 point back into the fluid
            let (u_x, u_y) = (0.0, 0.0);

            let f_4 = f[2] + 2.0 / 3.0 * rho / u_y;
            let f_7 =
                f[5] + 1.0 / 2.0 * (f[1] - f[3]) - 1.0 / 2.0 * rho * u_x + 1.0 / 6.0 * rho * u_y;
            let f_8 =
                f[6] - 1.0 / 2.0 * (f[1] - f[3]) + 1.0 / 2.0 * rho * u_x + 1.0 / 6.0 * rho * u_y;

            [0.0, 0.0, 0.0, 0.0, f_4, 0.0, 0.0, f_7, f_8]
        } else if y == y_limit {
            // Bottom Wall
            // Only 2, 5, 6 point back into the fluid
            let (u_x, u_y) = (0.0, 0.0);

            let f_2 = f[4] + 2.0 / 3.0 * rho / u_y;
            let f_5 =
                f[7] - 1.0 / 2.0 * (f[1] - f[3]) + 1.0 / 2.0 * rho * u_x + 1.0 / 6.0 * rho * u_y;
            let f_6 =
                f[8] + 1.0 / 2.0 * (f[1] - f[3]) - 1.0 / 2.0 * rho * u_x + 1.0 / 6.0 * rho * u_y;

            [0.0, 0.0, f_2, 0.0, 0.0, f_5, f_6, 0.0, 0.0]
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

    fn fixed_velocity_boundary(i: usize, f: &[f64; 9], rho: f64, u: (f64, f64)) -> [f64; 9]{
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
                    f[7] + 1.0 / 2.0 * (f[2] - f[4]) - 1.0 / 2.0 * rho * u_y
                        + 1.0 / 6.0 * rho * u_x,
                    0.0,
                    0.0,
                    f[6] - 1.0 / 2.0 * (f[4] - f[2])
                        + 1.0 / 2.0 * rho * u_y
                        + 1.0 / 6.0 * rho * u_x,
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
                    f[7] - 1.0 / 2.0 * (f[1] - f[3])
                        + 1.0 / 2.0 * rho * u_x
                        + 1.0 / 6.0 * rho * u_y,
                    f[6] + 1.0 / 2.0 * (f[1] - f[3]) - 1.0 / 2.0 * rho * u_x
                        + 1.0 / 6.0 * rho * u_y,
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
                    f[7] + 1.0 / 2.0 * (f[2] - f[4])
                        - 1.0 / 2.0 * rho * u_y
                        - 1.0 / 6.0 * rho * u_x,
                    0.0,
                    0.0,
                    f[6] - 1.0 / 2.0 * (f[4] - f[2]) + 1.0 / 2.0 * rho * u_y
                        - 1.0 / 6.0 * rho * u_x,
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
                    f[7] - 1.0 / 2.0 * (f[1] - f[3])
                        + 1.0 / 2.0 * rho * u_x
                        + 1.0 / 6.0 * rho * u_y,
                    f[8] + 1.0 / 2.0 * (f[1] - f[3]) - 1.0 / 2.0 * rho * u_x
                        + 1.0 / 6.0 * rho * u_y,
                    0.0,
                    0.0,
                ]
            }
            _ => panic!("Invalid used of boundary condition function"), // Default: No bounce-back for other directions
        }
    }
}

// Apply collision step to the mesh
//
// Takes the previous mesh, returning the new mesh.
//
// => f_i(x, t + delta_t) = f_i(x, t) + f_i^eq(x, t) - f_i(x, t) / tau_f
// => f_i^eq(x, t) = omega_i * rho * (1 + 3 u_i)
// where
// - func f is the density distribution
// - i \in [0, 8] and is the index for the D2Q9 Lattice
// - tau_f is characteristic timescale
// - omega_i is the lattice weight
// - e_i has same dimension as u (i.e. 2d)
//
// Note: e_i * u => u_i in this implementation.
// fn collision(prev: &Mesh) -> Mesh {
//     // refer back to wikipedia article (actually quite good)
//     // <https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods>
//     let mut new = prev.clone();
//
//     for (x, y) in coords!(prev.dimensions) {
//         let e = &DIRS;
//         let u = &prev[x][y].velocity;
//         let w = &prev[x][y].density;
//         let rho = &prev[x][y].rho;
//
//         let collect: [f64; L_COUNT] = (0..L_COUNT)
//             .map(|i| {
//                 let e_dot_u = e[i].0 as f64 * u.0 + e[i].1 as f64 * u.1;
//                 let u_dot_u = u.0 * u.0 + u.1 * u.1;
//
//                 w[i] * rho
//                     * (1.0
//                         + 3.0 / C.powi(2) * e_dot_u
//                         + 9.0 / (2.0 * C.powi(4)) * e_dot_u.powi(2)
//                         + 3.0 / (2.0 * C.powi(2)) * u_dot_u.powi(2))
//             })
//             .collect::<Vec<f64>>()
//             .try_into()
//             .unwrap();
//     }
//     todo!()
// }

// Applies the streaming step
//
// Takes the new mesh (after appling the collision step) and modifies it to stream.
// fn streaming(prev: &Mesh, new: &mut Mesh, geometry: &Vec<bool>) {
//     for (x, y) in coords!(new.dimensions) {
//         for (l) in (0..L_COUNT) {
//             let (dx, dy) = DIRS[l];
//             let nx = x as isize + dx;
//             let ny = y as isize + dy;
//
//             let (x_limit, y_limit) = (new.dimensions.0 as isize, new.dimensions.1 as isize);
//             // TODO handle boundary conditions
//             if (nx >= 0 && nx < x_limit) && (ny >= 0 && ny < y_limit) {
//                 // Check geometry
//                 if geometry[y * x_limit as usize + x] {
//                     // flip direction
//                     let bounce_l = (l + 4) % 8 + 1;
//                     new[x][y].density[l] = prev[x][y].density[bounce_l];
//                 }
//
//                 if (nx == 0) {}
//
//                 if (nx == x_limit) {}
//             }
//         }
//     }
// }
