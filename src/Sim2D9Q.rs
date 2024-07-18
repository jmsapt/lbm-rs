use core::{f64, panic};
use std::{isize, usize};

use rand::Rng;

macro_rules! coords {
    ($x:expr, $y:expr) => {
        (0..$x).flat_map(|x| (0..$y).map(move |y| (x, y)))
    };
    ($z:expr) => {
        (0..$z.0).flat_map(|x| (0..$z.1).map(move |y| (x, y)))
    };
}
// const X_COUNT: usize = 800;
// const Y_COUNT: usize = 800;
// const L_COUNT: usize = 9;
// const TAU: f64 = 0.53;
// const RHO_0: f64 = 100.0;
// const C: f64 = 1.0;
//
// const N: usize = 9;
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
    (1, 0),
    (0, 1), // source of all my pain :( unironically a couple of days of my life gone
    (-1, 0),
    (0, -1),
    (1, 1),
    (-1, 1),
    (-1, -1),
    (1, -1),
];

// pub trait LBM<T> {
//     fn step(&mut self, buffer: &mut Vec<T>);
//
//     fn collision(&self, buffer: &mut impl IntoIterator<Item = T>);
//
//     fn stream(&self, buffer: &mut impl IntoIterator<Item = T>);
// }

#[derive(Clone, Debug)]
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
    pub tau: f64,
    /// Total number of timesteps
    n_steps: usize,
    /// Time step
    pub delta_time: f64,
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
        let mut mesh = (0..x_count * y_count)
            .map(|_| Cell2D9Q::new())
            .collect::<Vec<_>>();

        // As defined here
        // <https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods#Lattice_units_conversion>
        let delta_t = cell_size / inlet_c;

        // Hardcoded circle
        let mut geometry = vec![false; x_count * y_count];
        // {
        //     let (x_0, y_0, r): (isize, isize, isize) = (
        //         X_COUNT as isize / 2,
        //         Y_COUNT as isize / 2,
        //         Y_COUNT as isize / 4,
        //     );
        //     coords!(x_count, y_count).for_each(|(x, y)| {
        //         // circle boundary
        //         geometry[y * x_count + x] =
        //             (x as isize - x_0).pow(2) + (y as isize - y_0).pow(2) < r.pow(2);
        //     })
        // }

        for x in 0..x_count {
            for y in 0..y_count {
                let i = y * x_count + x;
                // let f = 2.0 / 3.0 * inlet_density * inlet_velocity.0;
                mesh[i].f = incompressible_eq(inlet_density, inlet_velocity);
                update_cell_properties(&mut mesh[i]);
            }
        }

        coords!(x_count, y_count).for_each(|(x, y)| {
            let i = y * x_count + x;
            println!("Density values: {:?}", mesh[i].f);
        });

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
        // 3. Apply boundary conditions
        // 4. Calculate updated cell properties (u, rho)
        let coords = (0..self.dimensions.0)
            .flat_map(|x| (0..self.dimensions.1).map(move |y| (x, y)))
            .collect::<Vec<_>>();

        let collision = coords
            .par_iter()
            .zip(self.mesh.par_iter())
            .map(|(_pos, cell)| {
                // perform collision
                let f_eq = incompressible_eq(cell.rho, cell.velocity);
                let mut collision = Cell2D9Q::new();
                for i in 0..9 {
                    collision.f[i] = cell.f[i] - (cell.f[i] - f_eq[i]) / self.tau;
                }
                collision
            })
            .collect::<Vec<Cell2D9Q>>();

        let new_timestep = coords
            .par_iter()
            .map(|pos| {
                // Stream based on directions
                // NOTE: This is considering streaming othercells into the current cell, which is
                // the opposite of the common formula
                //
                // f(x, t + dt) = f_collision(x - e, t);
                let f_new = self.stream(&collision, *pos);
                f_new
            })
            .map(|f| {
                let mut new_cell = Cell2D9Q::new();
                new_cell.f = f;
                update_cell_properties(&mut new_cell);

                new_cell
            })
            .collect::<Vec<Cell2D9Q>>();

        self.mesh = new_timestep;
        self.n_steps += 1;
    }

    pub fn dt(&self) -> f64 {
        self.delta_time
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
            fixed_velocity_boundary(WallBoundary::Left, f, rho_in, u)
        } else if x == x_limit {
            // 0 Pressure case
            // Outlet (right wall)

            // Solve for outlet velocity
            let u_x = 1.0 - (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) / rho_in;
            let u = (u_x, 0.0);

            // pass to boundary condition
            fixed_velocity_boundary(WallBoundary::Right, f, rho_in, u)
        } else if y == 0 {
            // Top Wall
            // Only 2, 5, 6 point back into the fluid
            // noraml is 4
            fixed_velocity_boundary(WallBoundary::Top, f, rho, u_wall)
        } else if y == y_limit {
            // Bottom Wall
            // Only 2, 5, 6 point back into the fluid
            // normal is 2
            fixed_velocity_boundary(WallBoundary::Bottom, f, rho, u_wall)
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
enum WallBoundary {
    Top,
    Bottom,
    Left,
    Right,
}

fn fixed_velocity_boundary(wall: WallBoundary, f: &[f64; 9], rho: f64, u: (f64, f64)) -> [f64; 9] {
    use WallBoundary::*;
    // TODO add in velocity terms for top and bottom boundaries, though they are a non-factor with zero wall velocities
    // This would require further validation of the formulas, something I am not prepared
    // to do right now

    let (u_x, u_y) = u;
    match wall {
        Left => {
            // Left Wall: Only 1, 5, 8 point back into the fluid
            // [
            //     0.0,
            //     f[3] + 2.0 / 3.0 * rho * u_x,
            //     0.0,
            //     0.0,
            //     0.0,
            //     // f[7] + 1.0 / 2.0 * (f[2] - f[4]) - 1.0 / 2.0 * rho * u_y + 1.0 / 6.0 * rho * u_x,
            //     f[7] + 1.0 / 2.0 * (f[2] - f[4]) + 1.0 / 6.0 * rho * u_x,
            //     0.0,
            //     0.0,
            //     // f[6] - 1.0 / 2.0 * (f[2] - f[4]) + 1.0 / 2.0 * rho * u_y + 1.0 / 6.0 * rho * u_x,
            //     f[6] - 1.0 / 2.0 * (f[2] - f[4]) + 1.0 / 6.0 * rho * u_x,
            // ]
            let mut f_bound = *f;
            f_bound[1] = f[3] + 2.0 / 3.0 * rho * u_x;
            f_bound[5] = f[7] + 1.0 / 2.0 * (f[2] - f[4]) + 1.0 / 6.0 * rho * u_x;
            f_bound[8] = f[6] - 1.0 / 2.0 * (f[2] - f[4]) + 1.0 / 6.0 * rho * u_x;

            f_bound
        }
        Top => {
            let mut f_bound = *f;
            f_bound[4] = f[2];
            f_bound[7] = f[5] + 1.0 / 2.0 * (f[1] - f[3]); // + 1.0 / 2.0 * rho * u_x - 1.0 / 6.0 * rho * u_y;
            f_bound[8] = f[6] - 1.0 / 2.0 * (f[1] - f[3]); // - 1.0 / 2.0 * rho * u_x - 1.0 / 6.0 * rho * u_y;
            f_bound
        }
        Right => {
            // Right Wall: Only 1, 5, 8 point back into the fluid
            // [
            //     0.0,
            //     f[3] - 2.0 / 3.0 * rho * u_x,
            //     0.0,
            //     0.0,
            //     0.0,
            //     f[7] + 1.0 / 2.0 * (f[2] - f[4]) - 1.0 / 2.0 * rho * u_y - 1.0 / 6.0 * rho * u_x,
            //     0.0,
            //     0.0,
            //     f[6] - 1.0 / 2.0 * (f[2] - f[4]) + 1.0 / 2.0 * rho * u_y - 1.0 / 6.0 * rho * u_x,
            // ]
            let mut f_bound = *f;
            f_bound[1] = f[3] - 2.0 / 3.0 * rho * u_x;
            f_bound[5] = f[7] + 1.0 / 2.0 * (f[2] - f[4]) - 1.0 / 2.0 * rho * u_y - 1.0 / 6.0 * rho * u_x;
            f_bound[8] = f[6] - 1.0 / 2.0 * (f[2] - f[4]) + 1.0 / 2.0 * rho * u_y - 1.0 / 6.0 * rho * u_x;

            f_bound
        }
        Bottom => {
            // Bottom Wall: Only 2, 5, 6 point back into the fluid
            let mut f_bound = *f;
            f_bound[2] = f[4];
            f_bound[5] = f[7] - 1.0 / 2.0 * (f[1] - f[3]); // + 1.0 / 2.0 * rho * u_x + 1.0 / 6.0 * rho * u_y;
            f_bound[6] = f[8] + 1.0 / 2.0 * (f[1] - f[3]); // - 1.0 / 2.0 * rho * u_x + 1.0 / 6.0 * rho * u_y;

            f_bound
        }
        // _ => panic!("Invalid used of boundary condition function"), // Default: No bounce-back for other directions
    }
}

fn incompressible_eq(rho: f64, u: (f64, f64)) -> [f64; 9] {
    let mut f_eq = [0.0; 9];
    for i in 0..DIRS.len() {
        let t_i = WEIGHTS[i];
        let e = (DIRS[i].0 as f64, DIRS[i].1 as f64);
        let e_dot_u = e.0 * u.0 + e.1 * u.1;
        let u_dot_u = u.0 * u.0 + u.1 * u.1;

        f_eq[i] = nofmt::pls! {
            t_i * (
                rho + 3.0 * (e_dot_u)
                + 9.0 / 2.0 * (e_dot_u.powi(2))
                - 3.0 / 2.0 * (u_dot_u)
            )
        };
    }

    // println!("Incompressible eq: {:?}", f_eq);
    f_eq
}

fn update_cell_properties(cell: &mut Cell2D9Q) {
    let f = cell.f;

    // calculate updated properties
    // sum result & calculate new properties
    let mut rho = 0.0;
    let mut u = (0.0, 0.0);
    for i in 0..9 {
        // F_i(x_i + v_i \Delta t, t + \Delta t) - F_i(x_i, t)
        //      = -\frac{F_i(x_i, t) - F_i^eq(x_i, t)}{\tau}
        // let simga = -(f_stream[i] - f_eq[i]) / self.tau;
        // f[i] = f_eq[i] + f_stream[i];

        // rho = \Sigma F_i
        rho += f[i];

        // rho = \Sigma F_i
        u.0 += f[i] * DIRS[i].0 as f64;
        u.1 += f[i] * DIRS[i].1 as f64;
    }
    
    u = (u.0 / rho, u.1 / rho);
    dbg!(u);

    cell.rho = rho;
    cell.velocity = u;
}
