macro_rules! coords {
    ($x:expr, $y:expr) => {
        (0..$x).flat_map(|x| (0..$y).map(move |y| (x, y)))
    };
    ($z:expr) => {
        (0..$z.0).flat_map(|x| (0..$z.1).map(move |y| (x, y)))
    };
}
use core::f64;
use std::{
    isize,
    ops::{Index, IndexMut}, usize,
};

use rand::Rng;

const X_COUNT: usize = 800;
const Y_COUNT: usize = 800;
const L_COUNT: usize = 9;
const TAU: f64 = 0.53;
const RHO_0: f64 = 100.0;
const C: f64 = 1.0;

const WEIGHTS: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 9.0,
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

#[derive(Clone)]
struct Cell2D {
    velocity: (f64, f64),
    // Cell density vector (func f)
    density: [f64; 9],
    rho: f64,
}
impl Cell2D {
    fn new() -> Cell2D {
        Self {
            velocity: (0.0, 0.0),
            density: [0.0; 9],
            rho: 0.0,
        }
    }
}

#[derive(Clone)]
struct Mesh {
    cells: Vec<Cell2D>,
    dimensions: (usize, usize),
}
impl Mesh {
    fn new(x_count: usize, y_count: usize) -> Self {
        let cells = vec![Cell2D::new(); x_count * y_count];

        Self {
            cells,
            dimensions: (x_count, y_count),
        }
    }
}
impl IndexMut<usize> for Mesh {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.cells[index * self.dimensions.1..(index + 1) * self.dimensions.1]
    }
}
impl Index<usize> for Mesh {
    type Output = [Cell2D];

    fn index(&self, index: usize) -> &Self::Output {
        &self.cells[index * self.dimensions.0..(index + 1) * self.dimensions.0]
    }
}

struct LBM2D {
    timesteps: Vec<(f64, Mesh)>,
    geometry: Vec<bool>,
}
impl LBM2D {
    /// Construct and intialise the simulation
    /// TODO: currently hardcoded, make adjustments to parameterise this instead
    fn new(x_count: usize, y_count: usize) -> Self {
        let mut mesh = Mesh::new(x_count, y_count);
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
                geometry[y * x_count + x] = (x as isize - x_0).pow(2) + (y as isize - y_0).pow(2) < r.pow(2);
            })
        }

        // hardcode initial conditions
        coords!(x_count, y_count).for_each(|(x, y)| {
            let mut rng = rand::thread_rng();
            for l in 0..L_COUNT {
                // add random pertubations, and set flow to the right
                mesh[y][x].density[l] = 1.0 + 0.01 * rng.gen::<f64>();
            }
            // TODO
            mesh[y][x].density[3] = 2.0;
        });

        Self {
            timesteps: vec![(0.0, mesh)],
            geometry,
        }
    }

    fn timestep(&mut self, delta_t: f64) {}
}

fn main() {
    let time = 5.0; // seconds
    let n_steps = 1000;

    println!("Building meshes...");

    let mut sim = LBM2D::new(400, 400);

    for i in 0..n_steps {
        let t = time / (n_steps as f64) * i as f64;
        if i % 100 == 0 {
            println!("Time: {t}")
        }

        // collision step

        // streaming step
    }

    // ---- OUTPUT IMAGES ----
    // debug images
    // let mut img = Image::new(X_COUNT as u32, Y_COUNT as u32);
    // for (x, y) in img.coordinates() {
    //     let color = if geometry[y as usize][x as usize] {
    //         WHITE
    //     } else {
    //         BLACK
    //     };

    //     img.set_pixel(x, y, color);
    // }
    // let _ = img.save("geo.bmp");
}

/// Apply collision step to the mesh
///
/// Takes the previous mesh, returning the new mesh.
///
/// => f_i(x, t + delta_t) = f_i(x, t) + f_i^eq(x, t) - f_i(x, t) / tau_f
/// => f_i^eq(x, t) = omega_i * rho * (1 + 3 u_i)
/// where
/// - func f is the density distribution
/// - i \in [0, 8] and is the index for the D2Q9 Lattice
/// - tau_f is characteristic timescale
/// - omega_i is the lattice weight
/// - e_i has same dimension as u (i.e. 2d)
///
/// Note: e_i * u => u_i in this implementation.
fn collision(prev: &Mesh) -> Mesh {
    // refer back to wikipedia article (actually quite good)
    // <https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods>
    let mut new = prev.clone();

    for (x, y) in coords!(prev.dimensions) {
        let e = &DIRS;
        let u = &prev[x][y].velocity;
        let w = &prev[x][y].density;
        let rho = &prev[x][y].rho;

        let collect: [f64; L_COUNT] = (0..L_COUNT)
            .map(|i| {
                let e_dot_u = e[i].0 as f64 * u.0 + e[i].1 as f64 * u.1;
                let u_dot_u = u.0 * u.0 + u.1 * u.1;

                w[i] * rho
                    * (1.0
                        + 3.0 / C.powi(2) * e_dot_u
                        + 9.0 / (2.0 * C.powi(4)) * e_dot_u.powi(2)
                        + 3.0 / (2.0 * C.powi(2)) * u_dot_u.powi(2))
            })
            .collect::<Vec<f64>>()
            .try_into()
            .unwrap();
    }
    todo!()
}

/// Applies the streaming step
///
/// Takes the new mesh (after appling the collision step) and modifies it to stream.
fn streaming(prev: &Mesh, new: &mut Mesh, geometry: &Vec<bool>) {
    for (x, y) in coords!(new.dimensions) {
        for (l) in (0..L_COUNT) {
            let (dx, dy) = DIRS[l];
            let nx = x as isize + dx;
            let ny = y as isize + dy;

            let (x_limit, y_limit) = (new.dimensions.0 as isize, new.dimensions.1 as isize);
            // TODO handle boundary conditions
            if (nx >= 0 && nx < x_limit) && (ny >= 0 && ny < y_limit) {
                // Check geometry
                if geometry[y * x_limit as usize + x] {
                    // flip direction
                    let bounce_l = (l + 4) % 8 + 1;
                    new[x][y].density[l] = prev[x][y].density[bounce_l];
                }

                if (nx == 0) {}

                if (nx == x_limit) {}
            }
        }
    }
}
