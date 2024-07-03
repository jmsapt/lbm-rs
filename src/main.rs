
use core::f64;
use std::{
    isize,
    ops::{Index, IndexMut}, usize,
};

mod Sim2D9Q;

use Sim2D9Q::*;






fn main() {
    let time = 5.0; // seconds
    let n_steps = 1000;
    let delta_t = time / (n_steps as f64);

    println!("Building meshes...");
    let sim = LBM2D::new(500, 250, delta_t);


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


