use core::f64;
use plotters::prelude::*;
use std::error::Error;
mod Sim2D9Q;

use Sim2D9Q::*;

const OUTPUT_FILE: &str = "plot.png";

fn main() -> Result<(), Box<dyn Error>> {
    let time = 5.0; // seconds
    let n_steps = 1;
    let _delta_t = time / (n_steps as f64);

    println!("Building meshes...");
    let dx = 2.0e-3; // m
    let c_air = 100.0; // m/s
    let rho_air = 1.0; // kg/m^3
    let u_inlet = (0.010, 0.0); // m/s
    let size = (10, 10);
    let mut sim = LBM2D::new(size.0, size.1, dx, c_air, u_inlet, rho_air);

    // Initial state
    for x in 0..size.0 {
        for y in 0..size.1 {
            let u = sim.get_velocity(x, y);
            print!("{:.3} ", u.0);
        }
        println!("");
    }
    println!("-------");

    for i in 0..n_steps {
        // if i % (n_steps / 10) == 0 {
        //     println!("{}", i as f64 * sim.dt());
        // }
        sim.step();
    }
    for x in 0..size.0 {
        for y in 0..size.1 {
            let u = sim.get_velocity(x, y);
            print!("{:.3} ", u.0);
        }
        println!();
    }

    println!("Timestep size: {}", sim.delta_time);
    println!("Tau: {}", sim.tau);
    return Ok(());

    let root = BitMapBackend::new(OUTPUT_FILE, (1000, 1000)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(10)
        .y_label_area_size(10)
        .build_cartesian_2d(0f64..size.0 as f64 * dx, 0f64..size.1 as f64 * dx)?;

    chart
        .configure_mesh()
        // .disable_x_mesh()
        // .disable_y_mesh()
        .draw()?;

    let plotting_area = chart.plotting_area();

    let u_inlet_mag = (u_inlet.0.powi(2) + u_inlet.1.powi(2)).sqrt();
    for x in 0..size.0 {
        for y in 0..size.1 {
            let pos = (x as f64 * dx, y as f64 * dx);

            let u = sim.get_velocity(x, y);
            // print!("{:.3} ", u.0);
            let u_mag = (u.0.powi(2) + u.1.powi(2)).sqrt();
            let x = (255_f64 * (u_mag / u_inlet_mag)) as u8;

            let colour = RGBColor(x, x, x).filled();
            plotting_area.draw(&Rectangle::new(
                [(pos.0, pos.1), (pos.0 + dx, pos.1 + dx)],
                colour,
            ))?;
        }
        println!();
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
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    Ok(())
}
