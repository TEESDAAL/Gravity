use macroquad::prelude::*;
use macroquad::time;
#[derive(Copy, Clone, Debug, PartialEq)]
struct Planet {
    x: f32,
    y: f32,
    radius: f32,
    mass: f32,
    vx: f32,
    vy: f32,
    color: Color,
    fixed: bool,
}

const PLANET_DENSITY: f32 = 5520.0; // kg/m^3
const G: f32 = 6.67408e-11; // m^3/(kg*s^2)
const SECONDS_PER_TICK: f32 = 1000. * 15.; // s

#[derive(PartialEq)]
enum Status {
    Normal,
    Planet,
    Velocity,
}

fn draw(satellites: &mut Vec<Planet>, planet: &mut Planet, status: &mut Status, scale: &f32) {
    if is_key_pressed(KeyCode::F) {
        planet.fixed = !planet.fixed;
    }
    if *status == Status::Normal {
        // If the mouse is set the centre of the circle
        if is_mouse_button_pressed(MouseButton::Left) {
            planet.x = mouse_position().0;
            planet.y = mouse_position().1;
            *status = Status::Planet;
        }
    }

    if *status == Status::Planet {
        // Allow the user to set the radius by the position of the mouse
        if is_mouse_button_down(MouseButton::Left) {
            let x1 = mouse_position().0;
            let y1 = mouse_position().1;
            planet.radius = ((x1 - planet.x).powf(2.) + (y1 - planet.y).powf(2.)).sqrt();
            draw_circle(planet.x, planet.y, planet.radius, planet.color);
        }

        // Set the radius of the planet to the mouse position when the mouse is released
        if is_mouse_button_released(MouseButton::Left) {
            planet.mass = planet.radius * planet.radius * planet.radius * PLANET_DENSITY;
            *status = Status::Velocity;
        }
    }
    if *status == Status::Velocity {
        if planet.fixed {
            satellites.push(planet.clone());
            *status = Status::Normal;
            return;
        }
        clear_background(BLACK);
        draw_circle(planet.x, planet.y, planet.radius, planet.color);
        let (x, y) = mouse_position();
        draw_line(planet.x, planet.y, x, y, 10., WHITE);
        if is_mouse_button_pressed(MouseButton::Left) {
            let velocity = dist((planet.x, planet.y), (x, y)) / 100.;
            let mut angle = (y - planet.y).atan2(x - planet.x);
            if angle < 0. {
                angle += 2. * std::f32::consts::PI;
            }
            planet.vx = velocity * angle.cos();
            planet.vy = velocity * angle.sin();

            if planet.radius != 0. {
                satellites.push(planet.clone());
            }
            *status = Status::Normal;
        }
    }
    // Draw all the planets
    for satellite in satellites {
        draw_circle(
            (satellite.x) * scale,
            (satellite.y) * scale,
            satellite.radius * scale,
            satellite.color,
        );
    }
}

fn dist(point1: (f32, f32), point2: (f32, f32)) -> f32 {
    ((point1.0 - point2.0).powf(2.) + (point1.1 - point2.1).powf(2.)).sqrt()
}

fn update(satellites: &mut Vec<Planet>, paused: &mut bool, scale: &mut f32) {
    if is_key_pressed(KeyCode::Space) {
        *paused = !*paused;
    }

    if is_key_pressed(KeyCode::C) {
        satellites.clear();
        return;
    }
    if is_key_down(KeyCode::LeftControl) {
        if is_key_pressed(KeyCode::Equal) {
            *scale *= 1.1;
        }
        if is_key_pressed(KeyCode::Minus) {
            *scale /= 1.1;
        }
    }
    if *paused {
        return;
    }
    for i in 0..satellites.len() {
        satellites[i].x += satellites[i].vx;
        satellites[i].y += satellites[i].vy;
    }
    for i in 0..satellites.len() {
        for j in 0..satellites.len() {
            if i != j {
                let planet = satellites[i];
                let satellite = satellites[j];
                // println!("{}", satellite.fixed);

                let d = dist((planet.x, planet.y), (satellite.x, satellite.y));
                if d <= planet.radius + satellite.radius {
                    let new_area = (std::f32::consts::PI * (planet.radius.powf(2.)))
                        + (std::f32::consts::PI * (satellite.radius.powf(2.)));
                    let new_radius = (new_area / std::f32::consts::PI).sqrt();

                    let new_x_vel = (planet.mass * planet.vx + satellite.mass * satellite.vx)
                        / (planet.mass + satellite.mass);

                    let new_y_vel = (planet.mass * planet.vy + satellite.mass * satellite.vy)
                        / (planet.mass + satellite.mass);

                    let (larger, smaller) = match satellites[i].mass > satellites[j].mass {
                        true => (i, j),
                        false => (j, i),
                    };
                    if satellites[smaller].fixed {
                        satellites[larger].fixed = true;
                    }
                    satellites[larger].mass += satellites[smaller].mass;
                    satellites[larger].radius = new_radius;
                    if !satellites[larger].fixed {
                        satellites[larger].vx = new_x_vel;
                        satellites[larger].vy = new_y_vel;
                    }
                    satellites.remove(smaller);
                    return;
                }

                let gravitational_force = (G * satellite.mass * planet.mass) / d.powf(2.);

                let mut angle_between_planets =
                    (planet.y - satellite.y).atan2(planet.x - satellite.x);

                if angle_between_planets < 0. {
                    angle_between_planets += 2. * std::f32::consts::PI
                }

                let fy = gravitational_force * angle_between_planets.sin();
                let fx = gravitational_force * angle_between_planets.cos();

                satellites[i].vx += -fx * SECONDS_PER_TICK / planet.mass;
                satellites[i].vy += -fy * SECONDS_PER_TICK / planet.mass;

                satellites[j].vx += fx * SECONDS_PER_TICK / satellite.mass;
                satellites[j].vy += fy * SECONDS_PER_TICK / satellite.mass;

                if satellites[i].fixed {
                    satellites[i].vx = 0.;
                    satellites[i].vy = 0.;
                }
                if satellites[j].fixed {
                    satellites[j].vx = 0.;
                    satellites[j].vy = 0.;
                }
            }
        }
    }
}

#[macroquad::main("Gravity")]
async fn main() {
    println!("{}", time::get_fps());

    let mut scale = 1.;
    let mut paused = false;
    let mut planet = Planet {
        x: 0.0,
        y: 0.0,
        radius: 0.0,
        mass: 0.0,
        vx: 0.0,
        vy: 0.0,
        color: PURPLE,
        fixed: false,
    };
    let mut status = Status::Normal;
    let mut satellites: Vec<Planet> = Vec::new();

    // set_target_fps(1000.);
    loop {
        // let screen_centre = (screen_width() / 2., screen_height() / 2.);
        draw(&mut satellites, &mut planet, &mut status, &mut scale);
        update(&mut satellites, &mut paused, &mut scale);
        next_frame().await;
    }
}
