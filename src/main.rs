use fast_math::atan2;
use libm::sincos;
use macroquad::prelude::*;
use std::f32::consts::PI;

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
const LOOK_AHEAD: i32 = 10000;
const FIXED_COLOR: Color = RED;
const MOVING_COLOR: Color = PURPLE;

enum Status {
    Normal,
    Planet,
    Velocity,
}

fn look_ahead(
    planets: &Vec<Planet>,
    planet: &Planet,
    x_vel: &mut f32,
    y_vel: &mut f32,
    collisions_on: &bool,
) {
    let mut location = (planet.x, planet.y);
    for _ in 0..LOOK_AHEAD {
        let (mut fx, mut fy) = (0., 0.);
        for satellite in planets {
            if satellite == planet {
                continue;
            }
            let d = dist((location.0, location.1), (satellite.x, satellite.y));
            let gravitational_force = (G * satellite.mass * planet.mass) / (d * d);

            let mut angle_between_planets =
                atan2(location.1 - satellite.y, location.0 - satellite.x);

            if angle_between_planets < 0. {
                angle_between_planets += 2. * PI
            }

            let (sin, cos) = sincos(angle_between_planets as f64);
            fy += gravitational_force * sin as f32;
            fx += gravitational_force * cos as f32;
        }

        *x_vel += -fx * SECONDS_PER_TICK / planet.mass;
        *y_vel += -fy * SECONDS_PER_TICK / planet.mass;
        draw_line(
            location.0,
            location.1,
            location.0 + *x_vel,
            location.1 + *y_vel,
            2.,
            match collisions_on {
                true => WHITE,
                false => FIXED_COLOR,
            },
        );
        location = (location.0 + *x_vel, location.1 + *y_vel);
    }
}

fn planet_collision(
    planet1: &Planet,
    planet2: &Planet,
    planets: &mut Vec<Planet>,
    i: &usize,
    j: &usize,
) {
    let new_area = (PI * (planet1.radius.powf(2.))) + (PI * (planet2.radius.powf(2.)));
    let new_radius = (new_area / PI).sqrt();

    let new_x_vel =
        (planet1.mass * planet1.vx + planet2.mass * planet2.vx) / (planet1.mass + planet2.mass);

    let new_y_vel =
        (planet1.mass * planet1.vy + planet2.mass * planet2.vy) / (planet1.mass + planet2.mass);

    let (larger, smaller) = match planet1.mass > planet2.mass {
        true => (*i, *j),
        false => (*j, *i),
    };
    if planets[smaller].fixed {
        planets[larger].fixed = true;
        planets[larger].color = FIXED_COLOR;
    }
    planets[larger].mass += planets[smaller].mass;
    planets[larger].radius = new_radius;
    if !planets[larger].fixed {
        planets[larger].vx = new_x_vel;
        planets[larger].vy = new_y_vel;
    }
    planets.remove(smaller);
}

fn dist(point1: (f32, f32), point2: (f32, f32)) -> f32 {
    let (x_dist, y_dist) = (point1.0 - point2.0, point1.1 - point2.1);
    (x_dist * x_dist + y_dist * y_dist).sqrt()
}

fn key_events(
    planets: &mut Vec<Planet>,
    planet: &mut Planet,
    paused: &mut bool,
    scale: &mut f32,
    collisions_on: &mut bool,
) {
    if is_key_pressed(KeyCode::F) {
        planet.fixed = !planet.fixed;
    }
    if is_key_pressed(KeyCode::Space) {
        *paused = !*paused;
    }

    if is_key_pressed(KeyCode::D) {
        *collisions_on = !*collisions_on;
    }

    if is_key_pressed(KeyCode::C) {
        planets.clear();
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
}

fn draw(
    satellites: &mut Vec<Planet>,
    planet: &mut Planet,
    status: &mut Status,
    scale: &f32,
    collisions_on: &bool,
) {
    if planet.fixed {
        planet.color = FIXED_COLOR;
    } else {
        planet.color = MOVING_COLOR;
    }
    match *status {
        Status::Normal => {
            // If the mouse is set the centre of the circle
            if is_mouse_button_pressed(MouseButton::Left) {
                planet.x = mouse_position().0;
                planet.y = mouse_position().1;
                *status = Status::Planet;
            }
        }

        Status::Planet => {
            // Allow the user to set the radius by the position of the mouse
            if is_mouse_button_down(MouseButton::Left) {
                let x1 = mouse_position().0 * scale;
                let y1 = mouse_position().1 * scale;
                planet.radius = ((x1 - planet.x).powf(2.) + (y1 - planet.y).powf(2.)).sqrt();
                draw_circle(
                    (planet.x - (screen_width() / 2.) * (1. - scale)) / scale,
                    (planet.y - (screen_height() / 2.) * (1. - scale)) / scale,
                    planet.radius / scale,
                    match collisions_on {
                        true => planet.color,
                        false => BLACK,
                    },
                );
            }

            // Set the radius of the planet to the mouse position when the mouse is released
            if is_mouse_button_released(MouseButton::Left) {
                planet.mass = planet.radius * planet.radius * planet.radius * PLANET_DENSITY;
                *status = Status::Velocity;
            }
        }

        Status::Velocity => {
            if planet.fixed {
                satellites.push(planet.clone());
                *status = Status::Normal;
                return;
            }
            draw_circle(
                planet.x * scale,
                planet.y * scale,
                planet.radius * scale,
                planet.color,
            );
            let (x, y) = mouse_position();
            draw_line(
                planet.x * scale,
                planet.y * scale,
                x,
                y,
                10.,
                match collisions_on {
                    true => WHITE,
                    false => BLACK,
                },
            );

            let velocity = dist((planet.x, planet.y), (x, y)) / (100. * scale);
            let mut angle = atan2(y - planet.y, x - planet.x);
            if angle < 0. {
                angle += 2. * PI;
            }
            let (sin, cos) = sincos(angle as f64);
            planet.vx = velocity * cos as f32;
            planet.vy = velocity * sin as f32;
            look_ahead(
                satellites,
                planet,
                &mut planet.vx.clone(),
                &mut planet.vy.clone(),
                collisions_on,
            );

            if is_mouse_button_pressed(MouseButton::Left) {
                if planet.radius != 0. {
                    satellites.push(planet.clone());
                }
                *status = Status::Normal;
            }
        }
    }

    // Draw all the planets
    for satellite in satellites {
        draw_circle(
            (satellite.x) * scale + ((screen_width() / 2.) * (1. - scale)),
            (satellite.y) * scale + ((screen_height() / 2.) * (1. - scale)),
            satellite.radius * scale,
            match collisions_on {
                true => satellite.color,
                false => BLACK,
            },
        );
    }
}

fn update(satellites: &mut Vec<Planet>, paused: &bool, collisions_on: &bool) {
    if *paused {
        return;
    }
    for i in 0..satellites.len() {
        let planet = satellites[i];
        if planet.fixed {
            continue;
        }
        let (mut fx, mut fy) = (0., 0.);
        for j in 0..satellites.len() {
            if i == j {
                continue;
            }

            let satellite = satellites[j];

            let d = dist((planet.x, planet.y), (satellite.x, satellite.y));
            if d <= planet.radius + satellite.radius && *collisions_on {
                planet_collision(&planet, &satellite, satellites, &i, &j);
                return;
            }

            let gravitational_force = (G * satellite.mass * planet.mass) / (d * d);

            let mut angle_between_planets = atan2(planet.y - satellite.y, planet.x - satellite.x);

            if angle_between_planets < 0. {
                angle_between_planets += 2. * PI
            }
            let (sin, cos) = sincos(angle_between_planets as f64);
            fy += gravitational_force * sin as f32;
            fx += gravitational_force * cos as f32;
        }
        satellites[i].vx += -fx * SECONDS_PER_TICK / satellites[i].mass;
        satellites[i].vy += -fy * SECONDS_PER_TICK / satellites[i].mass;
    }
    for i in 0..satellites.len() {
        if satellites[i].fixed {
            satellites[i].vx = 0.;
            satellites[i].vy = 0.;
        } else {
            satellites[i].x += satellites[i].vx;
            satellites[i].y += satellites[i].vy;
        }
    }
}

#[macroquad::main("Gravity")]
async fn main() {
    let mut scale = 1.;
    let mut paused = false;
    let mut collisions_on: bool = true;
    let mut planet = Planet {
        x: 0.0,
        y: 0.0,
        radius: 0.0,
        mass: 0.0,
        vx: 0.0,
        vy: 0.0,
        color: MOVING_COLOR,
        fixed: false,
    };
    let mut status = Status::Normal;
    let mut satellites: Vec<Planet> = Vec::new();

    loop {
        match collisions_on {
            true => clear_background(BLACK),
            false => clear_background(WHITE),
        }
        key_events(
            &mut satellites,
            &mut planet,
            &mut paused,
            &mut scale,
            &mut collisions_on,
        );
        update(&mut satellites, &paused, &collisions_on);
        draw(
            &mut satellites,
            &mut planet,
            &mut status,
            &mut scale,
            &collisions_on,
        );
        next_frame().await;
    }
}
