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

fn scaled_coords(coord: &f32, scale: &f32) -> f32 {
    (coord) * scale + ((screen_width() / 2.) * (1. - scale))
}

fn _look_ahead(
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

            if d < planet.radius + satellite.radius && *collisions_on {
                return;
            }
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
        if location == (planet.x, planet.y) {
            return;
        }
    }
}

fn _true_look_ahead(planets: &mut Vec<Planet>, planet: &Planet, collisions_on: &bool) {
    planets.push(*planet);
    let index = planets.len() - 1;
    for _ in 0..LOOK_AHEAD {
        let location = (planets[index].x, planets[index].y);
        for i in 0..planets.len() {
            let planet = planets[i];
            if planet.fixed {
                continue;
            }
            let (mut fx, mut fy) = (0., 0.);
            for j in 0..planets.len() {
                if i == j {
                    continue;
                }

                let satellite = planets[j];

                let d = dist((planet.x, planet.y), (satellite.x, satellite.y));
                if d <= planet.radius + satellite.radius && *collisions_on {
                    draw_circle(planet.x, planet.y, 10., RED);
                    return;
                }

                let gravitational_force = (G * satellite.mass * planet.mass) / (d * d);

                let mut angle_between_planets =
                    atan2(planet.y - satellite.y, planet.x - satellite.x);

                if angle_between_planets < 0. {
                    angle_between_planets += 2. * PI
                }
                let (sin, cos) = sincos(angle_between_planets as f64);
                fy += gravitational_force * sin as f32;
                fx += gravitational_force * cos as f32;
            }
            planets[i].vx += -fx * SECONDS_PER_TICK / planets[i].mass;
            planets[i].vy += -fy * SECONDS_PER_TICK / planets[i].mass;
        }
        for i in 0..planets.len() {
            if planets[i].fixed {
                planets[i].vx = 0.;
                planets[i].vy = 0.;
            } else {
                planets[i].x += planets[i].vx;
                planets[i].y += planets[i].vy;
            }
        }
        draw_line(
            location.0,
            location.1,
            planets[index].x,
            planets[index].y,
            2.,
            match collisions_on {
                true => WHITE,
                false => RED,
            },
        );
    }
}

fn true_detailed_look_ahead(
    planets: &Vec<Planet>,
    planet: &Planet,
    collisions_on: &bool,
    scale: &f32,
) {
    let mut planets = planets.clone();
    planets.push(*planet);
    for _ in 0..LOOK_AHEAD {
        'outer: for i in 0..planets.len() {
            let planet = planets[i];
            let location = (planet.x, planet.y);

            if planet.fixed {
                continue;
            }
            let (mut fx, mut fy) = (0., 0.);
            for j in 0..planets.len() {
                if i == j {
                    continue;
                }

                let satellite = planets[j];

                let d = dist((planet.x, planet.y), (satellite.x, satellite.y));
                if d <= planet.radius + satellite.radius && *collisions_on {
                    draw_circle(planet.x, planet.y, 10., RED);
                    planet_collision(&planet, &satellite, &mut planets, &i, &j);
                    break 'outer;
                }

                let gravitational_force = (G * satellite.mass * planet.mass) / (d * d);

                let mut angle_between_planets =
                    atan2(planet.y - satellite.y, planet.x - satellite.x);

                if angle_between_planets < 0. {
                    angle_between_planets += 2. * PI
                }
                let (sin, cos) = sincos(angle_between_planets as f64);
                fy += gravitational_force * sin as f32;
                fx += gravitational_force * cos as f32;
            }
            planets[i].vx += -fx * SECONDS_PER_TICK / planets[i].mass;
            planets[i].vy += -fy * SECONDS_PER_TICK / planets[i].mass;
            if planets[i].fixed {
                planets[i].vx = 0.;
                planets[i].vy = 0.;
            } else {
                planets[i].x += planets[i].vx;
                planets[i].y += planets[i].vy;
            }
            draw_line(
                scaled_coords(&location.0, scale),
                scaled_coords(&location.1, scale),
                scaled_coords(&planets[i].x, scale),
                scaled_coords(&planets[i].y, scale),
                2.,
                match collisions_on {
                    true => WHITE,
                    false => RED,
                },
            );
        }
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

    let (x, y) = mouse_position();
    let real_x = (x - ((screen_width() / 2.) * (1. - scale))) / scale;
    let real_y = (y - ((screen_height() / 2.) * (1. - scale))) / scale;

    match *status {
        Status::Normal => {
            // If the mouse is set the centre of the circle
            if is_mouse_button_pressed(MouseButton::Left) {
                planet.x = real_x;
                planet.y = real_y;
                *status = Status::Planet;
            }
        }

        Status::Planet => {
            // Allow the user to set the radius by the position of the mouse
            if is_mouse_button_down(MouseButton::Left) {
                planet.radius =
                    ((real_x - planet.x).powf(2.) + (real_y - planet.y).powf(2.)).sqrt() * scale;
                draw_circle(
                    (planet.x) * scale + ((screen_width() / 2.) * (1. - scale)),
                    (planet.y) * scale + ((screen_height() / 2.) * (1. - scale)),
                    planet.radius,
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
                (planet.x) * scale + ((screen_width() / 2.) * (1. - scale)),
                (planet.y) * scale + ((screen_height() / 2.) * (1. - scale)),
                planet.radius,
                planet.color,
            );

            let velocity = dist((planet.x, planet.y), (x, y)) / (100. * scale);
            let angle = atan2(y - planet.y, x - planet.x);
            let (sin, cos) = sincos(angle as f64);
            planet.vx = velocity * cos as f32;
            planet.vy = velocity * sin as f32;

            true_detailed_look_ahead(&satellites, planet, collisions_on, scale);

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
            scaled_coords(&satellite.x, scale),
            scaled_coords(&satellite.y, scale),
            satellite.radius * scale,
            match collisions_on {
                true => satellite.color,
                false => BLACK,
            },
        );
    }
}

fn update(planets: &mut Vec<Planet>, paused: &bool, collisions_on: &bool) {
    if *paused {
        return;
    }
    'outer: for i in 0..planets.len() {
        let planet = planets[i];

        if planet.fixed {
            continue;
        }
        let (mut fx, mut fy) = (0., 0.);
        for j in 0..planets.len() {
            if i == j {
                continue;
            }

            let satellite = planets[j];

            let d = dist((planet.x, planet.y), (satellite.x, satellite.y));
            if d <= planet.radius + satellite.radius && *collisions_on {
                draw_circle(planet.x, planet.y, 10., RED);
                planet_collision(&planet, &satellite, planets, &i, &j);
                break 'outer;
            }

            let gravitational_force = (G * satellite.mass * planet.mass) / (d * d);

            let angle_between_planets = atan2(planet.y - satellite.y, planet.x - satellite.x);

            let (sin, cos) = sincos(angle_between_planets as f64);
            fy += gravitational_force * sin as f32;
            fx += gravitational_force * cos as f32;
        }
        planets[i].vx += -fx * SECONDS_PER_TICK / planets[i].mass;
        planets[i].vy += -fy * SECONDS_PER_TICK / planets[i].mass;
        if planets[i].fixed {
            planets[i].vx = 0.;
            planets[i].vy = 0.;
        } else {
            planets[i].x += planets[i].vx;
            planets[i].y += planets[i].vy;
        }
    }
}

#[macroquad::main("Gravity")]
async fn main() {
    let mut scale = 1.;
    let mut paused: bool = false;
    let mut collisions_on = true;
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
