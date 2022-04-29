use macroquad::prelude::*;

#[derive(Copy, Clone, Debug, PartialEq)]
struct Planet {
    x: f32,
    y: f32,
    radius: f32,
    mass: f32,
    vx: f32,
    vy: f32,
    color: Color,
}

const PLANET_DENSITY: f32 = 5520.0;
const G: f32 = 6.67408e-8;
// const G: f32 = 6.67408e-11;

fn draw(satellites: &mut Vec<Planet>, planet: &mut Planet) {
    clear_background(BLACK);
    for satellite in &*satellites {
        draw_circle(satellite.x, satellite.y, satellite.radius, satellite.color);
    }

    if is_mouse_button_pressed(MouseButton::Left) {
        planet.x = mouse_position().0;
        planet.y = mouse_position().1;
    }
    if is_mouse_button_down(MouseButton::Left) {
        let x1 = mouse_position().0;
        let y1 = mouse_position().1;
        planet.radius = ((x1 - planet.x).powf(2.) + (y1 - planet.y).powf(2.)).sqrt();
        draw_circle(planet.x, planet.y, planet.radius, planet.color);
    }
    if is_mouse_button_released(MouseButton::Left) {
        planet.mass = planet.radius * planet.radius * planet.radius * PLANET_DENSITY;
        if planet.radius != 0. {
            satellites.push(planet.clone());
        }
    }
}

fn dist(point1: (f32, f32), point2: (f32, f32)) -> f32 {
    ((point1.0 - point2.0).powf(2.) + (point1.1 - point2.1).powf(2.)).sqrt()
}

fn update(satellites: &mut Vec<Planet>) {
    for i in 0..satellites.len() {
        satellites[i].x += satellites[i].vx;
        satellites[i].y += satellites[i].vy;
    }
    for i in 0..satellites.len() {
        for j in 0..satellites.len() {
            if i != j {
                let planet = satellites[i];
                let satellite = satellites[j];

                let d = dist((planet.x, planet.y), (satellite.x, satellite.y));
                if d <= planet.radius + satellite.radius {
                    let new_area = (std::f32::consts::PI * (planet.radius.powf(2.)))
                        + (std::f32::consts::PI * (satellite.radius.powf(2.)));
                    let new_radius = (new_area / std::f32::consts::PI).sqrt();

                    if satellites[i].mass > satellites[j].mass {
                        satellites[i].mass += satellites[j].mass;
                        satellites[i].radius = new_radius;
                        satellites.remove(j);
                    } else {
                        satellites[j].mass += satellites[i].mass;
                        satellites[j].radius = new_radius;
                        satellites.remove(i);
                    }
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

                satellites[i].vx += -fx / planet.mass;
                satellites[i].vy += -fy / planet.mass;

                satellites[j].vx += fx / satellite.mass;
                satellites[j].vy += fy / satellite.mass;
            }
        }
    }
}
#[macroquad::main("Gravity")]
async fn main() {
    let mut planet = Planet {
        x: 0.0,
        y: 0.0,
        radius: 0.0,
        mass: 0.0,
        vx: 0.0,
        vy: 0.0,
        color: PURPLE,
    };
    let mut satellites: Vec<Planet> = Vec::new();
    loop {
        draw(&mut satellites, &mut planet);
        update(&mut satellites);
        next_frame().await;
    }
}
