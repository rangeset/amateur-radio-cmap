const DEGREE_TO_RADIAN_CONSTANT: f64 = std::f64::consts::PI / 180_f64;
const RADIAN_TO_DEGREE_CONSTANT: f64 = 180_f64 / std::f64::consts::PI;

pub fn latlon_to_amateur_radio_great_circle_map(
    latitude_delta: f64,
    longitude_delta: f64,
) -> (f64, f64) {
    let square = |x: f64| x * x;
    let degree_to_radian = |degree: f64| degree * DEGREE_TO_RADIAN_CONSTANT;
    let radian_to_degree = |radian: f64| radian * RADIAN_TO_DEGREE_CONSTANT;
    let latitude_delta_radian: f64 = degree_to_radian(latitude_delta);
    let longitude_delta_radian: f64 = degree_to_radian(longitude_delta);
    let hemispheres_anterior_or_posterior: bool = longitude_delta.abs() < 90_f64;

    #[cfg(test)]
    dbg!(latitude_delta_radian, longitude_delta_radian);

    let latitude_delta_radian_sine: f64 = latitude_delta_radian.sin();
    let longitude_delta_radian_sine: f64 = longitude_delta_radian.sin();

    #[cfg(test)]
    dbg!(latitude_delta_radian_sine, longitude_delta_radian_sine);

    let distance: f64 =
        (square(latitude_delta_radian_sine) + square(longitude_delta_radian_sine)).sqrt();
    let k: f64 = match distance {
        0_f64 => 0_f64,
        _ => {
            (if hemispheres_anterior_or_posterior {
                radian_to_degree(distance.asin())
            } else {
                180_f64 - radian_to_degree(distance.asin())
            }) / distance
        }
    };

    #[cfg(test)]
    dbg!(distance, k);

    (
        longitude_delta_radian_sine * k,
        latitude_delta_radian_sine * k,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transprojection() {
        let result = latlon_to_amateur_radio_great_circle_map(0_f64, 0_f64);
        assert_eq!(result, (0_f64, 0_f64));

        let result = latlon_to_amateur_radio_great_circle_map(-90_f64, 0_f64);
        assert_eq!(result, (0_f64, -90_f64));

        let result = latlon_to_amateur_radio_great_circle_map(0_f64, 180_f64);
        assert_eq!(result, (180_f64, 0_f64));
    }
}
