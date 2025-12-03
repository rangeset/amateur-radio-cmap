use rayon::prelude::*;
use shapefile::{Shape, ShapeReader};
use std::io::Cursor;
use std::iter::Zip;
use std::mem::ManuallyDrop;
use wide::CmpLt;
use wide::{f64x4, f64x8};

const DEGREE_TO_RADIAN_CONSTANT: f64 = std::f64::consts::PI / 180_f64;
const RADIAN_TO_DEGREE_CONSTANT: f64 = 180_f64 / std::f64::consts::PI;
const FRAC_PI_DEGREE: f64x8 = f64x8::new([180f64; 8]);

// TODO C support of ltargcm
//pub fn latlon_to_azimnth_isometric_csupport() {}
pub fn latlon_to_azimnth_isometric(latitude_delta: f64, longitude_delta: f64) -> (f64, f64) {
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

// TODO used wide to accelerate multipoint,line,polygon
//pub fn latlon_to_azimnth_isometric_simd_csupport() {}

pub fn latlon_to_azimnth_isometric_simd(
    latitude_delta_vec: Vec<f64>,
    longitude_delta_vec: Vec<f64>,
) -> (Vec<f64>, Vec<f64>) {
    let len = std::cmp::min(latitude_delta_vec.len(), longitude_delta_vec.len());
    let latitude_delta_vec_chunks = latitude_delta_vec[..len].par_chunks_exact(8);
    let longitude_delta_vec_chunks = longitude_delta_vec[..len].par_chunks_exact(8);
    let chunks_remainder = latitude_delta_vec_chunks
        .remainder()
        .iter()
        .copied()
        .zip(longitude_delta_vec_chunks.remainder().iter().copied());
    let chunks_bodys: Vec<([f64; 8], [f64; 8])> = latitude_delta_vec_chunks
        .zip(longitude_delta_vec_chunks)
        .map(
            |(latitude_slice, longitude_slice): (&[f64], &[f64])| -> ([f64; 8], [f64; 8]) {
                let latitude_array = latitude_slice;
                let longitude_array = longitude_slice;

                let square = |x: f64x8| x * x;

                let latitude = f64x8::from(latitude_array);
                let longitude = f64x8::from(longitude_array);

                let posistion_array: [f64; 8] = longitude.abs().to_array();

                let posistion = longitude.abs().simd_lt(f64x8::FRAC_PI_2);

                let latitude_radian = latitude.to_radians();
                let longitude_radian = longitude.to_radians();

                let latitude_radian_sine = latitude_radian.sin();
                let longitude_radian_sine = longitude_radian.sin();

                let distance =
                    (square(latitude_radian_sine) + square(longitude_radian_sine)).sqrt();

                let distance_arcsine = distance.asin().to_degrees();
                let distance_arcsine_back = FRAC_PI_DEGREE - distance_arcsine;

                let k_uncheck = posistion.blend(distance_arcsine, distance_arcsine_back) / distance;

                let k_infinite_mask = k_uncheck.is_inf();
                let k_nan_mask = k_uncheck.is_nan();

                let k_bad = k_infinite_mask | k_nan_mask;

                let k = k_bad.blend(f64x8::ZERO, k_uncheck);

                ((latitude * k).to_array(), (longitude * k).to_array())
            },
        )
        .collect();

    let mut xs: Vec<f64> = Vec::with_capacity(len);
    let mut ys: Vec<f64> = Vec::with_capacity(len);

    for (x, y) in chunks_bodys {
        xs.extend(x);
        ys.extend(y);
    }

    for (latitude, longitude) in chunks_remainder {
        let result = latlon_to_azimnth_isometric(latitude, longitude);
        xs.push(result.0);
        ys.push(result.1);
    }

    (xs, ys)
}

pub struct ReturnContent {
    pub status: bool,
    pub ptr: *const u8,
    pub len: usize,
}

impl ReturnContent {
    fn new(data: Vec<u8>, status: bool) -> Self {
        let data = ManuallyDrop::new(data);
        ReturnContent {
            status: status,
            ptr: data.as_ptr(),
            len: data.len(),
        }
    }
}

pub struct ColorData {
    pub color_point: u8,
    pub color_multipoint: u8,
    pub color_line: u8,
    pub color_polygon_line: u8,
    pub color_polygon_fill: u8,
}

pub fn shapefile_generate(
    buffer_ptr: *const u8,
    buffer_len: usize,
    color: ColorData,
) -> ReturnContent {
    let buffer = unsafe { std::slice::from_raw_parts(buffer_ptr, buffer_len) };
    let cursor = Cursor::new(buffer);
    let mut reader = match ShapeReader::new(cursor) {
        Ok(data) => data,
        Err(e) => return ReturnContent::new(e.to_string().into_bytes(), false),
    };
    // TODO Draw the picture used par_iter
    //    reader.for_each(|shape| )
    ReturnContent::new(Vec::new(), true)
}

// TODO Draw fuction
//fn shapefile_draw() -> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transprojection() {
        let result = latlon_to_azimnth_isometric(0_f64, 0_f64);
        assert_eq!(result, (0_f64, 0_f64));

        let result = latlon_to_azimnth_isometric(-90_f64, 0_f64);
        assert_eq!(result, (0_f64, -90_f64));

        let result = latlon_to_azimnth_isometric(0_f64, 180_f64);
        assert_eq!(result, (180_f64, 0_f64));
    }
}
