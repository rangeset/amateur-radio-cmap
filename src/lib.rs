use image::GenericImage;
use image::codecs::png::PngEncoder;
use image::{ColorType, ImageEncoder};
use imageproc::drawing::{
    draw_filled_circle_mut, draw_line_segment_mut, draw_polygon, draw_polygon_mut,
};
use imageproc::image::{GrayImage, Rgba, RgbaImage};
use rayon::prelude::*;
use shapefile::PolygonRing::{Inner, Outer};
use shapefile::record::{
    multipoint::GenericMultipoint, polygon::GenericPolygon, polyline::GenericPolyline,
    traits::HasXY,
};
use shapefile::{Shape, ShapeReader};
use std::collections::HashSet;
use std::io::Cursor;
use std::mem::ManuallyDrop;

const CHUNKS_SIZE: usize = 4096;
const INSERT_LIMIT: isize = 256;
const INSERT_THRESHOLD: f64 = 0.1;
const INSERT_DENSITY: f64 = 0.01;

const DEGREE_TO_RADIAN_CONSTANT: f64 = std::f64::consts::PI / 180_f64;
const RADIAN_TO_DEGREE_CONSTANT: f64 = 180_f64 / std::f64::consts::PI;
const PI: f64 = std::f64::consts::PI;

#[warn(no_mangle_generic_items)]
pub extern "C" fn free<T>(ptr: *const T, len: usize) {
    unsafe {
        let _ = Box::from_raw(std::slice::from_raw_parts_mut(ptr as *mut T, len));
    }
}

#[repr(C)]
pub struct LatlonToAzimnthIsometricCsupport {
    pub x: f64,
    pub y: f64,
}

#[unsafe(no_mangle)]
pub extern "C" fn latlon_to_azimuthal_equidistant_csupport(
    lat: f64,
    lat_base: f64,
    lon_delta: f64,
) -> LatlonToAzimnthIsometricCsupport {
    let (x, y, _): (f64, f64, f64) = latlon_to_azimuthal_equidistant(lat, lon_delta, lat_base);
    LatlonToAzimnthIsometricCsupport { x, y }
}

pub fn latlon_to_azimuthal_equidistant(lat: f64, lon_delta: f64, lat_base: f64) -> (f64, f64, f64) {
    let square = |x: f64| x * x;
    let square_overflow = |x: f64| square(x).clamp(-1f64, 1f64);
    let lat_rad = lat * DEGREE_TO_RADIAN_CONSTANT;
    let lon_delta_rad = lon_delta * DEGREE_TO_RADIAN_CONSTANT;
    let (lat_sin, lat_cos) = lat_rad.sin_cos();
    let (lon_delta_sin, lon_delta_cos) = lon_delta_rad.sin_cos();
    let p_x = lat_cos * lon_delta_cos;
    let p_y = lat_sin;
    let p_angle = p_y.atan2(p_x);
    let p_distance = (square(p_x) + square(p_y)).sqrt();
    let angle_delta = p_angle - lat_base * DEGREE_TO_RADIAN_CONSTANT;
    let (angle_delta_sin, angle_delta_cos) = angle_delta.sin_cos();
    let x = angle_delta_cos * p_distance;
    if x.abs() == 1f64 {
        return (0f64, 0f64, x);
    }
    let y = angle_delta_sin * p_distance;
    let distance = x.acos();
    let angle_cos = y / (1f64 - square(x)).sqrt();
    let angle_sin = if lon_delta_sin > 0f64 {
        (1f64 - square_overflow(angle_cos)).sqrt()
    } else {
        -(1f64 - square_overflow(angle_cos)).sqrt()
    };
    (angle_sin * distance, -angle_cos * distance, x)
}

#[derive(Copy, Clone)]
struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub fn new<T: HasXY>(p: &T) -> Self {
        Point { x: p.x(), y: p.y() }
    }
}

fn latlon_to_azimuthal_equidistant_array(
    points: Vec<Point>,
    lat_base: f64,
    lon_base: f64,
) -> (Vec<f64>, Vec<f64>) {
    let points_len = points.len();
    let (xs, ys, ks): (Vec<f64>, Vec<f64>, Vec<f64>) = points
        .par_chunks(CHUNKS_SIZE)
        .map(|chunk| {
            chunk
                .iter()
                .map(|point| latlon_to_azimuthal_equidistant(point.y, point.x - lon_base, lat_base))
                .fold(
                    (
                        Vec::with_capacity(chunk.len()),
                        Vec::with_capacity(chunk.len()),
                        Vec::with_capacity(chunk.len()),
                    ),
                    |(mut xs, mut ys, mut ks), (x, y, k)| {
                        xs.push(x);
                        ys.push(y);
                        ks.push(k);
                        (xs, ys, ks)
                    },
                )
        })
        .reduce(
            || {
                (
                    Vec::with_capacity(points_len),
                    Vec::with_capacity(points_len),
                    Vec::with_capacity(points_len),
                )
            },
            |(mut xs1, mut ys1, mut ks1), (xs2, ys2, ks2)| {
                xs1.extend(xs2);
                ys1.extend(ys2);
                ks1.extend(ks2);
                (xs1, ys1, ks1)
            },
        );
    let mut add_sum: isize = 0;
    let add_map: Vec<isize> = ks
        .into_iter()
        .map(|k| {
            let distance = 1f64 + k;
            if distance < INSERT_THRESHOLD {
                if k == -1f64 {
                    add_sum -= 1;
                    return -1isize;
                }
                let n = ((INSERT_DENSITY / distance).ceil() as isize).clamp(0isize, INSERT_LIMIT);
                add_sum += n;
                n
            } else {
                0isize
            }
        })
        .collect();

    let len_added = ((add_map.len() as isize) + (add_sum * 2isize)) as usize;
    let mut xs_added: Vec<f64> = Vec::with_capacity(len_added);
    let mut ys_added: Vec<f64> = Vec::with_capacity(len_added);
    let mut last_point = Point { x: 0f64, y: 0f64 };
    let mut k_cache = 0;
    for (index, k) in add_map.iter().enumerate() {
        if *k == -1isize {
            continue;
        }
        let x = xs[index];
        let y = ys[index];
        let point = points[index];
        k_cache += *k;
        if (index != 0) && (k_cache != 0) {
            let delta_x = point.x - last_point.x;
            let delta_x = if delta_x > 180f64 {
                360f64 - delta_x
            } else if delta_x < -180f64 {
                360f64 + delta_x
            } else {
                delta_x
            };
            let t_x = delta_x / ((k_cache + 1) as f64);
            let t_y = (point.y - last_point.y) / ((k_cache + 1) as f64);
            let last_lon_delta = last_point.x - lon_base;
            for i in 1..k_cache {
                let (x_add, y_add, _) = latlon_to_azimuthal_equidistant(
                    last_point.y + (t_y * i as f64),
                    last_lon_delta + (t_x * i as f64),
                    lat_base,
                );
                xs_added.push(x_add);
                ys_added.push(y_add);
            }
        }
        last_point = point;
        xs_added.push(x);
        ys_added.push(y);
        k_cache = *k;
    }
    k_cache += add_map[0];
    let point = points[0];
    if k_cache != 0 {
        let t_x = (point.x - last_point.x) / ((k_cache + 1) as f64);
        let t_y = (point.y - last_point.y) / ((k_cache + 1) as f64);
        let last_lon_delta = last_point.x - lon_base;
        for i in 1..k_cache {
            let (x_add, y_add, _) = latlon_to_azimuthal_equidistant(
                last_point.y + (t_y * i as f64),
                last_lon_delta + (t_x * i as f64),
                lat_base,
            );
            xs_added.push(x_add);
            ys_added.push(y_add);
        }
    }
    (xs_added, ys_added)
}

#[repr(C)]
pub struct ReturnContent {
    pub status: bool,
    pub ptr: *const u8,
    pub len: usize,
}

impl ReturnContent {
    fn new(data: Vec<u8>, status: bool) -> Self {
        let data = ManuallyDrop::new(data.into_boxed_slice());
        ReturnContent {
            status: status,
            ptr: data.as_ptr(),
            len: data.len(),
        }
    }

    fn new_result(result: Result<Vec<u8>, shapefile::Error>) -> Self {
        match result {
            Ok(data) => Self::new(data, true),
            Err(e) => Self::new(e.to_string().into_bytes(), false),
        }
    }
}

#[repr(C)]
pub struct GenerateParameters {
    pub color_point: u32,
    pub color_multipoint: u32,
    pub color_line: u32,
    pub color_polygon: u32,
    pub width_point: i32,
    pub width_multipoint: i32,
    pub width_line: i32,
    // TODO add jump point
    pub fineness: u8,
    pub radius: u32,
    pub base_lat: f64,
    pub base_lon: f64,
}

#[unsafe(no_mangle)]
pub fn shapefile_generate_csupport(
    buffer_ptr: *const u8,
    buffer_len: usize,
    para: GenerateParameters,
) -> ReturnContent {
    ReturnContent::new_result(shapefile_generate(
        unsafe { std::slice::from_raw_parts(buffer_ptr, buffer_len) },
        para,
    ))
}

pub fn shapefile_generate(
    buffer: &[u8],
    para: GenerateParameters,
) -> Result<Vec<u8>, shapefile::Error> {
    let cursor = Cursor::new(buffer);
    let mut reader = match ShapeReader::new(cursor) {
        Ok(data) => data,
        Err(e) => return Err(e),
    };
    let size = (para.radius * 2) + 1;
    let mut img: RgbaImage = RgbaImage::from_pixel(size, size, Rgba([0, 0, 0, 0]));
    let linear_transformation_constant = para.radius as f64 / PI;
    let mut img_all = GrayImage::new(size, size);
    let radius = para.radius as i32;
    draw_filled_circle_mut(&mut img_all, (radius, radius), radius, image::Luma([255]));
    for shape in reader.iter_shapes() {
        shapefile_draw(
            &mut img,
            &mut img_all,
            &para,
            shape?,
            linear_transformation_constant,
        );
    }
    let mut result = Vec::new();
    let encoder = PngEncoder::new(&mut result);
    encoder
        .write_image(&img, size, size, ColorType::Rgba8.into())
        .unwrap();
    Ok(result)
}

pub fn latlon_line_draw(para: GenerateParameters, spacing: usize, fineness: usize) -> Vec<u8> {
    let size = (para.radius * 2) + 1;
    let mut img: RgbaImage = RgbaImage::from_pixel(size, size, Rgba([0, 0, 0, 0]));
    let linear_transformation_constant = para.radius as f64 / PI;
    let mut img_all = GrayImage::new(size, size);
    let radius = para.radius as i32;
    draw_filled_circle_mut(&mut img_all, (radius, radius), radius, image::Luma([255]));
    let mut shape_vec_latlon = Vec::<Vec<shapefile::Point>>::new();
    for x in -17..19 {
        let mut shape_vec_cache = Vec::<shapefile::Point>::new();
        for y in -899..900 {
            shape_vec_cache.push(shapefile::Point::new((x * 10) as f64, y as f64 / 10_f64));
        }
        shape_vec_cache.push(shapefile::Point::new((x * 10) as f64, -89.9_f64));
        shape_vec_latlon.push(shape_vec_cache);
    }
    for y in -8..9 {
        let mut shape_vec_cache = Vec::<shapefile::Point>::new();
        for x in -1799..1800 {
            shape_vec_cache.push(shapefile::Point::new(x as f64 / 10_f64, (y * 10) as f64));
        }
        shape_vec_cache.push(shapefile::Point::new(-179.9_f64, (y * 10) as f64));
        shape_vec_latlon.push(shape_vec_cache);
    }
    shapefile_draw(
        &mut img,
        &mut img_all,
        &para,
        shapefile::Shape::Polyline(shapefile::record::polyline::Polyline::with_parts(
            shape_vec_latlon,
        )),
        linear_transformation_constant,
    );
    let mut result = Vec::new();
    let encoder = PngEncoder::new(&mut result);
    encoder
        .write_image(&img, size, size, ColorType::Rgba8.into())
        .unwrap();
    result
}

fn shapefile_draw(
    img: &mut RgbaImage,
    img_all: &mut GrayImage,
    para: &GenerateParameters,
    shape: Shape,
    ltc: f64,
) {
    match shape {
        Shape::NullShape => return,
        Shape::Point(point) => shapefile_draw_point(
            img,
            para.color_point,
            para.width_point,
            para.base_lat,
            para.base_lon,
            ltc,
            point,
            para.radius as i32,
        ),
        Shape::PointM(point) => shapefile_draw_point(
            img,
            para.color_point,
            para.width_point,
            para.base_lat,
            para.base_lon,
            ltc,
            point,
            para.radius as i32,
        ),
        Shape::PointZ(point) => shapefile_draw_point(
            img,
            para.color_point,
            para.width_point,
            para.base_lat,
            para.base_lon,
            ltc,
            point,
            para.radius as i32,
        ),
        Shape::Multipoint(points) => shapefile_draw_multipoint(
            img,
            para.color_multipoint,
            para.width_multipoint,
            para.base_lat,
            para.base_lon,
            ltc,
            points,
            para.radius as i32,
        ),
        Shape::MultipointM(points) => shapefile_draw_multipoint(
            img,
            para.color_multipoint,
            para.width_multipoint,
            para.base_lat,
            para.base_lon,
            ltc,
            points,
            para.radius as i32,
        ),
        Shape::MultipointZ(points) => shapefile_draw_multipoint(
            img,
            para.color_multipoint,
            para.width_multipoint,
            para.base_lat,
            para.base_lon,
            ltc,
            points,
            para.radius as i32,
        ),
        Shape::Polyline(lines) => shapefile_draw_polyline(
            img,
            para.color_line,
            para.width_line,
            para.base_lat,
            para.base_lon,
            ltc,
            lines,
            para.radius as f32,
        ),
        Shape::PolylineM(lines) => shapefile_draw_polyline(
            img,
            para.color_line,
            para.width_line,
            para.base_lat,
            para.base_lon,
            ltc,
            lines,
            para.radius as f32,
        ),
        Shape::PolylineZ(lines) => shapefile_draw_polyline(
            img,
            para.color_line,
            para.width_line,
            para.base_lat,
            para.base_lon,
            ltc,
            lines,
            para.radius as f32,
        ),
        Shape::Polygon(rings) => shapefile_draw_polygon(
            img,
            img_all,
            para.color_polygon,
            para.base_lat,
            para.base_lon,
            ltc,
            rings,
            para.radius as i32,
        ),
        Shape::PolygonM(rings) => shapefile_draw_polygon(
            img,
            img_all,
            para.color_polygon,
            para.base_lat,
            para.base_lon,
            ltc,
            rings,
            para.radius as i32,
        ),
        Shape::PolygonZ(rings) => shapefile_draw_polygon(
            img,
            img_all,
            para.color_polygon,
            para.base_lat,
            para.base_lon,
            ltc,
            rings,
            para.radius as i32,
        ),
        _ => return,
    }
}

fn shapefile_draw_point<T: HasXY>(
    img: &mut RgbaImage,
    color: u32,
    width: i32,
    base_lat: f64,
    base_lon: f64,
    ltc: f64,
    point: T,
    r: i32,
) {
    let (x, y, _) = latlon_to_azimuthal_equidistant(point.y(), point.x() - base_lon, base_lat);
    let x = (x * ltc).round() as i32 + r;
    let y = (y * ltc).round() as i32 + r;
    draw_filled_circle_mut(img, (x, y), width, get_color_rgba(color));
}

fn shapefile_draw_multipoint<T: HasXY>(
    img: &mut RgbaImage,
    color_data: u32,
    width: i32,
    base_lat: f64,
    base_lon: f64,
    ltc: f64,
    points: GenericMultipoint<T>,
    r: i32,
) {
    let color = get_color_rgba(color_data);
    let (xs, ys) = latlon_to_azimuthal_equidistant_array(
        points.points().into_iter().map(Point::new).collect(),
        base_lat,
        base_lon,
    );
    for (x_origin, y_origin) in xs.iter().copied().zip(ys.iter().copied()) {
        let x = (x_origin * ltc).round() as i32 + r;
        let y = (y_origin * ltc).round() as i32 + r;
        draw_filled_circle_mut(img, (x, y), width, color);
    }
}

fn shapefile_draw_polyline<T: HasXY>(
    img: &mut RgbaImage,
    color_data: u32,
    width: i32,
    base_lat: f64,
    base_lon: f64,
    ltc: f64,
    lines: GenericPolyline<T>,
    r: f32,
) {
    let color = get_color_rgba(color_data);
    for part in lines.parts() {
        let (xs, ys) = latlon_to_azimuthal_equidistant_array(
            part.into_iter().map(Point::new).collect(),
            base_lat,
            base_lon,
        );
        let mut positions = xs.iter().copied().zip(ys.iter().copied());
        let (mut x_cache, mut y_cache) = match positions.next() {
            Some((x, y)) => ((x * ltc).round() as f32 + r, (y * ltc).round() as f32 + r),
            None => return,
        };
        for (x_origin, y_origin) in positions {
            let x = (x_origin * ltc).round() as f32 + r;
            let y = (y_origin * ltc).round() as f32 + r;
            draw_line_segment_mut(img, (x_cache, y_cache), (x, y), color);
            x_cache = x;
            y_cache = y;
        }
    }
}

fn shapefile_draw_polygon<T: HasXY>(
    img: &mut RgbaImage,
    img_all: &GrayImage,
    color_data: u32,
    base_lat: f64,
    base_lon: f64,
    ltc: f64,
    polygon: GenericPolygon<T>,
    r: i32,
) {
    let color = get_color_rgba(color_data);
    let size = ((r * 2) + 1) as u32;
    let mut img_out: GrayImage = GrayImage::new(size, size);
    let mut img_in: GrayImage = GrayImage::new(size, size);
    for ring in polygon.rings() {
        let (xs, ys) = latlon_to_azimuthal_equidistant_array(
            ring.points().into_iter().map(Point::new).collect(),
            base_lat,
            base_lon,
        );
        let is_cw = is_clockwise(&xs, &ys);
        let mut point_set = HashSet::<(i32, i32)>::new();
        let points: Vec<imageproc::point::Point<i32>> = xs
            .iter()
            .copied()
            .zip(ys.iter().copied())
            .map(|(x_origin, y_origin)| {
                let x = (x_origin * ltc).round() as i32 + r;
                let y = (y_origin * ltc).round() as i32 + r;
                let point = (x, y);
                if point_set.contains(&point) {
                    return None;
                }
                point_set.insert(point);
                Some(imageproc::point::Point::new(x, y))
            })
            .filter(|option| *option != None)
            .map(|some| match some {
                Some(s) => s,
                None => imageproc::point::Point::new(0i32, 0i32),
            })
            .collect();
        if points.len() < 3 {
            return;
        }
        let luma = image::Luma([255u8]);
        match (ring, is_cw) {
            (Outer(_), false) => draw_polygon_mut(&mut img_out, &points, luma),
            (Outer(_), true) => draw_merge(
                &draw_polygon(img_all, &points, image::Luma([0u8])),
                &mut img_out,
            ),
            (Inner(_), false) => draw_merge(
                &draw_polygon(img_all, &points, image::Luma([0u8])),
                &mut img_in,
            ),
            (Inner(_), true) => draw_polygon_mut(&mut img_in, &points, luma),
        }
        draw_with_mask_gray(&img_out, &img_in, img, color);
    }
}

pub fn is_clockwise(xs: &Vec<f64>, ys: &Vec<f64>) -> bool {
    let len = std::cmp::min(xs.len(), ys.len());
    if len < 3 {
        return false;
    }
    let mut sum = get_vector_product((xs[len - 1], ys[len - 1]), (xs[0], ys[0]));
    for i in 1..len {
        sum += get_vector_product((xs[i - 1], ys[i - 1]), (xs[i], ys[i]));
    }
    sum > 0_f64
}

fn get_vector_product((x_a, y_a): (f64, f64), (x_b, y_b): (f64, f64)) -> f64 {
    (x_b - x_a) * (y_a + y_b)
}

pub fn draw_merge(img_total: &GrayImage, img: &mut GrayImage) {
    for x in 0..img.width() {
        for y in 0..img.height() {
            let img_total_pixel = img_total.get_pixel(x, y);
            if img_total_pixel[0] > 0 {
                img.put_pixel(x, y, *img_total_pixel);
            }
        }
    }
}

pub fn draw_with_mask_gray(
    img_total: &GrayImage,
    mask: &GrayImage,
    img: &mut RgbaImage,
    color: Rgba<u8>,
) {
    for x in 0..img.width() {
        for y in 0..img.height() {
            if mask.get_pixel(x, y)[0] == 0 && img_total.get_pixel(x, y)[0] > 0 {
                img.put_pixel(x, y, color);
            }
        }
    }
}

pub fn draw_with_mask<T: GenericImage>(img_total: &T, mask: &GrayImage, img: &mut T) {
    for x in 0..img.width() {
        for y in 0..img.height() {
            if mask.get_pixel(x, y)[0] > 0 {
                img.put_pixel(x, y, img_total.get_pixel(x, y));
            }
        }
    }
}

fn get_color_rgba(color: u32) -> Rgba<u8> {
    Rgba([
        ((color >> 24) & 0xFF) as u8,
        ((color >> 16) & 0xFF) as u8,
        ((color >> 8) & 0xFF) as u8,
        (color & 0xFF) as u8,
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::{self, File};
    use std::io::Write;

    #[test]
    fn test_main() {
        let buffer: Vec<u8> = fs::read("./shapefile/ne_110m_land.shp").expect("IO Error");
        let latlon_parameter = GenerateParameters {
            color_point: 0x0000FFu32,
            color_multipoint: 0x0000FFu32,
            color_line: 0xFFFFFFFFu32,
            color_polygon: 0xF2E6C2FFu32,
            width_point: 3_i32,
            width_multipoint: 3_i32,
            width_line: 2_i32,
            fineness: 0x00u8,
            radius: 200u32,
            base_lat: 0090f64,
            base_lon: 0000f64,
        };
        let map_parameter = GenerateParameters {
            color_point: 0x0000FFu32,
            color_multipoint: 0x0000FFu32,
            color_line: 0xFFFFFFFFu32,
            color_polygon: 0xF2E6C2FFu32,
            width_point: 3_i32,
            width_multipoint: 3_i32,
            width_line: 2_i32,
            fineness: 0x00u8,
            radius: 200u32,
            base_lat: 0090f64,
            base_lon: 0000f64,
        };

        let mut out_file = File::create("./out/image.png").unwrap();
        let out_buffer = shapefile_generate(&buffer, map_parameter).expect("Shapefile Error");
        //let out_buffer = latlon_line_draw(latlon_parameter, 0, 0);
        let _ = out_file.write_all(&out_buffer);
    }
}
